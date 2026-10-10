
import re
import json
import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from io import StringIO

# ============================================================
# 1. CONFIGURATION
# ============================================================

SEED = 42
np.random.seed(SEED)

API_URL = "https://rest.uniprot.org/uniprotkb/search"

# Use a local Windows-compatible output directory.
OUTPUT_DIR = Path.cwd() / "protein_project"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

print("Configuration initialized successfully.")

# ============================================================
# 2. DOWNLOAD PILOT DATASET FROM UNIPROT
# ============================================================

params = {
    "query": "reviewed:true AND ec:*",
    "format": "tsv",
    "fields": (
        "accession,protein_name,organism_name,"
        "length,sequence,ec,annotation_score"
    ),
    "size": 500,
}

response = requests.get(
    API_URL,
    params=params,
    timeout=120
)
response.raise_for_status()

raw_df = pd.read_csv(
    StringIO(response.text),
    sep="\t"
)

print("\nRecords retrieved:", len(raw_df))
print("Available columns:", raw_df.columns.tolist())
print("\nFirst five records:")
print(raw_df.head())

if raw_df.empty:
    raise ValueError("UniProt returned no records.")

# Normalize expected column names for consistency.
raw_df.columns = raw_df.columns.str.strip()

required_columns = [
    "Entry",
    "Protein names",
    "Organism",
    "Sequence",
    "EC number",
]

missing_columns = [
    col for col in required_columns
    if col not in raw_df.columns
]

if missing_columns:
    raise ValueError(
        f"Expected columns are missing: {missing_columns}. "
        f"Actual columns: {raw_df.columns.tolist()}"
    )

# ============================================================
# 3. INITIAL DATA INSPECTION
# ============================================================

print("\nDataset shape:", raw_df.shape)

print("\nMissing values:")
print(raw_df.isna().sum())

raw_df["sequence_length_actual"] = (
    raw_df["Sequence"]
    .fillna("")
    .astype(str)
    .str.len()
)

print("\nSequence-length statistics:")
print(raw_df["sequence_length_actual"].describe())

# ============================================================
# 4. EXTRACT TOP-LEVEL EC CLASSES
# ============================================================

def extract_ec_classes(ec_text):
    """
    Extract distinct top-level classes from complete EC numbers.

    A complete four-component EC number is required.
    """
    if pd.isna(ec_text):
        return []

    matches = re.findall(
        r"(?<!\d)([1-7])\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)(?!\d)",
        str(ec_text)
    )

    return sorted(set(matches))


df = raw_df.copy()

df["top_level_classes"] = df["EC number"].apply(
    extract_ec_classes
)

df["n_top_level_classes"] = (
    df["top_level_classes"].str.len()
)

# Retain a single unambiguous top-level EC class.
df["label"] = df["top_level_classes"].apply(
    lambda classes: classes[0]
    if len(classes) == 1 else None
)

ambiguous_count = int(df["label"].isna().sum())

labelled_df = df.dropna(
    subset=["label"]
).copy()

labelled_df["label"] = (
    labelled_df["label"].astype(str)
)

print("\nRaw records:", len(df))
print("Missing or ambiguous class records:", ambiguous_count)
print("Unambiguously labelled records:", len(labelled_df))

print("\nLabel examples:")
print(
    labelled_df[
        [
            "Entry",
            "Protein names",
            "Organism",
            "EC number",
            "label",
        ]
    ].head()
)

# ============================================================
# 5. PROTEIN SEQUENCE QUALITY CONTROL
# ============================================================

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = set("BJZXUO")


def sequence_qc(sequence):
    """
    Calculate basic sequence-quality statistics.

    Missing sequences are represented as empty strings.
    Ambiguous residues are counted separately from invalid symbols.
    """
    if pd.isna(sequence):
        sequence = ""
    else:
        sequence = str(sequence).strip().upper()

    invalid = (
        set(sequence)
        - CANONICAL_AA
        - AMBIGUOUS_AA
    )

    return {
        "sequence": sequence,
        "sequence_length": len(sequence),
        "ambiguous_residues": sum(
            aa in AMBIGUOUS_AA for aa in sequence
        ),
        "invalid_residues": len(invalid),
        "valid_characters": len(invalid) == 0,
        "is_sequence_missing": len(sequence) == 0,
    }


# Avoid index alignment problems by constructing a new DataFrame.
qc = pd.DataFrame(
    labelled_df["Sequence"].apply(sequence_qc).tolist(),
    index=labelled_df.index,
)

# Assign each QC column explicitly.
for column in qc.columns:
    labelled_df[column] = qc[column]

# Explicitly enforce Boolean data types.
labelled_df["valid_characters"] = (
    labelled_df["valid_characters"].astype(bool)
)

labelled_df["is_sequence_missing"] = (
    labelled_df["is_sequence_missing"].astype(bool)
)

# Exact duplicate sequence flags.
labelled_df["is_exact_duplicate"] = (
    labelled_df["sequence"].duplicated(keep=False)
)

print("\nSequence quality report:")
print(
    "Missing or empty sequences:",
    int(labelled_df["is_sequence_missing"].sum())
)

print(
    "Sequences with invalid characters:",
    int((~labelled_df["valid_characters"]).sum())
)

print(
    "Sequences containing ambiguous residues:",
    int((labelled_df["ambiguous_residues"] > 0).sum())
)

print(
    "Records belonging to exact duplicate groups:",
    int(labelled_df["is_exact_duplicate"].sum())
)

print("\nSequence statistics:")
print(
    labelled_df[
        [
            "sequence_length",
            "ambiguous_residues",
            "invalid_residues",
        ]
    ].describe()
)

# ============================================================
# 6. PREPARE A QUALITY-CONTROLLED PILOT SUBSET
# ============================================================

# These are initial technical filters, not biological claims.
# The thresholds will be reviewed before final dataset construction.

MIN_SEQUENCE_LENGTH = 10
MAX_SEQUENCE_LENGTH = 2000

eligible_mask = (
    labelled_df["valid_characters"]
    & ~labelled_df["is_sequence_missing"]
    & labelled_df["sequence_length"].between(
        MIN_SEQUENCE_LENGTH,
        MAX_SEQUENCE_LENGTH
    )
)

# Exclude ambiguous residues from the initial model-ready subset.
# Preserve all records in labelled_df for auditing.
eligible_mask &= labelled_df["ambiguous_residues"].eq(0)

model_ready_df = labelled_df.loc[
    eligible_mask
].copy()

print("\nInitial filtering report:")
print("Labelled records:", len(labelled_df))
print("Records passing sequence filters:", len(model_ready_df))
print(
    "Records excluded by sequence filters:",
    len(labelled_df) - len(model_ready_df)
)

# ============================================================
# 7. EXACT DEDUPLICATION
# ============================================================

before_deduplication = len(model_ready_df)

# Keep the first occurrence of each exact sequence.
# Record the removed duplicate count.
model_ready_df = model_ready_df.drop_duplicates(
    subset=["sequence"],
    keep="first"
).copy()

removed_duplicates = (
    before_deduplication - len(model_ready_df)
)

print("\nDeduplication report:")
print("Records before deduplication:", before_deduplication)
print("Exact duplicate records removed:", removed_duplicates)
print("Unique sequences retained:", len(model_ready_df))

# ============================================================
# 8. CLASS DISTRIBUTION
# ============================================================

class_counts = (
    model_ready_df["label"]
    .value_counts()
    .sort_index()
)

print("\nClass distribution after filtering:")
print(
    class_counts.rename_axis("EC top-level class")
    .to_frame("protein_count")
)

if not class_counts.empty:
    ax = class_counts.plot(
        kind="bar",
        xlabel="EC top-level class",
        ylabel="Number of unique proteins",
        title="Pilot dataset: enzyme-class distribution",
        rot=0
    )

    plt.tight_layout()

    figure_path = OUTPUT_DIR / "class_distribution.png"
    plt.savefig(figure_path, dpi=200, bbox_inches="tight")
    plt.show()
    plt.close()

    print("Class-distribution plot saved:", figure_path)

# ============================================================
# 9. SAVE DATASETS AND METADATA
# ============================================================

# Preserve the unmodified API response.
raw_path = OUTPUT_DIR / "uniprot_raw_pilot.tsv"
raw_df.to_csv(raw_path, sep="\t", index=False)

# Preserve labelled records, including QC flags.
labelled_path = OUTPUT_DIR / "uniprot_labelled_pilot.csv"
labelled_df.to_csv(labelled_path, index=False)

# Save the initial model-ready subset.
model_ready_path = OUTPUT_DIR / "uniprot_model_ready_pilot.csv"
model_ready_df.to_csv(model_ready_path, index=False)

metadata = {
    "source": "UniProtKB REST API",
    "query": params["query"],
    "retrieval_date_utc": pd.Timestamp.now(
        tz="UTC"
    ).isoformat(),
    "random_seed": SEED,
    "maximum_requested_records": 500,
    "raw_record_count": int(len(raw_df)),
    "unambiguous_label_count": int(len(labelled_df)),
    "model_ready_record_count": int(len(model_ready_df)),
    "ambiguous_or_missing_class_count": ambiguous_count,
    "exact_duplicate_records_removed": int(removed_duplicates),
    "label_definition": (
        "Single unambiguous top-level class from a complete EC number"
    ),
    "sequence_length_range": [
        MIN_SEQUENCE_LENGTH,
        MAX_SEQUENCE_LENGTH,
    ],
    "ambiguous_residue_policy": (
        "Exclude from initial model-ready subset; preserve in audit data"
    ),
    "deduplication_policy": (
        "Exact amino-acid sequence; retain first occurrence"
    ),
    "limitations": [
        "Pilot size may be insufficient for reliable class-level evaluation",
        "Reviewed annotations do not guarantee experimental verification",
        "No sequence-similarity clustering or leakage-resistant split yet",
    ],
}

metadata_path = OUTPUT_DIR / "pilot_metadata.json"

with open(metadata_path, "w", encoding="utf-8") as file:
    json.dump(metadata, file, indent=2)

print("\nSaved artifacts:")
for path in sorted(OUTPUT_DIR.iterdir()):
    if path.is_file():
        print(f"{path.name}: {path.stat().st_size:,} bytes")

print("\nNotebook 01 data inspection completed.")

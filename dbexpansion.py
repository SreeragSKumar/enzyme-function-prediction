
import re
import json
import time
import hashlib
from pathlib import Path
from datetime import datetime, timezone
from io import StringIO

import requests
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# 1. CONFIGURATION
# ============================================================

SEED = 42
TARGET_TOTAL = 10_000
PAGE_SIZE = 500
MAX_PER_CLASS = (TARGET_TOTAL + 6) // 7

MIN_LENGTH = 10
MAX_LENGTH = 2000
ALLOW_AMBIGUOUS_RESIDUES = False

MAX_RETRIES = 5
REQUEST_DELAY_SECONDS = 0.4

BASE_URL = "https://rest.uniprot.org/uniprotkb/search"

OUTPUT_DIR = Path.cwd() / "protein_project"
CHECKPOINT_DIR = OUTPUT_DIR / "checkpoints"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
CHECKPOINT_DIR.mkdir(parents=True, exist_ok=True)

RAW_OUTPUT = OUTPUT_DIR / "uniprot_expanded_raw.tsv"
LABELLED_OUTPUT = OUTPUT_DIR / "uniprot_expanded_labelled.csv"
AUDIT_OUTPUT = OUTPUT_DIR / "uniprot_annotation_audit.csv"
READY_OUTPUT = OUTPUT_DIR / "uniprot_expanded_model_ready.csv"
METADATA_OUTPUT = OUTPUT_DIR / "expanded_dataset_metadata.json"
QUERY_LOG = OUTPUT_DIR / "uniprot_query_log.json"

SESSION = requests.Session()
SESSION.headers.update({
    "User-Agent": "EnzymeFunctionResearch/1.0"
})

FIELDS = [
    "accession",
    "id",
    "protein_name",
    "organism_name",
    "organism_id",
    "length",
    "sequence",
    "ec",
    "annotation_score",
    "protein_existence",
    "go",
    "go_f",
    "go_p",
    "go_c",
]

print("Output directory:", OUTPUT_DIR.resolve())


# ============================================================
# 2. HTTP REQUESTS WITH RETRIES
# ============================================================

def get_with_retries(url, params=None):
    for attempt in range(MAX_RETRIES):
        try:
            response = SESSION.get(
                url,
                params=params,
                timeout=(20, 120)
            )

            if response.status_code in (429, 500, 502, 503, 504):
                if attempt == MAX_RETRIES - 1:
                    response.raise_for_status()

                retry_after = response.headers.get("Retry-After")

                try:
                    wait_seconds = float(retry_after)
                except (TypeError, ValueError):
                    wait_seconds = min(2 ** attempt, 30)

                print(f"Retrying request in {wait_seconds} seconds...")
                time.sleep(wait_seconds)
                continue

            response.raise_for_status()
            return response

        except requests.RequestException as exc:
            if attempt == MAX_RETRIES - 1:
                raise

            wait_seconds = min(2 ** attempt, 30)
            print(
                f"Request failed: {exc}. "
                f"Retrying in {wait_seconds} seconds..."
            )
            time.sleep(wait_seconds)

    raise RuntimeError("Maximum request retries exceeded.")


# ============================================================
# 3. NORMALIZE UNIPROT COLUMN NAMES
# ============================================================

def normalize_columns(df):
    """
    Keep the actual UniProt headers but strip surrounding spaces.
    Add consistent internal names where necessary.
    """
    df = df.copy()
    df.columns = df.columns.astype(str).str.strip()

    aliases = {
        "Entry": "Entry",
        "Entry Name": "Entry Name",
        "Protein names": "Protein names",
        "Organism": "Organism",
        "Organism (ID)": "Organism (ID)",
        "Length": "Length",
        "Sequence": "Sequence",
        "EC number": "EC number",
        "Annotation": "Annotation",
        "Annotation score": "Annotation score",
        "Protein existence": "Protein existence",
        "Gene Ontology (GO)": "GO",
        "Gene Ontology (molecular function)":
            "Gene Ontology (molecular function)",
        "Gene Ontology (biological process)":
            "Gene Ontology (biological process)",
        "Gene Ontology (cellular component)":
            "Gene Ontology (cellular component)",
    }

    rename = {
        column: aliases[column]
        for column in df.columns
        if column in aliases
    }

    df = df.rename(columns=rename)

    # UniProt may not provide the optional annotation-score field.
    # Do not fabricate a score if the API omits it.
    return df


def get_next_url(response):
    """Read the next page URL from the HTTP Link header."""
    link_header = response.headers.get("Link", "")

    for item in link_header.split(","):
        if 'rel="next"' in item:
            match = re.search(r"<([^>]+)>", item)

            if match:
                return match.group(1)

    return None


# ============================================================
# 4. CHECKPOINTED, PAGINATED COLLECTION
# ============================================================

def load_json(path, default):
    if not path.exists():
        return default

    try:
        with open(path, "r", encoding="utf-8") as file:
            return json.load(file)
    except (json.JSONDecodeError, OSError):
        print(f"Warning: Could not read {path.name}.")
        return default


def save_json(path, data):
    temporary_path = path.with_suffix(path.suffix + ".tmp")

    with open(temporary_path, "w", encoding="utf-8") as file:
        json.dump(data, file, indent=2, ensure_ascii=False)

    temporary_path.replace(path)


def fetch_class(ec_class, max_records, query_log):
    """
    Retrieve one top-level EC class.

    A separate checkpoint and cursor state are maintained for each
    class. Progress is saved after every page.
    """
    class_dir = CHECKPOINT_DIR / f"ec_{ec_class}"
    class_dir.mkdir(parents=True, exist_ok=True)

    pages_file = class_dir / "pages.tsv"
    state_file = class_dir / "state.json"

    query = f"reviewed:true AND ec:{ec_class}.*"

    state = load_json(state_file, {
        "query": query,
        "next_url": None,
        "next_params": {
            "query": query,
            "format": "tsv",
            "fields": ",".join(FIELDS),
            "size": PAGE_SIZE,
        },
        "pages_retrieved": 0,
        "records_retrieved": 0,
        "complete": False,
    })

    if state.get("query") != query:
        raise RuntimeError(
            f"Checkpoint query mismatch for EC class {ec_class}."
        )

    if state.get("complete"):
        print(
            f"EC {ec_class}: using completed checkpoint "
            f"({state['records_retrieved']} records)."
        )

    else:
        next_url = state.get("next_url") or BASE_URL
        next_params = state.get("next_params")

        while True:
            if state["records_retrieved"] >= max_records:
                state["complete"] = True
                save_json(state_file, state)
                break

            response = get_with_retries(
                next_url,
                params=next_params
            )

            page_df = pd.read_csv(
                StringIO(response.text),
                sep="\t",
                dtype=str
            )

            page_df = normalize_columns(page_df)

            if page_df.empty:
                state["complete"] = True
                state["next_url"] = None
                state["next_params"] = None
                save_json(state_file, state)
                break

            if "Entry" not in page_df.columns:
                raise ValueError(
                    "UniProt response lacks the Entry column. "
                    f"Received: {page_df.columns.tolist()}"
                )

            # Add retrieval provenance.
            page_df["retrieval_query"] = query
            page_df["retrieval_ec_class"] = str(ec_class)
            page_df["retrieval_timestamp_utc"] = (
                datetime.now(timezone.utc).isoformat()
            )
            page_df["source_url"] = response.url

            remaining = max_records - state["records_retrieved"]
            page_df = page_df.head(remaining)

            # Write each downloaded page immediately.
            write_header = not pages_file.exists()

            page_df.to_csv(
                pages_file,
                sep="\t",
                index=False,
                mode="a",
                header=write_header
            )

            state["pages_retrieved"] += 1
            state["records_retrieved"] += len(page_df)

            next_page = get_next_url(response)

            state["next_url"] = next_page
            state["next_params"] = None

            if (
                next_page is None
                or state["records_retrieved"] >= max_records
            ):
                state["complete"] = True
                state["next_url"] = None
                state["next_params"] = None

            # Persist the cursor only after saving the page.
            save_json(state_file, state)

            print(
                f"EC {ec_class}: "
                f"page {state['pages_retrieved']}, "
                f"{state['records_retrieved']:,}/{max_records:,}"
            )

            if state["complete"]:
                break

            next_url = state["next_url"]
            time.sleep(REQUEST_DELAY_SECONDS)

    if pages_file.exists():
        class_df = pd.read_csv(
            pages_file,
            sep="\t",
            dtype=str
        )
        class_df = normalize_columns(class_df)
    else:
        class_df = pd.DataFrame()

    query_log[str(ec_class)] = {
        "query": query,
        "records_retrieved": int(len(class_df)),
        "pages_retrieved": int(state["pages_retrieved"]),
        "complete": bool(state["complete"]),
        "retrieved_at_utc": datetime.now(
            timezone.utc
        ).isoformat(),
    }

    save_json(QUERY_LOG, query_log)

    return class_df


# ============================================================
# 5. COLLECT ALL SEVEN CLASSES
# ============================================================

query_log = load_json(QUERY_LOG, {})
class_frames = []

for ec_class in range(1, 8):
    print(f"\nCollecting top-level EC class {ec_class}...")

    class_df = fetch_class(
        ec_class=ec_class,
        max_records=MAX_PER_CLASS,
        query_log=query_log
    )

    if not class_df.empty:
        class_frames.append(class_df)

if not class_frames:
    raise RuntimeError(
        "No records found. Check the internet connection and API response."
    )

raw_df = pd.concat(
    class_frames,
    ignore_index=True,
    sort=False
)

# UniProt can return the same accession for multiple class queries.
raw_df = raw_df.drop_duplicates(
    subset=["Entry"],
    keep="first"
).reset_index(drop=True)

raw_df.to_csv(
    RAW_OUTPUT,
    sep="\t",
    index=False
)

print("\nRaw unique accessions:", len(raw_df))
print("Raw columns:", raw_df.columns.tolist())


# ============================================================
# 6. STANDARDIZE REQUIRED COLUMNS
# ============================================================

# Keep the real downloaded data, and add missing optional fields as
# empty values so downstream auditing can proceed safely.
for column in [
    "Entry",
    "Entry Name",
    "Protein names",
    "Organism",
    "Organism (ID)",
    "Length",
    "Sequence",
    "EC number",
    "Annotation score",
    "Annotation",
    "Protein existence",
    "GO",
    "Gene Ontology (molecular function)",
    "Gene Ontology (biological process)",
    "Gene Ontology (cellular component)",
]:
    if column not in raw_df.columns:
        raw_df[column] = pd.NA

labelled_df = raw_df.copy()

# Resolve annotation score without assuming a single API header.
if labelled_df["Annotation score"].isna().all():
    if "Annotation" in labelled_df.columns:
        labelled_df["Annotation score"] = labelled_df["Annotation"]

labelled_df["Annotation score"] = pd.to_numeric(
    labelled_df["Annotation score"],
    errors="coerce"
)

labelled_df["Length"] = pd.to_numeric(
    labelled_df["Length"],
    errors="coerce"
)

print(
    "\nRecords with an available numeric annotation score:",
    int(labelled_df["Annotation score"].notna().sum())
)


# ============================================================
# 7. EXTRACT EC NUMBERS AND LABELS
# ============================================================

def extract_complete_ec_numbers(ec_text):
    if pd.isna(ec_text):
        return []

    # Only fully specified EC numbers are accepted as complete.
    pattern = (
        r"(?<![\d.])"
        r"([1-7]\.\d+\.\d+\.\d+)"
        r"(?![\d.])"
    )

    return sorted(set(re.findall(pattern, str(ec_text))))


def extract_top_level_classes(ec_text):
    return sorted({
        ec_number.split(".")[0]
        for ec_number in extract_complete_ec_numbers(ec_text)
    })


labelled_df["complete_ec_numbers"] = (
    labelled_df["EC number"].apply(extract_complete_ec_numbers)
)

labelled_df["top_level_classes"] = (
    labelled_df["EC number"].apply(extract_top_level_classes)
)

labelled_df["n_top_level_classes"] = (
    labelled_df["top_level_classes"].str.len()
)

labelled_df["label"] = labelled_df[
    "top_level_classes"
].apply(
    lambda classes: classes[0]
    if len(classes) == 1 else None
)

labelled_df["annotation_class_status"] = (
    labelled_df["top_level_classes"].apply(
        lambda classes:
            "missing" if len(classes) == 0
            else "single_class" if len(classes) == 1
            else "multiple_classes"
    )
)

labelled_df["retrieval_class_matches_annotation"] = (
    labelled_df.apply(
        lambda row: (
            str(row["retrieval_ec_class"])
            in row["top_level_classes"]
        ),
        axis=1
    )
)

print("\nAnnotation class status:")
print(
    labelled_df["annotation_class_status"]
    .value_counts(dropna=False)
)

print(
    "\nRetrieval/annotation class mismatches:",
    int(
        (~labelled_df["retrieval_class_matches_annotation"]).sum()
    )
)


# ============================================================
# 8. SEQUENCE QUALITY CONTROL
# ============================================================

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = set("BJZXUO")


def sequence_qc(sequence):
    if pd.isna(sequence):
        sequence = ""
    else:
        sequence = str(sequence).strip().upper()

    invalid = set(sequence) - CANONICAL_AA - AMBIGUOUS_AA

    return {
        "sequence_clean": sequence,
        "sequence_length_actual": len(sequence),
        "ambiguous_residue_count": sum(
            aa in AMBIGUOUS_AA for aa in sequence
        ),
        "invalid_residue_count": len(invalid),
        "has_invalid_characters": len(invalid) > 0,
        "is_empty_sequence": len(sequence) == 0,
    }


qc_df = pd.DataFrame(
    labelled_df["Sequence"].apply(sequence_qc).tolist(),
    index=labelled_df.index
)

labelled_df = pd.concat(
    [labelled_df, qc_df],
    axis=1
)

labelled_df["exact_sequence_hash"] = (
    labelled_df["sequence_clean"].apply(
        lambda sequence: hashlib.sha256(
            sequence.encode("utf-8")
        ).hexdigest()
    )
)

labelled_df["is_exact_duplicate"] = (
    labelled_df.duplicated(
        subset=["exact_sequence_hash"],
        keep=False
    )
)

labelled_df["length_matches_annotation"] = (
    labelled_df["Length"]
    == labelled_df["sequence_length_actual"]
)

print("\nSequence QC summary:")
print(
    "Empty sequences:",
    int(labelled_df["is_empty_sequence"].sum())
)
print(
    "Invalid-character sequences:",
    int(labelled_df["has_invalid_characters"].sum())
)
print(
    "Sequences with ambiguous residues:",
    int((labelled_df["ambiguous_residue_count"] > 0).sum())
)
print(
    "Records in exact-duplicate groups:",
    int(labelled_df["is_exact_duplicate"].sum())
)
print(
    "Reported/actual length mismatches:",
    int((~labelled_df["length_matches_annotation"]).sum())
)


# ============================================================
# 9. SAVE LABELLED DATA AND AUDIT TABLE
# ============================================================

labelled_df.to_csv(
    LABELLED_OUTPUT,
    index=False
)

audit_columns = [
    "Entry",
    "Entry Name",
    "Protein names",
    "Organism",
    "Organism (ID)",
    "EC number",
    "Annotation score",
    "Protein existence",
    "GO",
    "Gene Ontology (molecular function)",
    "Gene Ontology (biological process)",
    "Gene Ontology (cellular component)",
    "retrieval_query",
    "retrieval_ec_class",
    "retrieval_timestamp_utc",
    "source_url",
    "complete_ec_numbers",
    "top_level_classes",
    "annotation_class_status",
    "retrieval_class_matches_annotation",
    "sequence_length_actual",
    "ambiguous_residue_count",
    "invalid_residue_count",
    "has_invalid_characters",
    "is_exact_duplicate",
    "length_matches_annotation",
]

audit_columns = [
    column for column in audit_columns
    if column in labelled_df.columns
]

labelled_df[audit_columns].to_csv(
    AUDIT_OUTPUT,
    index=False
)


# ============================================================
# 10. DEFINE THE MODEL-READY DATASET
# ============================================================

valid_labels = [str(i) for i in range(1, 8)]

eligible = (
    labelled_df["label"].isin(valid_labels)
    & labelled_df["retrieval_class_matches_annotation"]
    & ~labelled_df["is_empty_sequence"]
    & ~labelled_df["has_invalid_characters"]
    & labelled_df["sequence_length_actual"].between(
        MIN_LENGTH,
        MAX_LENGTH
    )
)

if not ALLOW_AMBIGUOUS_RESIDUES:
    eligible &= labelled_df["ambiguous_residue_count"].eq(0)

# The initial seven-class experiment uses only proteins with
# one unambiguous top-level EC class.
eligible &= labelled_df["n_top_level_classes"].eq(1)

model_ready_df = labelled_df.loc[eligible].copy()

# Retain one occurrence of each exact sequence.
model_ready_df = model_ready_df.drop_duplicates(
    subset=["exact_sequence_hash"],
    keep="first"
).copy()

model_ready_df["label"] = model_ready_df["label"].astype(str)

print("\nModel-ready records:", len(model_ready_df))
print(
    "Unique model-ready sequences:",
    model_ready_df["sequence_clean"].nunique()
)

model_ready_df.to_csv(
    READY_OUTPUT,
    index=False
)


# ============================================================
# 11. CLASS COUNTS AND DISTRIBUTION REPORTS
# ============================================================

def class_count_table(data):
    counts = (
        data["label"]
        .dropna()
        .astype(str)
        .value_counts()
        .reindex(valid_labels, fill_value=0)
    )

    report = pd.DataFrame({
        "EC_class": valid_labels,
        "protein_count": counts.values,
    })

    total = int(report["protein_count"].sum())

    report["percentage"] = (
        report["protein_count"] / total * 100
        if total else 0.0
    )

    return report


raw_class_report = class_count_table(labelled_df)
ready_class_report = class_count_table(model_ready_df)

raw_class_report.to_csv(
    OUTPUT_DIR / "class_counts_before_filtering.csv",
    index=False
)

ready_class_report.to_csv(
    OUTPUT_DIR / "class_counts_after_filtering.csv",
    index=False
)

print("\nClass counts before filtering:")
print(raw_class_report.to_string(index=False))

print("\nClass counts after filtering:")
print(ready_class_report.to_string(index=False))


# ============================================================
# 12. PLOT CLASS DISTRIBUTIONS
# ============================================================

plots = [
    (
        "Expanded dataset: before filtering",
        raw_class_report,
        "class_distribution_before_filtering.png",
    ),
    (
        "Expanded dataset: model-ready",
        ready_class_report,
        "class_distribution_after_filtering.png",
    ),
]

for title, report, filename in plots:
    plt.figure(figsize=(9, 5))

    plt.bar(
        report["EC_class"],
        report["protein_count"]
    )

    plt.xlabel("Top-level EC class")
    plt.ylabel("Number of proteins")
    plt.title(title)

    plt.tight_layout()
    plt.savefig(
        OUTPUT_DIR / filename,
        dpi=200,
        bbox_inches="tight"
    )
    plt.close()


# ============================================================
# 13. SAVE METADATA
# ============================================================

metadata = {
    "source": "UniProtKB REST API",
    "retrieval_date_utc": datetime.now(
        timezone.utc
    ).isoformat(),
    "query_pattern": "reviewed:true AND ec:<class>.*",
    "fields_requested": FIELDS,
    "target_total_records": TARGET_TOTAL,
    "maximum_per_class": MAX_PER_CLASS,
    "raw_unique_accessions": int(len(raw_df)),
    "labelled_records": int(len(labelled_df)),
    "model_ready_records": int(len(model_ready_df)),
    "class_counts_before_filtering": {
        str(row["EC_class"]): int(row["protein_count"])
        for _, row in raw_class_report.iterrows()
    },
    "class_counts_after_filtering": {
        str(row["EC_class"]): int(row["protein_count"])
        for _, row in ready_class_report.iterrows()
    },
    "sequence_filter": {
        "minimum_length": MIN_LENGTH,
        "maximum_length": MAX_LENGTH,
        "ambiguous_residues_allowed": ALLOW_AMBIGUOUS_RESIDUES,
        "invalid_characters_allowed": False,
        "single_top_level_class_required": True,
        "exact_sequence_deduplication": True,
    },
    "limitations": [
        "Reviewed status alone does not guarantee experimental evidence.",
        "Annotation score is not a probability of functional correctness.",
        "Sequence similarity clustering has not yet been performed.",
        "The final class distribution is not guaranteed to be balanced.",
        "The dataset is not yet ready for leakage-resistant model evaluation.",
    ],
}

save_json(METADATA_OUTPUT, metadata)

print("\nExpansion and audit completed.")
print("Files saved in:", OUTPUT_DIR.resolve())

for path in sorted(OUTPUT_DIR.iterdir()):
    if path.is_file():
        print(f"{path.name}: {path.stat().st_size:,} bytes")

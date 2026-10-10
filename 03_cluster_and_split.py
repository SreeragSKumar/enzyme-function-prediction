
from pathlib import Path
import json
import subprocess
import hashlib
import sys

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold


# ============================================================
# CONFIGURATION
# ============================================================

PROJECT_DIR = Path(__file__).resolve().parent
DATA_PATH = (
    PROJECT_DIR
    / "protein_project"
    / "processed"
    / "uniprot_expanded_model_ready.csv"
)

MMSEQS = Path(r"C:\tools\mmseqs2\mmseqs\mmseqs.bat")

OUTPUT_DIR = Path(r"C:\protein_project_work\clustered")
FASTA_DIR = OUTPUT_DIR / "fasta"
MMSEQS_TMP = Path(r"C:\protein_project_work\mmseqs_tmp")

# Initial clustering parameters; review sensitivity later.
MIN_SEQ_ID = 0.30
COVERAGE = 0.80
COVERAGE_MODE = 0
THREADS = 4

# Target split proportions: train / validation / test
TARGET_PROPORTIONS = np.array([0.64, 0.16, 0.20])
RANDOM_SEEDS = range(100)

CANONICAL_AA = set("ACDEFGHIKLMNPQRSTVWY")


# ============================================================
# HELPERS
# ============================================================

def stop(message):
    raise RuntimeError(message)


def run_command(command):
    print("\nRunning:")
    print(subprocess.list2cmdline([str(x) for x in command]))
    subprocess.run(command, check=True)


def read_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def split_score(parts, overall_class_props):
    """
    Lower is better.
    Penalize split-size deviations and class-distribution deviations.
    """
    target_sizes = TARGET_PROPORTIONS
    actual_sizes = np.array([len(p) for p in parts], dtype=float)
    actual_sizes /= actual_sizes.sum()

    size_error = np.abs(actual_sizes - target_sizes).sum()

    class_error = 0.0
    for part in parts:
        proportions = (
            part["label"].value_counts(normalize=True)
            .reindex(range(1, 8), fill_value=0)
            .to_numpy()
        )
        class_error += np.abs(
            proportions - overall_class_props
        ).sum()

    # Missing a class in any split is strongly discouraged.
    missing_class_penalty = sum(
        (part["label"].nunique() < 7) * 10 for part in parts
    )

    return size_error + class_error + missing_class_penalty


# ============================================================
# 1. CHECK INPUTS
# ============================================================

if not DATA_PATH.exists():
    stop(f"Dataset not found: {DATA_PATH}")

if not MMSEQS.exists():
    stop(f"MMseqs2 executable not found: {MMSEQS}")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FASTA_DIR.mkdir(parents=True, exist_ok=True)
MMSEQS_TMP.parent.mkdir(parents=True, exist_ok=True)

print("Loading dataset:", DATA_PATH)
df = pd.read_csv(DATA_PATH)

required = {"Entry", "Sequence", "label"}
missing = required - set(df.columns)
if missing:
    stop(f"Missing required columns: {sorted(missing)}")

print("Original shape:", df.shape)

# Keep rows with valid identifiers, sequences and labels.
df = df.copy()
df["Entry"] = df["Entry"].astype(str).str.strip()
df["label"] = pd.to_numeric(df["label"], errors="coerce")
df["Sequence"] = (
    df["Sequence"]
    .fillna("")
    .astype(str)
    .str.upper()
    .str.replace(r"\s+", "", regex=True)
)

df = df[
    df["Entry"].ne("")
    & df["Sequence"].ne("")
    & df["label"].isin(range(1, 8))
].copy()

df["label"] = df["label"].astype(int)

# Prefer the prepared sequence column when it is complete.
if (
    "sequence_clean" in df.columns
    and df["sequence_clean"].notna().all()
    and df["sequence_clean"].astype(str).str.len().gt(0).all()
):
    df["_cluster_sequence"] = (
        df["sequence_clean"]
        .astype(str)
        .str.upper()
        .str.replace(r"\s+", "", regex=True)
    )
else:
    df["_cluster_sequence"] = df["Sequence"]

# Fail loudly rather than silently removing unexpected residues.
bad = ~df["_cluster_sequence"].map(
    lambda s: set(s).issubset(CANONICAL_AA)
)
if bad.any():
    examples = df.loc[bad, ["Entry", "_cluster_sequence"]].head(5)
    stop(
        f"{bad.sum()} sequences contain non-canonical residues. "
        "Inspect these records before clustering.\n"
        + examples.to_string(index=False)
    )

# Remove duplicate protein identifiers, then exact sequence duplicates.
before = len(df)
df = df.drop_duplicates(subset="Entry", keep="first").copy()
removed_duplicate_ids = before - len(df)

before = len(df)
df = df.drop_duplicates(subset="_cluster_sequence", keep="first").copy()
removed_exact_sequences = before - len(df)

if df.empty:
    stop("No usable records remain after input validation.")

df = df.reset_index(drop=True)
df["_fasta_id"] = [f"seq{i:06d}" for i in range(len(df))]

print("Records used:", len(df))
print("Duplicate identifiers removed:", removed_duplicate_ids)
print("Exact duplicate sequences removed:", removed_exact_sequences)
print("\nLabel counts:")
print(df["label"].value_counts().sort_index())

if set(df["label"].unique()) != set(range(1, 8)):
    stop("The usable dataset does not contain all seven labels.")


# ============================================================
# 2. EXPORT FASTA AND ID MAPPING
# ============================================================

fasta_path = FASTA_DIR / "proteins_for_mmseqs.fasta"
mapping_path = FASTA_DIR / "fasta_id_mapping.csv"

df[["_fasta_id", "Entry", "label", "_cluster_sequence"]].rename(
    columns={"_cluster_sequence": "clustering_sequence"}
).to_csv(mapping_path, index=False)


with fasta_path.open("w", encoding="utf-8", newline="\n") as handle:
    for fasta_id, sequence in zip(
        df["_fasta_id"],
        df["_cluster_sequence"],
    ):
        handle.write(f">{fasta_id}\n")
        handle.write(f"{sequence}\n")


print("\nFASTA written:", fasta_path)
print("FASTA sequences:", len(df))


# ============================================================
# 3. RUN MMSEQS2 CLUSTERING
# ============================================================

cluster_prefix = OUTPUT_DIR / "mmseqs_clusters"
cluster_tsv = Path(str(cluster_prefix) + "_cluster.tsv")

# Avoid accidentally reusing a stale result after changing input/settings.
if cluster_tsv.exists():
    cluster_tsv.unlink()

if MMSEQS_TMP.exists():
    import shutil
    shutil.rmtree(MMSEQS_TMP)
MMSEQS_TMP.mkdir(parents=True, exist_ok=True)

command = [
    str(MMSEQS),
    "easy-cluster",
    str(fasta_path),
    str(cluster_prefix),
    str(MMSEQS_TMP),
    "--min-seq-id", str(MIN_SEQ_ID),
    "-c", str(COVERAGE),
    "--cov-mode", str(COVERAGE_MODE),
    "--threads", str(THREADS),
]

run_command(command)

if not cluster_tsv.exists():
    stop(f"MMseqs2 did not produce the expected file: {cluster_tsv}")

clusters = pd.read_csv(
    cluster_tsv,
    sep="\t",
    header=None,
    names=["cluster_representative", "_fasta_id"],
)

if clusters.empty:
    stop("The MMseqs2 cluster TSV is empty.")

clusters = clusters.merge(
    df[["_fasta_id", "Entry", "label"]],
    on="_fasta_id",
    how="left",
    validate="many_to_one",
)

if clusters["Entry"].isna().any():
    stop("Some cluster FASTA IDs could not be mapped to dataset entries.")

if clusters["_fasta_id"].nunique() != len(df):
    stop(
        "The number of clustered sequences differs from the input count. "
        "Check MMseqs2 output and FASTA mapping."
    )

clusters = clusters.rename(
    columns={"cluster_representative": "cluster_id"}
)

cluster_map = clusters[["_fasta_id", "Entry", "cluster_id"]].copy()
cluster_map.to_csv(OUTPUT_DIR / "sequence_cluster_membership.csv", index=False)

df = df.merge(
    cluster_map[["Entry", "cluster_id"]],
    on="Entry",
    how="left",
    validate="one_to_one",
)

if df["cluster_id"].isna().any():
    stop("Some proteins have no cluster assignment.")

n_clusters = df["cluster_id"].nunique()
cluster_sizes = df.groupby("cluster_id").size()

print("\nMMseqs2 clustering completed.")
print("Sequences:", len(df))
print("Clusters:", n_clusters)
print("Largest cluster:", int(cluster_sizes.max()))
print("Singleton clusters:", int((cluster_sizes == 1).sum()))


# ============================================================
# 4. CREATE CLUSTER-GROUPED TRAIN / VALIDATION / TEST SPLITS
# ============================================================

# Each cluster is kept entirely within one split.
y = df["label"].to_numpy()
groups = df["cluster_id"].to_numpy()
overall_props = (
    df["label"].value_counts(normalize=True)
    .reindex(range(1, 8), fill_value=0)
    .to_numpy()
)

best = None
best_score = float("inf")

for seed in RANDOM_SEEDS:
    try:
        outer = StratifiedGroupKFold(
            n_splits=5, shuffle=True, random_state=seed
        )
        trainval_idx, test_idx = next(
            outer.split(np.zeros(len(df)), y, groups)
        )

        trainval = df.iloc[trainval_idx].copy().reset_index(drop=True)
        test = df.iloc[test_idx].copy().reset_index(drop=True)

        inner = StratifiedGroupKFold(
            n_splits=5, shuffle=True, random_state=seed + 1000
        )
        train_idx, val_idx = next(
            inner.split(
                np.zeros(len(trainval)),
                trainval["label"].to_numpy(),
                trainval["cluster_id"].to_numpy(),
            )
        )

        train = trainval.iloc[train_idx].copy()
        val = trainval.iloc[val_idx].copy()

        parts = [train, val, test]

        # Every class should be present in all splits.
        if any(set(p["label"].unique()) != set(range(1, 8)) for p in parts):
            continue

        score = split_score(parts, overall_props)
        if score < best_score:
            best_score = score
            best = (train.copy(), val.copy(), test.copy(), seed)

    except ValueError:
        continue

if best is None:
    stop(
        "Could not find splits containing all seven classes in all three "
        "sets. We may need to adjust clustering or use a different "
        "group-allocation method."
    )

train, val, test, chosen_seed = best

# Verify that no cluster crosses split boundaries.
train_clusters = set(train["cluster_id"])
val_clusters = set(val["cluster_id"])
test_clusters = set(test["cluster_id"])

assert train_clusters.isdisjoint(val_clusters)
assert train_clusters.isdisjoint(test_clusters)
assert val_clusters.isdisjoint(test_clusters)

assert len(train) + len(val) + len(test) == len(df)
assert set(train["Entry"]).isdisjoint(set(val["Entry"]))
assert set(train["Entry"]).isdisjoint(set(test["Entry"]))
assert set(val["Entry"]).isdisjoint(set(test["Entry"]))


# ============================================================
# 5. SAVE OUTPUTS AND REPORTS
# ============================================================

# Remove internal helper columns from the saved CSVs.
internal_cols = ["_cluster_sequence", "_fasta_id"]
for part in (df, train, val, test):
    part.drop(columns=internal_cols, errors="ignore", inplace=True)

df.to_csv(OUTPUT_DIR / "uniprot_clustered.csv", index=False)
train.to_csv(OUTPUT_DIR / "train.csv", index=False)
val.to_csv(OUTPUT_DIR / "validation.csv", index=False)
test.to_csv(OUTPUT_DIR / "test.csv", index=False)

count_table = pd.DataFrame(index=range(1, 8))
count_table.index.name = "label"

for name, part in [
    ("all", df),
    ("train", train),
    ("validation", val),
    ("test", test),
]:
    count_table[name] = (
        part["label"].value_counts()
        .reindex(range(1, 8), fill_value=0)
        .astype(int)
    )

count_table.to_csv(OUTPUT_DIR / "split_class_counts.csv")

split_summary = {
    "input_file": str(DATA_PATH),
    "mmseqs_executable": str(MMSEQS),
    "records_used": int(len(df)),
    "number_of_clusters": int(n_clusters),
    "min_sequence_identity": MIN_SEQ_ID,
    "coverage": COVERAGE,
    "coverage_mode": COVERAGE_MODE,
    "split_method": "StratifiedGroupKFold; clusters kept intact",
    "chosen_random_seed": int(chosen_seed),
    "split_score": float(best_score),
    "split_sizes": {
        "train": int(len(train)),
        "validation": int(len(val)),
        "test": int(len(test)),
    },
    "split_proportions": {
        "train": float(len(train) / len(df)),
        "validation": float(len(val) / len(df)),
        "test": float(len(test) / len(df)),
    },
    "cluster_counts": {
        "train": int(train["cluster_id"].nunique()),
        "validation": int(val["cluster_id"].nunique()),
        "test": int(test["cluster_id"].nunique()),
    },
    "class_counts": {
        name: {
            str(int(k)): int(v)
            for k, v in part["label"].value_counts().sort_index().items()
        }
        for name, part in [
            ("train", train),
            ("validation", val),
            ("test", test),
        ]
    },
    "cluster_overlap_checks": {
        "train_validation_overlap": len(train_clusters & val_clusters),
        "train_test_overlap": len(train_clusters & test_clusters),
        "validation_test_overlap": len(val_clusters & test_clusters),
    },
    "duplicate_identifiers_removed": int(removed_duplicate_ids),
    "exact_duplicate_sequences_removed": int(removed_exact_sequences),
}

with (OUTPUT_DIR / "split_metadata.json").open("w", encoding="utf-8") as f:
    json.dump(split_summary, f, indent=2)

print("\nCLASS COUNTS BY SPLIT")
print(count_table)

print("\nSPLIT SIZES")
print(f"Train:      {len(train):5d} ({len(train)/len(df):.1%})")
print(f"Validation: {len(val):5d} ({len(val)/len(df):.1%})")
print(f"Test:       {len(test):5d} ({len(test)/len(df):.1%})")

print("\nCLUSTER OVERLAP CHECKS")
print("Train/validation:", len(train_clusters & val_clusters))
print("Train/test:", len(train_clusters & test_clusters))
print("Validation/test:", len(val_clusters & test_clusters))

print("\nSUCCESS: clustering and grouped splitting completed.")
print("Outputs saved in:", OUTPUT_DIR)

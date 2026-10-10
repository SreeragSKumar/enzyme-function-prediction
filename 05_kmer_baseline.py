
"""
Experiment 2: Protein 3-mer features + Logistic Regression.

Fits the vectorizer and classifier on training data only.
Evaluates on validation data only.
The held-out test set is intentionally not loaded.
"""

import json
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    classification_report,
    confusion_matrix,
    ConfusionMatrixDisplay,
    f1_score,
)

# ---------------------------------------------------------
# 1. CONFIGURATION
# ---------------------------------------------------------

SEED = 42

DATA_DIR = Path(r"C:\protein_project_work\clustered")
OUTPUT_DIR = Path(r"C:\protein_project_work\kmer_baseline")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
KMER_SIZE = 3

TRAIN_PATH = DATA_DIR / "train.csv"
VALIDATION_PATH = DATA_DIR / "validation.csv"


# ---------------------------------------------------------
# 2. LOAD AND VALIDATE DATA
# ---------------------------------------------------------

def load_split(path):
    if not path.exists():
        raise FileNotFoundError(
            f"Dataset file not found: {path}\n"
            "Check that the cluster-and-split script "
            "has completed successfully."
        )

    df = pd.read_csv(path)

    if "label" not in df.columns:
        raise ValueError(
            f"'label' column missing from {path}. "
            f"Available columns: {df.columns.tolist()}"
        )

    # Prefer the cleaned sequence if available.
    if "sequence_clean" in df.columns:
        sequence_column = "sequence_clean"
    elif "Sequence" in df.columns:
        sequence_column = "Sequence"
    elif "sequence" in df.columns:
        sequence_column = "sequence"
    else:
        raise ValueError(
            f"No sequence column found in {path}. "
            f"Available columns: {df.columns.tolist()}"
        )

    sequences = (
        df[sequence_column]
        .fillna("")
        .astype(str)
        .str.upper()
        .str.strip()
    )

    labels = pd.to_numeric(df["label"], errors="raise").astype(int)

    # The current experiment uses the 20 standard amino acids.
    invalid_rows = [
        i for i, seq in enumerate(sequences)
        if not seq or not set(seq).issubset(AMINO_ACIDS)
    ]

    if invalid_rows:
        raise ValueError(
            f"{path.name} contains {len(invalid_rows)} empty "
            "or non-standard sequences. Inspect sequence "
            "cleaning before training; do not silently discard them."
        )

    print(f"\nLoaded: {path.name}")
    print(f"Records: {len(df)}")
    print(f"Sequence column: {sequence_column}")
    print(f"Labels: {sorted(labels.unique().tolist())}")

    return sequences.tolist(), labels.to_numpy()


X_train_seq, y_train = load_split(TRAIN_PATH)
X_val_seq, y_val = load_split(VALIDATION_PATH)

if set(np.unique(y_train)) != set(np.unique(y_val)):
    raise ValueError(
        "Training and validation label sets differ. "
        "Inspect the split files before proceeding."
    )

print("\nTest set will NOT be loaded or evaluated.")


# ---------------------------------------------------------
# 3. EXTRACT 3-MER FEATURES
# ---------------------------------------------------------

print("\nFitting the 3-mer TF-IDF vectorizer...")

vectorizer = TfidfVectorizer(
    analyzer="char",
    ngram_range=(KMER_SIZE, KMER_SIZE),
    lowercase=False,
    use_idf=True,
    sublinear_tf=True,
    norm="l2",
    dtype=np.float32,
)

# Fit the vocabulary and IDF weights on training sequences only.
X_train = vectorizer.fit_transform(X_train_seq)

# Transform validation data using the existing training vocabulary.
X_val = vectorizer.transform(X_val_seq)

print(f"Training matrix shape: {X_train.shape}")
print(f"Validation matrix shape: {X_val.shape}")
print(f"Number of learned 3-mer features: {len(vectorizer.vocabulary_)}")


# ---------------------------------------------------------
# 4. TRAIN LOGISTIC REGRESSION
# ---------------------------------------------------------

print("\nTraining Logistic Regression...")

model = LogisticRegression(
    max_iter=3000,
    class_weight="balanced",
    random_state=SEED,
    solver="lbfgs",
)

model.fit(X_train, y_train)

print("Training completed.")


# ---------------------------------------------------------
# 5. VALIDATION EVALUATION
# ---------------------------------------------------------

print("\nEvaluating on the validation set...")

y_pred = model.predict(X_val)

accuracy = accuracy_score(y_val, y_pred)
balanced_accuracy = balanced_accuracy_score(y_val, y_pred)
macro_f1 = f1_score(y_val, y_pred, average="macro", zero_division=0)
weighted_f1 = f1_score(
    y_val, y_pred, average="weighted", zero_division=0
)

metrics = {
    "experiment": "3-mer TF-IDF + Logistic Regression",
    "kmer_size": KMER_SIZE,
    "training_records": len(y_train),
    "validation_records": len(y_val),
    "number_of_features": int(X_train.shape[1]),
    "accuracy": float(accuracy),
    "balanced_accuracy": float(balanced_accuracy),
    "macro_f1": float(macro_f1),
    "weighted_f1": float(weighted_f1),
    "test_set_used": False,
    "random_seed": SEED,
}

report = classification_report(
    y_val,
    y_pred,
    labels=sorted(np.unique(y_train)),
    output_dict=True,
    zero_division=0,
)

report_df = pd.DataFrame(report).transpose()

labels = sorted(np.unique(y_train))
cm = confusion_matrix(y_val, y_pred, labels=labels)
cm_df = pd.DataFrame(
    cm,
    index=[f"True_{label}" for label in labels],
    columns=[f"Predicted_{label}" for label in labels],
)


# ---------------------------------------------------------
# 6. SAVE RESULTS
# ---------------------------------------------------------

with open(OUTPUT_DIR / "metrics.json", "w", encoding="utf-8") as f:
    json.dump(metrics, f, indent=4)

report_df.to_csv(OUTPUT_DIR / "classification_report.csv")
cm_df.to_csv(OUTPUT_DIR / "confusion_matrix.csv")

fig, ax = plt.subplots(figsize=(9, 7))

display = ConfusionMatrixDisplay(
    confusion_matrix=cm,
    display_labels=labels,
)

display.plot(
    ax=ax,
    cmap="Blues",
    values_format="d",
    colorbar=False,
)

ax.set_title("3-mer TF-IDF + Logistic Regression\nValidation Confusion Matrix")
fig.tight_layout()
fig.savefig(
    OUTPUT_DIR / "confusion_matrix.png",
    dpi=200,
    bbox_inches="tight",
)
plt.close(fig)

joblib.dump(model, OUTPUT_DIR / "logistic_regression.joblib")
joblib.dump(vectorizer, OUTPUT_DIR / "tfidf_vectorizer.joblib")

print("\n" + "=" * 55)
print("3-MER BASELINE — VALIDATION RESULTS")
print("=" * 55)

for key, value in metrics.items():
    print(f"{key}: {value}")

print("\nClassification report:")
print(report_df.round(4).to_string())

print(f"\nAll outputs saved to:\n{OUTPUT_DIR}")
print("The held-out test set remains untouched.")

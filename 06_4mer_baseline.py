
"""
Experiment 3: 4-mer TF-IDF + Logistic Regression.

Uses the existing cluster-based training and validation splits.
Fits TF-IDF on training sequences only.
Does not load or evaluate the held-out test set.
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
KMER_SIZE = 4

DATA_DIR = Path(r"C:\protein_project_work\clustered")
OUTPUT_DIR = Path(r"C:\protein_project_work\4mer_baseline")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TRAIN_PATH = DATA_DIR / "train.csv"
VALIDATION_PATH = DATA_DIR / "validation.csv"

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------
# 2. LOAD AND VALIDATE DATA
# ---------------------------------------------------------

def load_split(path):
    if not path.exists():
        raise FileNotFoundError(f"Dataset file not found: {path}")

    df = pd.read_csv(path)

    if "label" not in df.columns:
        raise ValueError(
            f"'label' column missing from {path}. "
            f"Available columns: {df.columns.tolist()}"
        )

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

    invalid = [
        i for i, seq in enumerate(sequences)
        if len(seq) < KMER_SIZE
        or not set(seq).issubset(AMINO_ACIDS)
    ]

    if invalid:
        raise ValueError(
            f"{path.name} contains {len(invalid)} empty, short, "
            "or non-standard sequences. Check the data before training."
        )

    print(f"\nLoaded: {path.name}")
    print(f"Records: {len(df)}")
    print(f"Sequence column: {sequence_column}")
    print(f"Labels: {sorted(labels.unique().tolist())}")

    return sequences.tolist(), labels.to_numpy()


X_train_seq, y_train = load_split(TRAIN_PATH)
X_val_seq, y_val = load_split(VALIDATION_PATH)

if set(np.unique(y_train)) != set(np.unique(y_val)):
    raise ValueError("Training and validation label sets differ.")

print("\nHeld-out test set will not be loaded.")


# ---------------------------------------------------------
# 3. EXTRACT 4-MER TF-IDF FEATURES
# ---------------------------------------------------------

print("\nFitting 4-mer TF-IDF vectorizer...")

vectorizer = TfidfVectorizer(
    analyzer="char",
    ngram_range=(KMER_SIZE, KMER_SIZE),
    lowercase=False,
    use_idf=True,
    sublinear_tf=True,
    norm="l2",
    dtype=np.float32,
)

# Learn vocabulary and IDF weights from training data only.
X_train = vectorizer.fit_transform(X_train_seq)

# Apply the fitted vectorizer to validation data.
X_val = vectorizer.transform(X_val_seq)

print(f"Training matrix shape: {X_train.shape}")
print(f"Validation matrix shape: {X_val.shape}")
print(f"Learned 4-mer features: {len(vectorizer.vocabulary_)}")


# ---------------------------------------------------------
# 4. TRAIN CLASSIFIER
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

print("\nEvaluating on validation data...")

y_pred = model.predict(X_val)

metrics = {
    "experiment": "4-mer TF-IDF + Logistic Regression",
    "kmer_size": KMER_SIZE,
    "training_records": int(len(y_train)),
    "validation_records": int(len(y_val)),
    "number_of_features": int(X_train.shape[1]),
    "accuracy": float(accuracy_score(y_val, y_pred)),
    "balanced_accuracy": float(
        balanced_accuracy_score(y_val, y_pred)
    ),
    "macro_f1": float(
        f1_score(y_val, y_pred, average="macro", zero_division=0)
    ),
    "weighted_f1": float(
        f1_score(y_val, y_pred, average="weighted", zero_division=0)
    ),
    "test_set_used": False,
    "random_seed": SEED,
}

labels = sorted(np.unique(y_train))

report = classification_report(
    y_val,
    y_pred,
    labels=labels,
    output_dict=True,
    zero_division=0,
)

report_df = pd.DataFrame(report).transpose()

cm = confusion_matrix(y_val, y_pred, labels=labels)

cm_df = pd.DataFrame(
    cm,
    index=[f"True_{label}" for label in labels],
    columns=[f"Predicted_{label}" for label in labels],
)


# ---------------------------------------------------------
# 6. SAVE OUTPUTS
# ---------------------------------------------------------

with open(
    OUTPUT_DIR / "metrics.json", "w", encoding="utf-8"
) as f:
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

ax.set_title("4-mer TF-IDF + Logistic Regression\nValidation Confusion Matrix")
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
print("4-MER BASELINE — VALIDATION RESULTS")
print("=" * 55)

for key, value in metrics.items():
    print(f"{key}: {value}")

print("\nClassification report:")
print(report_df.round(4).to_string())

print(f"\nOutputs saved to: {OUTPUT_DIR}")
print("The held-out test set remains untouched.")

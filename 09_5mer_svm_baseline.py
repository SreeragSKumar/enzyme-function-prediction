
import json
import time
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.svm import LinearSVC
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    f1_score,
    classification_report,
    confusion_matrix,
    ConfusionMatrixDisplay,
)

# --------------------------------------------------
# 1. CONFIGURATION
# --------------------------------------------------

DATA_DIR = Path(r"C:\protein_project_work\clustered")
OUTPUT_DIR = Path(r"C:\protein_project_work\5mer_svm")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

KMER_SIZE = 5
MAX_FEATURES = 100_000
C_VALUE = 1.0
RANDOM_SEED = 42

# --------------------------------------------------
# 2. LOAD EXISTING CLUSTER-BASED SPLITS
# --------------------------------------------------

train_df = pd.read_csv(DATA_DIR / "train.csv")
val_df = pd.read_csv(DATA_DIR / "validation.csv")

required_columns = {"Sequence", "label"}
for split_name, frame in [
    ("training", train_df),
    ("validation", val_df),
]:
    missing = required_columns - set(frame.columns)
    if missing:
        raise ValueError(
            f"{split_name} CSV is missing columns: {sorted(missing)}"
        )

def clean_sequence(sequence):
    """Keep only standard amino-acid letters."""
    sequence = str(sequence).upper()
    return "".join(
        aa for aa in sequence
        if aa in "ACDEFGHIKLMNPQRSTVWY"
    )

X_train_text = train_df["Sequence"].map(clean_sequence)
X_val_text = val_df["Sequence"].map(clean_sequence)

y_train = train_df["label"].astype(int)
y_val = val_df["label"].astype(int)

if (X_train_text.str.len() == 0).any():
    raise ValueError("At least one training sequence is empty after cleaning.")

if (X_val_text.str.len() == 0).any():
    raise ValueError("At least one validation sequence is empty after cleaning.")

print("=" * 65)
print("5-MER TF-IDF + LINEAR SVM")
print("=" * 65)
print(f"Training records:   {len(train_df)}")
print(f"Validation records: {len(val_df)}")
print(f"K-mer size:         {KMER_SIZE}")
print(f"Maximum features:   {MAX_FEATURES}")
print(f"SVM C:              {C_VALUE}")
print("Test set:           NOT USED")

# --------------------------------------------------
# 3. CONVERT SEQUENCES INTO 5-MER FEATURES
# --------------------------------------------------

# A k-mer is a consecutive substring of k amino acids.
# For example, ACDEFG produces ACDEF and CDEFG for k=5.
# TF-IDF weights k-mers according to their usefulness
# across the training sequences.

vectorizer = TfidfVectorizer(
    analyzer="char",
    ngram_range=(KMER_SIZE, KMER_SIZE),
    lowercase=False,
    max_features=MAX_FEATURES,
    sublinear_tf=True,
    dtype="float32",
)

start = time.time()

X_train = vectorizer.fit_transform(X_train_text)
X_val = vectorizer.transform(X_val_text)

print(f"\nTraining feature matrix: {X_train.shape}")
print(f"Validation feature matrix: {X_val.shape}")
print(f"Feature extraction time: {time.time() - start:.2f} seconds")

# --------------------------------------------------
# 4. TRAIN THE CLASSIFIER
# --------------------------------------------------

model = LinearSVC(
    C=C_VALUE,
    class_weight="balanced",
    random_state=RANDOM_SEED,
    max_iter=10000,
)

start = time.time()
model.fit(X_train, y_train)
training_seconds = time.time() - start

# --------------------------------------------------
# 5. VALIDATION EVALUATION
# --------------------------------------------------

predictions = model.predict(X_val)

accuracy = accuracy_score(y_val, predictions)
balanced_accuracy = balanced_accuracy_score(y_val, predictions)
macro_f1 = f1_score(y_val, predictions, average="macro")
weighted_f1 = f1_score(y_val, predictions, average="weighted")

report = classification_report(
    y_val,
    predictions,
    labels=[1, 2, 3, 4, 5, 6, 7],
    output_dict=True,
    zero_division=0,
)

print("\nVALIDATION RESULTS")
print("=" * 65)
print(f"Accuracy:          {accuracy:.4f}")
print(f"Balanced accuracy: {balanced_accuracy:.4f}")
print(f"Macro F1:          {macro_f1:.4f}")
print(f"Weighted F1:       {weighted_f1:.4f}")
print(f"Training time:     {training_seconds:.2f} seconds")
print("\nClassification report:")
print(
    classification_report(
        y_val,
        predictions,
        labels=[1, 2, 3, 4, 5, 6, 7],
        zero_division=0,
    )
)

# --------------------------------------------------
# 6. SAVE METRICS AND CLASSIFICATION REPORT
# --------------------------------------------------

metrics = {
    "experiment": "5-mer TF-IDF + Linear SVM",
    "kmer_size": KMER_SIZE,
    "max_features": MAX_FEATURES,
    "actual_features": int(X_train.shape[1]),
    "C": C_VALUE,
    "class_weight": "balanced",
    "training_records": int(len(train_df)),
    "validation_records": int(len(val_df)),
    "accuracy": float(accuracy),
    "balanced_accuracy": float(balanced_accuracy),
    "macro_f1": float(macro_f1),
    "weighted_f1": float(weighted_f1),
    "training_seconds": float(training_seconds),
    "test_set_used": False,
    "random_seed": RANDOM_SEED,
}

with open(OUTPUT_DIR / "metrics.json", "w", encoding="utf-8") as f:
    json.dump(metrics, f, indent=4)

pd.DataFrame(report).transpose().to_csv(
    OUTPUT_DIR / "classification_report.csv",
    index=True,
)

cm = confusion_matrix(
    y_val,
    predictions,
    labels=[1, 2, 3, 4, 5, 6, 7],
)

pd.DataFrame(
    cm,
    index=[f"True_{i}" for i in range(1, 8)],
    columns=[f"Predicted_{i}" for i in range(1, 8)],
).to_csv(OUTPUT_DIR / "confusion_matrix.csv")

# --------------------------------------------------
# 7. SAVE CONFUSION MATRIX FIGURE
# --------------------------------------------------

fig, ax = plt.subplots(figsize=(9, 8))

display = ConfusionMatrixDisplay(
    confusion_matrix=cm,
    display_labels=[1, 2, 3, 4, 5, 6, 7],
)

display.plot(
    ax=ax,
    cmap="Blues",
    values_format="d",
    colorbar=False,
)

ax.set_title("5-mer TF-IDF + Linear SVM\nValidation Confusion Matrix")
fig.tight_layout()
fig.savefig(
    OUTPUT_DIR / "confusion_matrix.png",
    dpi=200,
    bbox_inches="tight",
)
plt.close(fig)

# --------------------------------------------------
# 8. SAVE MODEL AND VECTORIZER
# --------------------------------------------------

joblib.dump(model, OUTPUT_DIR / "model.joblib")
joblib.dump(vectorizer, OUTPUT_DIR / "tfidf_vectorizer.joblib")

print(f"\nAll outputs saved to: {OUTPUT_DIR}")
print("The held-out test set remains untouched.")

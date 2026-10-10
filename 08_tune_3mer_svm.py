
"""
Experiment 5: Tune the regularization parameter C for a
3-mer TF-IDF + Linear SVM enzyme classification model.

Uses the existing cluster-based train/validation splits.
Fits TF-IDF on training sequences only.
Selects C using validation macro F1.
Does not load or evaluate the held-out test set.
"""

import json
import time
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.svm import LinearSVC
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
KMER_SIZE = 3

C_VALUES = [0.1, 0.5, 1.0, 2.0, 5.0]

DATA_DIR = Path(r"C:\protein_project_work\clustered")
OUTPUT_DIR = Path(r"C:\protein_project_work\svm_tuning")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

TRAIN_PATH = DATA_DIR / "train.csv"
VALIDATION_PATH = DATA_DIR / "validation.csv"

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------
# 2. LOAD DATA
# ---------------------------------------------------------

def load_split(path):
    if not path.exists():
        raise FileNotFoundError(f"Dataset file not found: {path}")

    df = pd.read_csv(path)

    if "label" not in df.columns:
        raise ValueError(f"'label' column missing from {path}")

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
            "or non-standard sequences."
        )

    return sequences.tolist(), labels.to_numpy()


print("Loading training and validation data...")

X_train_seq, y_train = load_split(TRAIN_PATH)
X_val_seq, y_val = load_split(VALIDATION_PATH)

if set(np.unique(y_train)) != set(np.unique(y_val)):
    raise ValueError("Training and validation label sets differ.")

print(f"Training records: {len(y_train)}")
print(f"Validation records: {len(y_val)}")
print("Held-out test set will not be loaded.")


# ---------------------------------------------------------
# 3. FIT 3-MER TF-IDF ON TRAINING DATA ONLY
# ---------------------------------------------------------

print("\nFitting 3-mer TF-IDF...")

vectorizer = TfidfVectorizer(
    analyzer="char",
    ngram_range=(KMER_SIZE, KMER_SIZE),
    lowercase=False,
    use_idf=True,
    sublinear_tf=True,
    norm="l2",
    dtype=np.float32,
)

X_train = vectorizer.fit_transform(X_train_seq)
X_val = vectorizer.transform(X_val_seq)

print(f"Training feature matrix: {X_train.shape}")
print(f"Validation feature matrix: {X_val.shape}")


# ---------------------------------------------------------
# 4. TRAIN AND COMPARE C VALUES
# ---------------------------------------------------------

results = []
trained_models = {}

for C in C_VALUES:
    print(f"\nTraining Linear SVM with C={C}...")

    start_time = time.time()

    model = LinearSVC(
        C=C,
        class_weight="balanced",
        random_state=SEED,
        max_iter=10000,
    )

    model.fit(X_train, y_train)
    y_pred = model.predict(X_val)

    elapsed = time.time() - start_time

    result = {
        "C": C,
        "accuracy": accuracy_score(y_val, y_pred),
        "balanced_accuracy": balanced_accuracy_score(y_val, y_pred),
        "macro_f1": f1_score(
            y_val, y_pred, average="macro", zero_division=0
        ),
        "weighted_f1": f1_score(
            y_val, y_pred, average="weighted", zero_division=0
        ),
        "training_seconds": elapsed,
        "iterations": int(np.max(model.n_iter_)),
    }

    results.append(result)
    trained_models[C] = model

    print(
        f"Accuracy={result['accuracy']:.4f} | "
        f"Balanced accuracy={result['balanced_accuracy']:.4f} | "
        f"Macro F1={result['macro_f1']:.4f} | "
        f"Weighted F1={result['weighted_f1']:.4f} | "
        f"Time={elapsed:.1f}s"
    )


# ---------------------------------------------------------
# 5. SELECT BEST C USING VALIDATION MACRO F1
# ---------------------------------------------------------

results_df = pd.DataFrame(results)

# In a tie, prefer the smaller C as a simple, deterministic rule.
results_df = results_df.sort_values(
    by=["macro_f1", "C"],
    ascending=[False, True],
).reset_index(drop=True)

best_C = float(results_df.iloc[0]["C"])
best_model = trained_models[best_C]

print("\n" + "=" * 60)
print("SVM HYPERPARAMETER TUNING RESULTS")
print("=" * 60)
print(results_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
print(f"\nBest C based on validation macro F1: {best_C}")


# ---------------------------------------------------------
# 6. FINAL VALIDATION REPORT FOR SELECTED MODEL
# ---------------------------------------------------------

y_pred = best_model.predict(X_val)
labels = sorted(np.unique(y_train))

metrics = {
    "experiment": "3-mer TF-IDF + tuned Linear SVM",
    "kmer_size": KMER_SIZE,
    "selected_C": best_C,
    "selection_metric": "validation_macro_f1",
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
# 7. SAVE RESULTS AND MODEL
# ---------------------------------------------------------

results_df.to_csv(
    OUTPUT_DIR / "c_value_comparison.csv",
    index=False,
)

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

ax.set_title(
    f"3-mer TF-IDF + Linear SVM (C={best_C})\n"
    "Validation Confusion Matrix"
)

fig.tight_layout()
fig.savefig(
    OUTPUT_DIR / "confusion_matrix.png",
    dpi=200,
    bbox_inches="tight",
)
plt.close(fig)

joblib.dump(best_model, OUTPUT_DIR / "best_linear_svm.joblib")
joblib.dump(vectorizer, OUTPUT_DIR / "tfidf_vectorizer.joblib")

print("\nSelected model validation metrics:")
for key, value in metrics.items():
    print(f"{key}: {value}")

print("\nClassification report:")
print(report_df.round(4).to_string())

print(f"\nAll outputs saved to: {OUTPUT_DIR}")
print("The held-out test set remains untouched.")

"""
Stage 2: Amino-acid-composition baseline for enzyme EC-class prediction.

Inputs:
  C:\\protein_project_work\\clustered\\train.csv
  C:\\protein_project_work\\clustered\\validation.csv

This first run deliberately does not load or evaluate test.csv.
It trains on train.csv and evaluates on validation.csv only.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import sklearn
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (
        accuracy_score,
        balanced_accuracy_score,
        classification_report,
        confusion_matrix,
        ConfusionMatrixDisplay,
        f1_score,
    )
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler
except ImportError as exc:
    print(f"\nMissing Python package: {exc}")
    print("Install dependencies with:")
    print("python -m pip install pandas numpy matplotlib scikit-learn")
    sys.exit(1)


DATA_DIR = Path(r"C:\protein_project_work\clustered")
OUTPUT_DIR = Path(r"C:\protein_project_work\baseline")

TRAIN_FILE = DATA_DIR / "train.csv"
VALIDATION_FILE = DATA_DIR / "validation.csv"

# The 20 standard amino acids. Other symbols are ignored in this baseline.
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
LABELS = [1, 2, 3, 4, 5, 6, 7]
RANDOM_STATE = 42


def load_split(path: Path, split_name: str) -> pd.DataFrame:
    """Load a split and check the columns and essential values."""
    if not path.exists():
        raise FileNotFoundError(
            f"Could not find {path}\n"
            "Check that your successful MMseqs2 split outputs are in "
            r"C:\protein_project_work\clustered"
        )

    df = pd.read_csv(path)

    if "label" not in df.columns:
        raise ValueError(
            f"{split_name} file has no 'label' column. "
            f"Available columns: {list(df.columns)}"
        )

    # Prefer the cleaned sequence column if available.
    sequence_col = "sequence_clean" if "sequence_clean" in df.columns else "Sequence"
    if sequence_col not in df.columns:
        raise ValueError(
            f"{split_name} file has neither 'sequence_clean' nor 'Sequence'. "
            f"Available columns: {list(df.columns)}"
        )

    df = df.copy()
    df["label"] = pd.to_numeric(df["label"], errors="coerce")
    df[sequence_col] = df[sequence_col].fillna("").astype(str).str.upper().str.strip()

    before = len(df)
    df = df.dropna(subset=["label"])
    df = df[df[sequence_col].str.len() > 0].copy()

    if len(df) != before:
        print(
            f"Warning: dropped {before - len(df)} rows from {split_name} "
            "because the label or sequence was missing/invalid."
        )

    df["label"] = df["label"].astype(int)
    unexpected = sorted(set(df["label"]) - set(LABELS))
    if unexpected:
        raise ValueError(
            f"Unexpected label(s) in {split_name}: {unexpected}. "
            f"Expected labels: {LABELS}"
        )

    df["_sequence_for_features"] = df[sequence_col]
    return df


def amino_acid_composition(sequence: str) -> np.ndarray:
    """Return 20 normalized amino-acid frequencies for one sequence."""
    sequence = str(sequence).upper()
    valid_residues = [residue for residue in sequence if residue in AMINO_ACIDS]

    if not valid_residues:
        return np.zeros(len(AMINO_ACIDS), dtype=np.float32)

    denominator = len(valid_residues)
    return np.array(
        [valid_residues.count(aa) / denominator for aa in AMINO_ACIDS],
        dtype=np.float32,
    )


def make_features(df: pd.DataFrame) -> np.ndarray:
    return np.vstack(
        [amino_acid_composition(seq) for seq in df["_sequence_for_features"]]
    )


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Python:", sys.version.split()[0])
    print("pandas:", pd.__version__)
    print("scikit-learn:", sklearn.__version__)
    print("\nLoading saved train and validation splits...")

    train_df = load_split(TRAIN_FILE, "train")
    val_df = load_split(VALIDATION_FILE, "validation")

    print(f"Training records:   {len(train_df)}")
    print(f"Validation records: {len(val_df)}")
    print("Training class counts:")
    print(train_df["label"].value_counts().sort_index().to_string())
    print("\nValidation class counts:")
    print(val_df["label"].value_counts().sort_index().to_string())

    # Guard against accidental overlap if a stable accession column exists.
    if "Entry" in train_df.columns and "Entry" in val_df.columns:
        overlap = set(train_df["Entry"].dropna()) & set(val_df["Entry"].dropna())
        if overlap:
            raise ValueError(
                f"Data leakage check failed: {len(overlap)} Entry values appear "
                "in both train and validation."
            )

    # Confirm the cluster-aware split remains intact.
    if "cluster_id" in train_df.columns and "cluster_id" in val_df.columns:
        overlap_clusters = (
            set(train_df["cluster_id"].dropna())
            & set(val_df["cluster_id"].dropna())
        )
        if overlap_clusters:
            raise ValueError(
                "Cluster leakage check failed: "
                f"{len(overlap_clusters)} cluster IDs appear in train and validation."
            )
        print("\nCluster overlap check: PASS (0 shared cluster IDs)")
    else:
        print(
            "\nNote: cluster_id column not found in both files; "
            "the script could not independently recheck cluster overlap."
        )

    print("\nExtracting amino-acid-composition features...")
    X_train = make_features(train_df)
    y_train = train_df["label"].to_numpy()
    X_val = make_features(val_df)
    y_val = val_df["label"].to_numpy()

    if not np.isfinite(X_train).all() or not np.isfinite(X_val).all():
        raise ValueError("Feature matrix contains NaN or infinite values.")

    # The scaler and classifier are fitted using training data only.
    model = Pipeline(
        steps=[
            ("scaler", StandardScaler()),
            (
                "classifier",
                LogisticRegression(
                    max_iter=3000,
                    class_weight="balanced",
                    random_state=RANDOM_STATE,
                ),
            ),
        ]
    )

    print("Training logistic regression on the training split only...")
    model.fit(X_train, y_train)

    print("Evaluating on validation split...")
    val_predictions = model.predict(X_val)

    metrics = {
        "model": "amino_acid_composition_logistic_regression",
        "feature_type": "normalized_frequency_of_20_standard_amino_acids",
        "number_of_features": len(AMINO_ACIDS),
        "amino_acid_order": list(AMINO_ACIDS),
        "random_state": RANDOM_STATE,
        "training_records": int(len(train_df)),
        "validation_records": int(len(val_df)),
        "test_set_used": False,
        "accuracy": float(accuracy_score(y_val, val_predictions)),
        "balanced_accuracy": float(balanced_accuracy_score(y_val, val_predictions)),
        "macro_f1": float(f1_score(y_val, val_predictions, average="macro", zero_division=0)),
        "weighted_f1": float(f1_score(y_val, val_predictions, average="weighted", zero_division=0)),
    }

    report = classification_report(
        y_val,
        val_predictions,
        labels=LABELS,
        output_dict=True,
        zero_division=0,
    )
    report_text = classification_report(
        y_val,
        val_predictions,
        labels=LABELS,
        digits=4,
        zero_division=0,
    )
    cm = confusion_matrix(y_val, val_predictions, labels=LABELS)

    (OUTPUT_DIR / "validation_metrics.json").write_text(
        json.dumps(metrics, indent=2), encoding="utf-8"
    )
    (OUTPUT_DIR / "validation_classification_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    (OUTPUT_DIR / "validation_classification_report.txt").write_text(
        report_text, encoding="utf-8"
    )
    pd.DataFrame(
        cm,
        index=[f"actual_{label}" for label in LABELS],
        columns=[f"predicted_{label}" for label in LABELS],
    ).to_csv(OUTPUT_DIR / "validation_confusion_matrix.csv", index=True)

    fig, ax = plt.subplots(figsize=(8, 7))
    display = ConfusionMatrixDisplay(
        confusion_matrix=cm,
        display_labels=LABELS,
    )
    display.plot(ax=ax, values_format="d", colorbar=False)
    ax.set_title("Validation Confusion Matrix\nAmino-acid Composition + Logistic Regression")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "validation_confusion_matrix.png", dpi=180)
    plt.close(fig)

    # Save validation predictions to support error analysis.
    prediction_df = pd.DataFrame(
        {
            "Entry": val_df["Entry"].to_numpy() if "Entry" in val_df.columns else np.arange(len(val_df)),
            "true_label": y_val,
            "predicted_label": val_predictions,
            "correct": y_val == val_predictions,
        }
    )
    prediction_df.to_csv(OUTPUT_DIR / "validation_predictions.csv", index=False)

    print("\n========== VALIDATION RESULTS ==========")
    print(f"Accuracy:          {metrics['accuracy']:.4f}")
    print(f"Balanced accuracy: {metrics['balanced_accuracy']:.4f}")
    print(f"Macro F1:           {metrics['macro_f1']:.4f}")
    print(f"Weighted F1:        {metrics['weighted_f1']:.4f}")
    print("\nPer-class report:")
    print(report_text)
    print(f"\nOutputs saved to: {OUTPUT_DIR}")
    print("The test set was not loaded or evaluated.")
    print("Next: inspect validation results before choosing the next experiment.")


if __name__ == "__main__":
    main()

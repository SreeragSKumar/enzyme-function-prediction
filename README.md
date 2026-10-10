# Low-Data Protein Language Model for Enzyme Function Prediction

## Project Overview

This project develops a machine-learning pipeline for predicting enzyme function from protein sequences. The current task classifies protein sequences into seven top-level Enzyme Commission (EC) classes.

The approach begins with interpretable sequence-feature baselines and will later compare them with protein language model representations such as ESM-2.

## Current Progress

- Collected and audited reviewed enzyme records across EC classes 1–7.
- Prepared a model-ready dataset containing 2,673 protein sequences.
- Clustered sequences using MMseqs2 at 30% minimum sequence identity and 80% coverage.
- Created cluster-aware training, validation, and test splits.
- Trained an amino-acid-composition logistic regression baseline.
- Evaluated the baseline on the validation split.

### Dataset Splits

| Split | Records |
|---|---:|
| Training | 1,757 |
| Validation | 377 |
| Test | 539 |
| **Total** | **2,673** |

Clusters do not overlap between the three splits.

### Initial Baseline Results

| Metric | Validation result |
|---|---:|
| Accuracy | 35.01% |
| Balanced accuracy | 35.46% |
| Macro F1 | 0.3427 |
| Weighted F1 | 0.3386 |

These results are from validation data only. The held-out test set has not been evaluated.

## Technology Stack

- Python
- pandas and NumPy
- scikit-learn
- Matplotlib
- MMseqs2
- Planned: ESM-2 protein embeddings

## Project Structure

```text
project sem 9/
├── .gitignore
├── README.md
├── 03_cluster_and_split.py
├── 04_baseline_model.py
├── protein_project/
│   ├── raw/
│   ├── processed/
│   ├── audit/
│   ├── metadata/
│   ├── checkpoints/
│   ├── figures/
│   └── archive/
└── protein_project_complete_progress.md
```

Additional data splits and experiment outputs are stored under `C:\protein_project_work\clustered` and `C:\protein_project_work\baseline`.

## Environment Setup

The working baseline environment uses Python 3.10.11 and scikit-learn 1.7.2.

Activate the virtual environment in PowerShell:

```powershell
.\.venv310\Scripts\Activate.ps1
```

Install the core Python dependencies if needed:

```powershell
python -m pip install pandas numpy matplotlib scikit-learn
```

MMseqs2 is required for the clustering stage and must be installed separately.

## Experimental Roadmap

1. Establish the amino-acid-composition baseline.
2. Implement a k-mer baseline.
3. Generate ESM-2 embeddings and train a downstream classifier.
4. Compare models using the same cluster-aware splits and validation metrics.
5. Select the final approach using validation data.
6. Evaluate the selected approach on the held-out test set.

## Reproducibility and Data Handling

- Preserve the saved cluster-aware splits across experiments.
- Fit feature transformations on training data only.
- Use validation data for model selection.
- Do not use the test set for tuning.
- Keep large datasets, virtual environments, temporary files, and credentials out of Git.
- Consult `protein_project_complete_progress.md` for the detailed experiment history, audit statistics, and baseline analysis.

## Current Status

**Completed:** dataset preparation, sequence clustering, split integrity checks, and the initial validation baseline.

**Next:** implement the k-mer baseline in a separate script.
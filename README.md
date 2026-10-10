# Low-Data Protein Language Model for Enzyme Function Prediction

## Project Overview

This project investigates machine-learning methods for predicting enzyme function from protein sequences. The current task classifies protein sequences into seven top-level Enzyme Commission (EC) classes (labels 1–7).

The project compares conventional sequence-feature baselines with embeddings from a pretrained protein language model, ESM-2. The current best reported model is an ESM-2 embedding classifier using a tuned Linear Support Vector Machine (SVM).

> **Current best result:** 71.43% held-out test accuracy and 0.7035 macro F1. The model was selected using validation data; the held-out test set was evaluated once after model selection.

## Dataset and Splits

The model-ready dataset contains **2,673 protein sequences** collected from reviewed UniProt enzyme records across top-level EC classes 1–7.

Sequences were clustered using MMseqs2 with a minimum sequence identity of 30% and coverage of 80%. The dataset was then split by cluster to prevent clusters from crossing between train, validation, and test partitions.

| Split | Records | Share |
|---|---:|---:|
| Training | 1,757 | 65.7% |
| Validation | 377 | 14.1% |
| Test | 539 | 20.2% |
| **Total** | **2,673** | **100%** |

The saved integrity checks reported zero cluster overlap between train, validation, and test. Because complete clusters were kept together, the final split proportions differ slightly from the original target proportions.

## Methodology

1. **Data collection and quality control:** review annotations, sequence validity, class labels, and sequence lengths.
2. **Cluster-aware splitting:** cluster related sequences with MMseqs2 and assign whole clusters to one split.
3. **Conventional baselines:** amino-acid composition and k-mer representations (including 3-mer, 4-mer, and 5-mer TF-IDF experiments) with scikit-learn classifiers.
4. **Protein language model embeddings:** use `facebook/esm2_t12_35M_UR50D` to generate 480-dimensional, mean-pooled embeddings. Long sequences are processed in windows and the window embeddings are averaged into a protein-level representation.
5. **Model selection:** compare downstream classifiers and tune the Linear SVM using validation macro F1.
6. **Final evaluation:** after selecting the configuration, fit the final classifier using the combined train and validation embeddings, then evaluate once on the held-out test set.

## Results

### ESM-2 validation comparison

The following models were compared on the validation split. The test set was not used for this comparison or for selecting the hyperparameter.

| Model | Accuracy | Balanced accuracy | Macro F1 | Weighted F1 |
|---|---:|---:|---:|---:|
| ESM-2 + Logistic Regression | 0.6419 | 0.6342 | 0.6333 | 0.6373 |
| ESM-2 + initial Linear SVM | 0.6525 | 0.6443 | 0.6399 | 0.6441 |
| ESM-2 + tuned Linear SVM (`C = 2.0`) | 0.6605 | 0.6522 | 0.6463 | 0.6506 |

The tuned Linear SVM used `class_weight="balanced"` and `random_state=42`. `C = 2.0` was selected based on validation macro F1.

### Final held-out test results

The selected SVM was fitted on **2,134 records** (training plus validation) and evaluated on **539 held-out test records**.

| Metric | Test result |
|---|---:|
| Accuracy | 0.7143 (71.43%) |
| Balanced accuracy | 0.7112 (71.12%) |
| Macro F1 | 0.7035 |
| Weighted F1 | 0.7093 |

#### Per-class test performance

| EC top-level label | Precision | Recall | F1-score | Support |
|---:|---:|---:|---:|---:|
| 1 | 0.63 | 0.58 | 0.60 | 72 |
| 2 | 0.82 | 0.82 | 0.82 | 85 |
| 3 | 0.66 | 0.79 | 0.72 | 62 |
| 4 | 0.56 | 0.59 | 0.58 | 75 |
| 5 | 0.56 | 0.43 | 0.49 | 81 |
| 6 | 0.76 | 0.86 | 0.81 | 71 |
| 7 | 0.90 | 0.90 | 0.90 | 93 |

Class 7 performed best by F1, while class 5 was the most difficult class to recover. These results measure top-level EC-class classification on this dataset and split; they do not experimentally verify enzyme activity or establish full EC-number correctness.

**Test-set protocol:** the held-out test set was not used for model selection, but it has now been evaluated once. Do not continue tuning this model against the same test results. Use validation data for further model development; a future independent study should ideally use a fresh external or temporally separated test set.

## Technology Stack

- Python
- pandas and NumPy
- scikit-learn
- Matplotlib
- MMseqs2
- PyTorch
- Hugging Face Transformers
- ESM-2 protein language model
- Kaggle GPU notebook (NVIDIA Tesla T4 was reported for embedding extraction)

## Project Structure

```text
project sem 9/
├── README.md
├── 03_cluster_and_split.py
├── 04_baseline_model.py
├── 06_4mer_baseline.py
├── 07_3mer_svm_baseline.py
├── 08_tune_3mer_svm.py
├── 09_5mer_svm_baseline.py
├── protein_project_complete_progress.md
└── protein_project/
    ├── raw/
    ├── processed/
    ├── audit/
    ├── metadata/
    ├── checkpoints/
    ├── figures/
    ├── archive/
    ├── notebooks/
    │   └── esm2-enzyme-function-prediction.ipynb
    └── results/
        └── esm2_final/
            ├── final_test_metrics.json
            ├── final_test_classification_report.csv
            ├── final_test_confusion_matrix.csv
            └── final_test_confusion_matrix.png
```

Additional clustering and baseline outputs have been generated under `C:\protein_project_work\clustered` and `C:\protein_project_work\baseline`. These are external working directories and are not automatically part of the Git repository.

Kaggle files saved under `/kaggle/working/` are temporary unless downloaded or copied to persistent storage. The repository notebook and saved result files are the durable references for the experiment; do not assume that all cached embedding arrays or model weights are committed.

## Environment Setup

The original local baseline environment used Python 3.10.11 and scikit-learn 1.7.2. The ESM-2 experiments were run in Kaggle with GPU acceleration and the package versions recorded in the experiment log.

Activate the local virtual environment in PowerShell:

```powershell
.\.venv310\Scripts\Activate.ps1
```

Install the core local baseline dependencies if needed:

```powershell
python -m pip install pandas numpy matplotlib scikit-learn
```

MMseqs2 is required for the clustering stage and must be installed separately. The ESM-2 notebook additionally requires PyTorch and Hugging Face Transformers, along with access to the pretrained model. For exact settings, refer to the notebook and the detailed progress log.

## Reproducibility and Data Handling

- Preserve the saved cluster-aware train, validation, and test splits across experiments.
- Fit feature transformations on training data only.
- Use validation data for model selection and hyperparameter tuning.
- Do not use the existing test results to tune the current model.
- Record the ESM-2 model identifier, windowing and pooling settings, classifier configuration, random seed, and package versions.
- Keep credentials, virtual environments, temporary files, large downloaded model weights, and raw data that should not be redistributed out of Git.
- Consult `protein_project_complete_progress.md` for the detailed experiment history, data-quality audit, baseline results, troubleshooting notes, and next steps.

## Limitations and Next Steps

- Audit the eight inherited `is_exact_duplicate=True` flags identified during data review, even though the current filtered sequences had no repeated sequence strings or hashes in the earlier integrity check.
- Perform error analysis, especially for class 5 and the confusions shown in the final test confusion matrix.
- Verify that the notebook and saved outputs can be reproduced from documented inputs and settings.
- Run a planned **learning-curve experiment** using progressively smaller, cluster-aware subsets of training data while keeping validation fixed. This is still future work; the current full-training result alone does not establish low-data performance.
- Consider a fresh external or temporally separated test set if further model development is performed.

## Current Status

**Completed:** dataset preparation and audit, MMseqs2 clustering, cluster-aware splits, amino-acid-composition baseline, k-mer experiments, ESM-2 embedding extraction, validation-based Linear SVM tuning, and one final held-out test evaluation.

**Remaining:** duplicate-flag audit, detailed error analysis, reproducibility review, final research write-up, and the planned low-data learning-curve experiment.

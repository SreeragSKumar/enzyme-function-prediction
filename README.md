# Low-Data Protein Language Model for Enzyme Function Prediction and Candidate Prioritization

## Project overview

This project evaluates protein-sequence representations for predicting one of seven top-level Enzyme Commission (EC) classes. It compares conventional sequence-feature baselines with embeddings from the pretrained ESM-2 model `facebook/esm2_t12_35M_UR50D` and studies performance under reduced labelled-data availability.

The project also explores how predictions might support candidate prioritization. Predictions are computational hypotheses, not experimental confirmation of enzyme activity or proof of novel function.

## Dataset and split design

- Model-ready reviewed UniProt sequences: **2,673**
- Top-level EC classes: 1–7
- MMseqs2 clustering parameters: minimum sequence identity 30%, coverage 80%
- Total clusters before splitting: 1,427
- Cluster-aware partitions: training 1,757; validation 377; test 539
- No cluster overlap between train, validation, and test partitions was detected.

| Split | Proteins | Share |
|---|---:|---:|
| Training | 1,757 | 65.7% |
| Validation | 377 | 14.1% |
| Test | 539 | 20.2% |
| **Total** | **2,673** | **100%** |

The split proportions differ slightly from nominal targets because clusters were kept intact.

## Main model and evaluation

ESM-2 embeddings have 480 features per protein. The selected downstream classifier is a scaled Linear SVM (`LinearSVC`, `C=2.0`, `class_weight="balanced"`, `max_iter=10000`). Model selection used validation macro F1. The final classifier was fit on the combined training and validation embeddings (2,134 proteins) and evaluated on the held-out test set. The test set was not used for hyperparameter selection.

| Test metric | Result |
|---|---:|
| Accuracy | 71.43% |
| Balanced accuracy | 71.12% |
| Macro F1 | 0.7035 |
| Weighted F1 | 0.7093 |

The per-class test F1 scores were: class 1 0.60, class 2 0.82, class 3 0.72, class 4 0.58, class 5 0.49, class 6 0.81, and class 7 0.90. Class 5 was the weakest by F1; its recall was 0.43. These results apply to this dataset and split and do not establish full EC-number correctness or experimental enzyme activity.

## Low-data experiments

Two learning-curve experiments were run using saved ESM-2 embeddings and the fixed 377-protein validation set. Each fraction was evaluated with seeds 42, 123, and 2026. Metrics below are mean validation macro F1 across seeds; the 100% subset is the complete training split, so its score is identical across seeds.

### Protein-level stratified sampling

| Target fraction | Proteins selected | Linear SVM | Logistic regression |
|---:|---:|---:|---:|
| 10% | 173 | 0.5281 | 0.5336 |
| 25% | 436 | 0.5642 | 0.5912 |
| 50% | 877 | 0.5681 | 0.5878 |
| 100% | 1,757 | 0.6463 | 0.6333 |

### Cluster-aware sampling

Whole clusters were selected rather than individual proteins. Because cluster sizes differ, the actual protein counts vary around the target fractions.

| Target fraction | Mean proteins selected | Mean clusters selected | Linear SVM | Logistic regression |
|---:|---:|---:|---:|---:|
| 10% | 187.3 | 103.0 | 0.4675 | 0.4336 |
| 25% | 440.0 | 234.3 | 0.5753 | 0.5743 |
| 50% | 879.7 | 466.0 | 0.5825 | 0.5982 |
| 100% | 1,757 | 914 | 0.6463 | 0.6333 |

These learning curves show an overall improvement with more training data, although gains are not linear at every fraction and small subsets vary across seeds. The cluster-aware experiment is more aligned with the low-data research question, but it is still evaluated on the same fixed validation set. It does not by itself prove generalization to entirely novel enzyme families.

### Experiment outputs

Kaggle output paths recorded during the experiments:

- `/kaggle/working/esm2_results/low_data/low_data_runs.csv`
- `/kaggle/working/esm2_results/low_data/low_data_summary.csv`
- `/kaggle/working/esm2_results/low_data/low_data_learning_curve.png`
- `/kaggle/working/esm2_results/cluster_low_data/cluster_low_data_runs.csv`
- `/kaggle/working/esm2_results/cluster_low_data/cluster_low_data_summary.csv`
- `/kaggle/working/esm2_results/cluster_low_data/cluster_low_data_learning_curve.png`

Kaggle `/kaggle/working` files are session-local unless downloaded or copied to persistent storage. These paths identify the outputs created in the notebook; they should not be assumed to be in Git unless explicitly committed.

## Data integrity audit

The model-ready CSV contains 2,673 unique accessions. A read-only audit found:

- 2,673 unique `exact_sequence_hash` values
- 2,673 unique raw `Sequence` values
- 2,673 unique `sequence_clean` values
- 8 records with inherited `is_exact_duplicate=True` flags
- 0 exact duplicate rows when recalculating duplicates from the current cleaned sequence values

The eight flags appear inconsistent with the current sequence values and may be stale metadata. The original CSV was not modified. Audit reports were saved locally under `protein_project/results/data_audit/` as `previously_flagged_duplicates.csv` and `duplicate_flag_audit.csv`. Recalculate or document the flag's provenance before using it as a quality-control exclusion.

## External UniProt evaluation attempt

A separate UniProt collection attempt produced 179 unique accessions after excluding the original dataset. Only 2 records had a single top-level EC class under the parsing rule used for the initial external evaluation; 177 had annotations spanning multiple top-level classes. The two unambiguous records were insufficient for meaningful external performance metrics. **No valid external performance evaluation was completed.** The 179-record pool should not be presented as a balanced benchmark, and the queried class must not be treated as ground truth when a protein has multiple EC annotations.

## Interpretation and limitations

- The model predicts a broad top-level EC class, not a complete EC number, substrate specificity, or experimentally validated function.
- The held-out test set has been evaluated; do not use it for further tuning. For a new final unbiased evaluation after substantial changes, obtain a fresh external or temporally separated test set.
- Cluster-aware splitting reduces leakage at the chosen MMseqs2 thresholds but cannot guarantee complete absence of remote homology or biological relatedness.
- The low-data curves use the existing validation set. Cluster-aware subsets preserve whole training clusters, but the experiment is not an independent test of entirely unseen protein families.
- The external collection attempt did not produce enough unambiguous records for evaluation.
- Candidate-ranking scores, confidence, and sequence novelty are computational evidence for prioritization, not proof of novelty or catalytic activity. A suitable separate candidate set and validated ranking procedure are still needed.
- The eight inherited duplicate flags require provenance review; do not silently delete those records.

## Technology and reproducibility

- Python, pandas, NumPy, scikit-learn, Matplotlib
- MMseqs2 for sequence clustering
- ESM-2: `facebook/esm2_t12_35M_UR50D`
- Baseline local environment previously recorded: Python 3.10.11, scikit-learn 1.7.2
- Kaggle experiment environment recorded: PyTorch 2.11.0+cu128, Transformers 5.16.1, pandas 2.3.3, scikit-learn 1.6.1

Keep the cluster-aware splits fixed, fit transformations on training data, use validation data for model selection, and preserve the test set from further tuning. Do not commit raw datasets, large embedding arrays, model weights, virtual environments, or credentials unless there is a deliberate, documented reason.

## Repository and result references

- ESM-2 notebook: `protein_project/notebooks/esm2-enzyme-function-prediction.ipynb`
- Final ESM-2 results: `protein_project/results/esm2_final/`
- Final test metric file: `final_test_metrics.json`
- Per-class test report: `final_test_classification_report.csv`
- Confusion matrix: `final_test_confusion_matrix.csv` and `final_test_confusion_matrix.png`
- Detailed history and experiment log: `protein_project_complete_progress.md`

## Current status and next steps

**Completed:** dataset preparation and audit; MMseqs2 clustering and cluster-aware splitting; amino-acid composition and k-mer baselines; ESM-2 embedding extraction; validation-based SVM selection; one final held-out test evaluation; protein-level and cluster-aware low-data experiments; initial duplicate-flag audit.

**Remaining:** reconcile the duplicate flag's provenance; ensure low-data CSVs/plots are downloaded or committed if they are to be part of the repository; finalize error analysis, methodology, limitations, and report; decide whether to implement a minimal candidate-ranking prototype using a separate suitable candidate set. The attempted external UniProt pool is not sufficient for external benchmarking.

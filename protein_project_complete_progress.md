# Complete Project Progress and Experiment Log
## Low-Data Protein Language Model for Enzyme Function Prediction and Novel Candidate Discovery

**Last updated:** 9 October 2026
**Status:** Dataset preparation, MMseqs2 clustering, cluster-aware splitting, amino-acid-composition and k-mer baselines, and ESM-2 embedding experiments are complete. The selected ESM-2 + tuned Linear SVM was evaluated on the held-out test split after model selection. Remaining work includes reproducibility cleanup, error analysis, and documenting limitations and future low-data experiments.

## 1. Project objective

Develop a low-data machine-learning pipeline to predict enzyme function from protein sequences. Begin with transparent sequence-feature baselines, then compare them with protein language model representations such as ESM-2. The longer-term goal is to help prioritize candidate proteins for further investigation.

The current task predicts one of seven top-level Enzyme Commission (EC) classes, represented by labels 1–7.

Planned sequence:
1. Collect and audit reviewed enzyme sequences.
2. Clean records and create a model-ready dataset.
3. Cluster similar sequences and create cluster-aware train/validation/test splits.
4. Train an amino-acid-composition baseline.
5. Train a k-mer baseline.
6. Compare with ESM-2 embeddings.
7. Select models using validation data and evaluate the selected approach on the held-out test set.
8. Document results, limitations, and reproducibility.

Dataset preparation, cluster-aware splitting, amino-acid-composition and k-mer baselines, and ESM-2 evaluation have now been completed. The detailed chronology below includes earlier planning statements; this addendum records the latest ESM-2 results and supersedes outdated statements that the test set has not been evaluated.

## 2. Development environment and troubleshooting

- OS: Windows
- Editor/terminal: VS Code and PowerShell
- Main project directory: `C:\Users\sreer\OneDrive\Desktop\project sem 9`
- Sequence clustering: MMseqs2

### Python DLL/Application Control issue

The original `.venv` used Python 3.14.2 and scikit-learn 1.9.1. Importing `sklearn.metrics` failed because Windows Code Integrity/Application Control blocked `sklearn\metrics\_dist_metrics.cp314-win_amd64.pyd`. Event ID 3077 reported that the component did not meet enterprise signing requirements or violated the active policy. Policy ID:

`0283ac0f-fff1-49ae-ada1-8a933130cad6`

Python 3.10 was tested and successfully imported scikit-learn metrics. A separate `.venv310` environment was created; the original `.venv` was left in place.

Successful baseline environment:
- Python 3.10.11
- pandas 2.3.3
- scikit-learn 1.7.2
- NumPy and Matplotlib available

Commands used:

```powershell
py -3.10 -m venv .venv310
.\.venv310\Scripts\Activate.ps1
python -m pip install pandas numpy matplotlib scikit-learn
python -c "import pandas, numpy, matplotlib, sklearn; from sklearn.metrics import accuracy_score; print('Python environment OK'); print('scikit-learn:', sklearn.__version__)"
```

Activate `.venv310` before future runs:

```powershell
.\.venv310\Scripts\Activate.ps1
```

## 3. Dataset collection and audit

The dataset was expanded from UniProt reviewed enzyme records across seven top-level EC classes. The collection target was up to 10,000 records but was not reached. The final model-ready dataset contains **2,673 records**.

### Earlier annotation-audit counts

| Audit item | Count |
|---|---:|
| Raw unique accessions | 3,393 |
| Single-class annotations | 2,759 |
| Missing class annotations | 448 |
| Multiple-class annotations | 186 |
| Retrieval/annotation mismatches | 475 |

These are earlier audit-stage counts, not counts from the final model-ready table.

### Sequence QC

| QC item | Result |
|---|---:|
| Empty sequences | 0 |
| Invalid-character sequences | 0 |
| Sequences with ambiguous residues | 13 |
| Records in exact-duplicate groups at earlier QC stage | 20 |
| Sequence-length/annotation mismatches | 0 |
| Final model-ready records | 2,673 |

### Class counts

| Top-level EC class | Before filtering | After filtering |
|---:|---:|---:|
| 1 | 400 | 385 |
| 2 | 437 | 412 |
| 3 | 319 | 306 |
| 4 | 379 | 376 |
| 5 | 395 | 385 |
| 6 | 420 | 410 |
| 7 | 409 | 399 |
| **Total** | **2,759** | **2,673** |

Important fields include `Entry`, `Sequence`, `EC number`, `retrieval_ec_class`, `label`, `annotation_class_status`, `retrieval_class_matches_annotation`, `sequence_clean`, `exact_sequence_hash`, `is_exact_duplicate`, and `length_matches_annotation`.

### Dataset directories

```text
protein_project/
├── archive/pilot/
├── figures/
├── metadata/
│   ├── expanded_dataset_metadata.json
│   └── uniprot_query_log.json
├── audit/
├── processed/
│   ├── uniprot_expanded_labelled.csv
│   └── uniprot_expanded_model_ready.csv
├── raw/
│   └── uniprot_expanded_raw.tsv
└── checkpoints/
    ├── ec_1 through ec_7
    └── uniprot_expansion_checkpoint.tsv
```

Pilot files, raw data, checkpoints, audit files, and metadata were retained for traceability.

Model-ready file: `protein_project\processed\uniprot_expanded_model_ready.csv`

## 4. MMseqs2 clustering

MMseqs2 was used before splitting to reduce the chance that highly similar proteins appear across training and evaluation splits.

- Executable: `C:\tools\mmseqs2\mmseqs\mmseqs.bat`
- Version string: `8cc5ce367b5638c4306c2d7fcf652dd099a4643f`
- Method: `easy-cluster`
- Minimum sequence identity: `--min-seq-id 0.3`
- Coverage: `-c 0.8`
- Coverage mode: `--cov-mode 0`
- Threads: `--threads 4`

Troubleshooting included a pandas row-attribute issue involving `_fasta_id` and a Windows path-with-spaces issue. The successful workflow used a separate working path:

- FASTA input: `C:\protein_project_work\clustered\fasta\proteins_for_mmseqs.fasta`
- Temporary directory: `C:\protein_project_work\mmseqs_tmp`
- Outputs: `C:\protein_project_work\clustered`

### Results

| Metric | Result |
|---|---:|
| Input sequences | 2,673 |
| Clusters | 1,427 |
| Largest cluster | 37 |
| Singleton clusters | 989 |

Important output files include `train.csv`, `validation.csv`, `test.csv`, `uniprot_clustered.csv`, `sequence_cluster_membership.csv`, `split_class_counts.csv`, `split_metadata.json`, `fasta_id_mapping.csv`, `mmseqs_clusters_cluster.tsv`, `proteins_for_mmseqs.fasta`, `mmseqs_clusters_all_seqs.fasta`, and `mmseqs_clusters_rep_seq.fasta`.

## 5. Cluster-aware splits and integrity checks

### Split sizes

| Split | Records | Proportion |
|---|---:|---:|
| Train | 1,757 | 65.7% |
| Validation | 377 | 14.1% |
| Test | 539 | 20.2% |
| **Total** | **2,673** | **100%** |

The proportions differ slightly from the approximate 64%/16%/20% target because whole clusters were kept together. Retain these splits rather than randomizing records across cluster boundaries.

### Class counts by split

| Label | Train | Validation | Test | Total |
|---:|---:|---:|---:|---:|
| 1 | 253 | 60 | 72 | 385 |
| 2 | 263 | 64 | 85 | 412 |
| 3 | 203 | 41 | 62 | 306 |
| 4 | 246 | 55 | 75 | 376 |
| 5 | 252 | 52 | 81 | 385 |
| 6 | 282 | 57 | 71 | 410 |
| 7 | 258 | 48 | 93 | 399 |
| **Total** | **1,757** | **377** | **539** | **2,673** |

### Checks that passed
- Split totals sum to 2,673.
- There are 2,673 unique `Entry` values across the combined splits.
- All labels 1–7 appear in each split.
- Cluster overlap: train/validation = 0; train/test = 0; validation/test = 0.
- Membership CSV and MMseqs2 TSV assignments agree.
- Entry-to-cluster mappings in each split agree with the saved mapping.
- Class-count reports agree with the split CSVs.
- The final clustered dataset was reported to have single-class annotations, matching retrieval/annotation labels, no invalid or empty sequences, zero ambiguous residues, and matching sequence lengths.

These checks support split integrity but do not eliminate every possible biological or annotation bias.

### Duplicate-flag caveat
Eight rows in `uniprot_clustered.csv` have `is_exact_duplicate=True`, but the final dataset has zero repeated `exact_sequence_hash` groups and zero duplicate raw `Sequence` or cleaned `sequence_clean` strings. The flag likely came from an earlier QC stage and was not recalculated after filtering. Investigate and document or recalculate it; do not assume these flags prove duplicates remain.

## 6. Experiment A: amino-acid composition + logistic regression

### Script and paths
- Script: `04_baseline_model.py`
- Inputs: `C:\protein_project_work\clustered\train.csv` and `validation.csv`
- Outputs: `C:\protein_project_work\baseline`

### Method
Each sequence is represented by normalized frequencies of the 20 standard amino acids:

`ACDEFGHIKLMNPQRSTVWY`

This creates 20 features per protein and ignores residue order. A scikit-learn pipeline applies `StandardScaler` and logistic regression with `max_iter=3000`, `class_weight="balanced"`, and `random_state=42`.

The model was fitted on training data only and evaluated on validation data. The script checked that train and validation had no shared cluster IDs. The test CSV was not loaded or evaluated.

### Overall validation metrics

| Metric | Result |
|---|---:|
| Accuracy | 0.3501 (35.01%) |
| Balanced accuracy | 0.3546 (35.46%) |
| Macro F1 | 0.3427 |
| Weighted F1 | 0.3386 |
| Validation records | 377 |
| Test evaluated? | **No** |

### Per-class results

| Label | Precision | Recall | F1-score | Support |
|---:|---:|---:|---:|---:|
| 1 | 0.3333 | 0.2667 | 0.2963 | 60 |
| 2 | 0.3333 | 0.4688 | 0.3896 | 64 |
| 3 | 0.3061 | 0.3659 | 0.3333 | 41 |
| 4 | 0.2258 | 0.1273 | 0.1628 | 55 |
| 5 | 0.3261 | 0.2885 | 0.3061 | 52 |
| 6 | 0.3036 | 0.2982 | 0.3009 | 57 |
| 7 | 0.5614 | 0.6667 | 0.6095 | 48 |
| Macro average | 0.3414 | 0.3546 | 0.3427 | 377 |
| Weighted average | 0.3382 | 0.3501 | 0.3386 | 377 |

### Confusion matrix

Rows are actual labels and columns are predicted labels.

| Actual \ Predicted | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 16 | 20 | 3 | 7 | 1 | 11 | 2 |
| 2 | 6 | 30 | 11 | 1 | 6 | 6 | 4 |
| 3 | 3 | 11 | 15 | 3 | 2 | 5 | 2 |
| 4 | 5 | 12 | 11 | 7 | 9 | 7 | 4 |
| 5 | 2 | 10 | 6 | 7 | 15 | 6 | 6 |
| 6 | 11 | 4 | 3 | 3 | 12 | 17 | 7 |
| 7 | 5 | 3 | 0 | 3 | 1 | 4 | 32 |

### Interpretation
- Class 7 is strongest: 32/48 correct; recall 66.7%; F1 0.6095.
- Class 4 is weakest: 7/55 correct; recall 12.7%; F1 0.1628.
- Class 1 is often predicted as class 2 (20 cases) and class 6 (11 cases).
- Class 4 is often predicted as class 2 (12) or class 3 (11).
- Class 6 is often predicted as class 5 (12) or class 1 (11).
- Class 2 and class 3 are confused with each other (11 cases in each direction).

These are observed prediction confusions; the matrix alone does not establish their biological cause.

The baseline has learned some signal but is limited because amino-acid composition ignores sequence order, motifs, domains, and long-range relationships. The 35.01% validation accuracy is above the 14.3% uniform-random reference for seven classes, but a majority-class or stratified baseline would add context. This result alone does not establish biological usefulness.

### Saved output files
Under `C:\protein_project_work\baseline`:
- `validation_metrics.json`
- `validation_classification_report.json`
- `validation_classification_report.txt`
- `validation_confusion_matrix.csv`
- `validation_confusion_matrix.png`
- `validation_predictions.csv`

## 7. Completed checklist

- [x] Expanded reviewed enzyme dataset across EC classes 1–7.
- [x] Recorded audit, query log, metadata, and checkpoints.
- [x] Performed sequence QC and filtering.
- [x] Created a 2,673-record model-ready dataset.
- [x] Installed and verified MMseqs2.
- [x] Clustered 2,673 sequences into 1,427 clusters.
- [x] Created cluster-aware train/validation/test splits.
- [x] Verified zero cluster overlap across splits and checked mapping/count consistency.
- [x] Created `04_baseline_model.py`.
- [x] Resolved the blocked scikit-learn import using Python 3.10 in `.venv310`.
- [x] Trained and evaluated the amino-acid-composition baseline on validation data.
- [x] Saved and inspected the validation metrics and confusion matrix.

## 8. Caveats and safeguards

1. Investigate the eight inherited `is_exact_duplicate=True` flags.
2. Preserve the cluster-aware splits and report the achieved proportions accurately.
3. Do not use the test set for feature selection, hyperparameter tuning, or model selection.
4. A top-level EC-class prediction is a broad label, not a complete function annotation or proof of catalytic activity.
5. Clustering reduces similarity leakage at the chosen settings but cannot guarantee elimination of all homology-related bias.
6. Record scripts, package versions, parameters, and metrics for each experiment.
7. Use the same saved splits and comparable metrics for each model.

## 9. Next steps

### Experiment B — k-mer baseline (next)
Create a new script, `05_kmer_baseline.py`; do not overwrite `04_baseline_model.py`.

Planned approach:
- Represent sequences by short consecutive amino-acid patterns, initially 3-mers.
- Use a sparse bag-of-k-mers or TF-IDF representation.
- Fit the vectorizer and classifier using training data only.
- Transform validation data using the training-fitted vectorizer.
- Use validation data for model comparison.
- Save metrics, reports, predictions, and a confusion matrix under a separate directory such as `C:\protein_project_work\kmer_baseline`.
- Compare macro F1, balanced accuracy, per-class recall, and confusion matrices against Experiment A.

### Experiment C — ESM-2
After the simple baselines:
- Select an ESM-2 model that fits available hardware and memory.
- Generate protein embeddings.
- Train a downstream classifier using the training split.
- Use validation data for model selection.
- Compare with both baselines using the same partitions and metrics.

### Final evaluation
Once model selection is complete, follow a documented final-training procedure and evaluate the held-out test split for the final report. Do not present validation results as test results.

## 10. Path reference

| Purpose | Path |
|---|---|
| Main project directory | `C:\Users\sreer\OneDrive\Desktop\project sem 9` |
| Python 3.10 environment | `.venv310` inside main project directory |
| Model-ready dataset | `protein_project\processed\uniprot_expanded_model_ready.csv` |
| Cluster-aware splits | `C:\protein_project_work\clustered\` |
| Baseline script | `04_baseline_model.py` |
| Baseline outputs | `C:\protein_project_work\baseline\` |
| Next planned script | `05_kmer_baseline.py` (not yet created or run) |

## 11. Handover summary

The project has 2,673 quality-checked enzyme sequences across seven top-level EC classes. MMseqs2 produced 1,427 clusters, and cluster-aware splits contain 1,757 training, 377 validation, and 539 test records. The cluster-overlap and other core split-integrity checks passed.

The first model, normalized 20-amino-acid composition with logistic regression, ran successfully in Python 3.10.11 and achieved validation accuracy 35.01%, balanced accuracy 35.46%, and macro F1 0.3427. Class 7 was strongest and class 4 weakest. The test set has not been used.

**Historical next action:** build and evaluate the k-mer baseline. This was subsequently completed; see the k-mer experiment entries and ESM-2 addendum below.

## Baseline Experiment 1: Amino-Acid Composition + Logistic Regression

**Status:** Completed — validation evaluation finished.

### Objective
Establish a baseline for seven-class enzyme function prediction using amino-acid composition features and a Logistic Regression classifier.

### Dataset and split
- Training records: 1,757
- Validation records: 377
- Number of enzyme classes: 7
- Split integrity: Passed; zero shared cluster IDs between training and validation.
- The existing cluster-based splits were retained to reduce sequence-similarity leakage.

### Methodology
1. Extracted amino-acid composition features from protein sequences using the 20 standard amino acids.
2. Trained a Logistic Regression classifier on the training split only.
3. Evaluated the trained model on the validation split.
4. Generated a validation confusion matrix and per-class classification report.
5. Kept the held-out test set unused.

### Validation results

| Metric | Score |
|---|---:|
| Accuracy | 0.3501 (35.01%) |
| Balanced accuracy | 0.3546 (35.46%) |
| Macro F1-score | 0.3427 (34.27%) |
| Weighted F1-score | 0.3386 (33.86%) |

### Initial observations
- The baseline provides a reference point for evaluating more expressive protein sequence representations.
- Performance varies considerably across enzyme classes.
- Class 7 achieved the strongest class-level F1-score (0.6095), while Class 4 had the weakest (0.1628).
- The confusion matrix shows substantial misclassification between several classes, including Classes 1 and 2.
- Amino-acid composition captures overall residue frequencies but does not represent the order of amino acids or detailed sequence context. This may limit its ability to distinguish enzyme functions.

### Limitations
- These are validation results, not final test results.
- The scores do not establish how well the model generalizes to unseen protein families beyond the current cluster-based split.
- No conclusion about the final model should be drawn from this baseline alone.

### Next steps
1. Preserve the current train, validation, and test splits.
2. Establish a stronger sequence representation, such as pretrained protein language-model embeddings.
3. Compare the new model against this baseline using the same validation split and metrics.
4. Select the model using validation results, then evaluate the chosen model on the held-out test set.

**Baseline output location:** `C:\protein_project_work\baseline`

**Test set:** Not evaluated.

## Baseline Experiment 2: 3-mer TF-IDF + Logistic Regression

**Status:** Completed — validation evaluation finished.

### Objective
Evaluate whether local amino-acid sequence patterns improve seven-class enzyme function prediction compared with amino-acid composition features.

### Dataset and methodology
- Training records: 1,757
- Validation records: 377
- Test set: Not used
- Feature representation: Character-level 3-mer TF-IDF
- Number of features: 8,000
- Classifier: Logistic Regression
- Class weighting: Balanced
- Random seed: 42

The TF-IDF vectorizer was fitted on the training sequences only. The fitted vectorizer transformed the validation sequences, and the classifier was trained using the training split. Evaluation was performed on the validation split.

### Validation results

| Metric | Score |
|---|---:|
| Accuracy | 41.91% |
| Balanced accuracy | 41.09% |
| Macro F1 | 42.07% |
| Weighted F1 | 42.63% |

### Comparison with Experiment 1

The amino-acid composition baseline achieved 35.01% accuracy, 35.46% balanced accuracy, 34.27% macro F1, and 33.86% weighted F1.

The 3-mer TF-IDF model improved accuracy by 6.90 percentage points and macro F1 by 7.80 percentage points on the validation set.

### Class-level observations

- Class 6 achieved the highest F1-score at approximately 0.560.
- Class 7 had precision of approximately 0.774 and recall of 0.500.
- Class 3 remained difficult to classify, with an F1-score of approximately 0.235.
- Confusion remained notable among classes 1, 2, 3 and 4.

### Interpretation and limitations

The results suggest that local amino-acid sequence patterns provide more useful predictive information than overall amino-acid composition alone for this dataset and evaluation setup.

However, validation performance remains moderate, and several classes are difficult to distinguish. These findings do not establish biological equivalence between classes or demonstrate generalization to unrelated protein families.

At the time of this earlier baseline, the held-out test set had not been used. Later, after model selection, the final ESM-2 model was evaluated on it; see the ESM-2 addendum at the end of this document.

### Next step

Run a 4-mer TF-IDF baseline using the same training and validation splits. Compare it against the amino-acid composition and 3-mer models before selecting a model for final test evaluation.



# Baseline Experiments: Comparative Evaluation

**Status:** Completed
**Dataset:** 2,673 quality-filtered enzyme protein sequences
**Task:** Seven-class enzyme function prediction
**Classes:** EC top-level classes 1–7

## 1. Experimental Setup

The dataset was divided into training, validation, and test sets using sequence-similarity clustering to reduce leakage between related proteins.

| Split | Records | Percentage |
|---|---:|---:|
| Training | 1,757 | 65.7% |
| Validation | 377 | 14.1% |
| Test | 539 | 20.2% |
| Total | 2,673 | 100% |

No cluster IDs were shared across the three splits. The held-out test set has not been used for model selection or evaluation.

The split proportions differ slightly from the originally intended proportions because complete sequence clusters were kept together.

## 2. Evaluation Metrics

- **Accuracy:** Proportion of all predictions that are correct.
- **Balanced accuracy:** Average recall across the seven enzyme classes.
- **Macro F1:** Unweighted mean of the F1 scores across classes; used as the primary model-selection metric.
- **Weighted F1:** F1 averaged according to the number of examples in each class.
- **Confusion matrix:** Used to inspect class-specific prediction errors.

## 3. Experiment Results

| Experiment | Accuracy | Balanced Accuracy | Macro F1 |
|---|---:|---:|---:|
| Amino-acid composition + Logistic Regression | 35.01% | 35.46% | 0.3427 |
| 3-mer TF-IDF + Logistic Regression | 41.91% | 41.09% | 0.4207 |
| 4-mer TF-IDF + Logistic Regression | 41.64% | 39.78% | 0.4007 |
| 3-mer TF-IDF + Linear SVM | 46.95% | 46.63% | 0.4691 |
| 5-mer TF-IDF + Linear SVM | 39.26% | 38.83% | 0.4036 |

*Note: Values are from validation experiments. The 3-mer SVM result shown is the tuned model with C = 1.0.*

## 4. Experiment 1: Amino-Acid Composition

**Method:** Calculate the relative frequency of each of the 20 standard amino acids and train a Logistic Regression classifier.

**Purpose:** Establish a simple baseline using global sequence composition without considering the order of amino acids.

**Validation accuracy:** 35.01%
**Validation macro F1:** 0.3427

This baseline provides a reference for evaluating whether local sequence patterns improve prediction.

## 5. Experiments 2–4: K-mer TF-IDF Features

K-mers are consecutive amino-acid substrings of length k. TF-IDF converts these patterns into weighted numerical features.

- A 3-mer contains three consecutive amino acids.
- A 4-mer contains four consecutive amino acids.
- A 5-mer contains five consecutive amino acids.

The experiments evaluated whether local sequence-pattern representations improve enzyme-class prediction.

### 3-mer TF-IDF + Logistic Regression

Validation accuracy was 41.91%, with a macro F1 of 0.4207.

### 4-mer TF-IDF + Logistic Regression

Validation accuracy was 41.64%, with a macro F1 of 0.4007.

The 4-mer model did not improve upon the 3-mer Logistic Regression model.

### 3-mer TF-IDF + Linear SVM

Validation accuracy reached 46.95%, and macro F1 reached 0.4691.

This was the strongest validation result among the baseline experiments completed so far.

## 6. Experiment 5: Linear SVM Hyperparameter Tuning

The regularization parameter C was tested using the validation macro F1 score for model selection.

| C | Accuracy | Balanced Accuracy | Macro F1 |
|---:|---:|---:|---:|
| 0.1 | 43.50% | 42.92% | 0.4274 |
| 0.5 | 46.15% | 45.92% | 0.4631 |
| 1.0 | 46.95% | 46.63% | 0.4691 |
| 2.0 | 44.03% | 43.92% | 0.4422 |
| 5.0 | 42.44% | 42.30% | 0.4274 |

**Selected parameter:** C = 1.0.

The model with C = 1.0 achieved the highest validation macro F1 among the tested values.

## 7. Experiment 6: 5-mer TF-IDF + Linear SVM

The 5-mer experiment used a maximum vocabulary size of 100,000 features, with C = 1.0 and balanced class weights.

| Metric | Result |
|---|---:|
| Training records | 1,757 |
| Validation records | 377 |
| Actual TF-IDF features | 100,000 |
| Accuracy | 39.26% |
| Balanced accuracy | 38.83% |
| Macro F1 | 0.4036 |
| Weighted F1 | 0.4004 |
| Training time | Approximately 0.20 seconds |

The 5-mer model performed worse than the 3-mer SVM model. The larger feature space and limited training data may contribute to this result, but these explanations have not yet been experimentally verified.

## 8. Error Analysis

The 3-mer SVM achieved its strongest F1 score for class 7, at 0.80. Class 3 was more difficult, with an F1 score of approximately 0.29.

The 5-mer SVM achieved high precision for class 7 (0.96) but lower recall (0.46), indicating that it missed many actual class-7 examples.

Confusion matrices and classification reports were saved for further analysis.

## 9. Current Conclusion

The 3-mer TF-IDF + Linear SVM model, with C = 1.0, is the strongest baseline evaluated so far. Its validation accuracy is 46.95%, with a macro F1 score of 0.4691.

These results suggest that local sequence-pattern features improve upon global amino-acid composition features for this dataset. However, prediction performance remains limited, and class-specific errors remain substantial.

The results do not establish that k-mer models are inferior to protein language models; that comparison has not yet been performed.

## 10. Next Experiment

The next planned experiment is to extract pretrained ESM-2 protein embeddings and train a lightweight classifier using the existing cluster-based training and validation splits.

The goal is to compare learned protein representations against the best k-mer baseline.

Model selection is complete for the ESM-2 + Linear SVM experiment, and the held-out test set has been evaluated once. Do not use these test results to tune the current model further; any new experiment should be clearly designated as a new study and should ideally use a fresh external test set for unbiased final evaluation.

---

# Addendum — ESM-2 Embeddings and Final Held-Out Test Evaluation

**Status:** Completed. This addendum reflects the later Kaggle notebook workflow and supersedes earlier statements in this document that ESM-2 was still planned or that the test set had not been evaluated.

## A. Environment and embedding extraction

- Environment: Kaggle notebook with GPU acceleration.
- GPU reported by the notebook: NVIDIA Tesla T4, approximately 15.64 GB GPU memory.
- PyTorch: `2.11.0+cu128`; CUDA available.
- Transformers: `5.16.1`.
- pandas: `2.3.3`.
- scikit-learn: `1.6.1`.
- Protein language model: `facebook/esm2_t12_35M_UR50D`.
- Embedding dimension: 480.
- Protein sequences longer than the model's practical token window were divided into windows; window-level embeddings were averaged to create one 480-dimensional vector per protein.
- Mean pooling excluded padding and special tokens.
- The train and validation embedding arrays had shapes `(1757, 480)` and `(377, 480)`. The test embedding array had shape `(539, 480)`.

The notebook initially encountered a test-embedding integrity problem: a diagnostic showed only one unique embedding row for 539 test proteins, and the classifier predicted class 4 for every test example. This was not accepted as a valid result. The test embeddings were regenerated using the same windowing and pooling workflow, after which the diagnostic showed 539 unique embedding rows. The final metrics below refer to the corrected embedding array, not the failed constant-prediction run.

## B. Validation comparison

The first downstream comparison used the existing cluster-aware training and validation splits. The held-out test set was not used to choose between these models.

| Model | Validation accuracy | Balanced accuracy | Macro F1 | Weighted F1 |
|---|---:|---:|---:|---:|
| ESM-2 embeddings + Logistic Regression | 0.6419 | 0.6342 | 0.6333 | 0.6373 |
| ESM-2 embeddings + initial Linear SVM | 0.6525 | 0.6443 | 0.6399 | 0.6441 |
| ESM-2 embeddings + tuned Linear SVM (`C = 2.0`) | 0.6605 | 0.6522 | 0.6463 | 0.6506 |

The tuned SVM used `class_weight="balanced"` and `random_state=42`. The selection metric was validation macro F1. Among the tested values of `C`, `C = 2.0` gave the best validation macro F1 (approximately 0.6463). The selected configuration was saved as `selected_esm2_svm_config.json`.

## C. Final training and held-out test results

After selecting the configuration, the downstream classifier was trained using the combined training and validation embeddings (2,134 proteins) and evaluated on the held-out test embeddings (539 proteins). The test set was not used during model selection.

| Metric | Held-out test result |
|---|---:|
| Accuracy | 0.7143 (71.43%) |
| Balanced accuracy | 0.7112 (71.12%) |
| Macro F1 | 0.7035 |
| Weighted F1 | 0.7093 |
| Records used to fit final classifier | 2,134 |
| Held-out test records | 539 |
| Selected SVM `C` | 2.0 |
| Test used for model selection? | No |

### Per-class test report

| EC top-level label | Precision | Recall | F1-score | Support |
|---:|---:|---:|---:|---:|
| 1 | 0.63 | 0.58 | 0.60 | 72 |
| 2 | 0.82 | 0.82 | 0.82 | 85 |
| 3 | 0.66 | 0.79 | 0.72 | 62 |
| 4 | 0.56 | 0.59 | 0.58 | 75 |
| 5 | 0.56 | 0.43 | 0.49 | 81 |
| 6 | 0.76 | 0.86 | 0.81 | 71 |
| 7 | 0.90 | 0.90 | 0.90 | 93 |
| **Macro average** | **0.70** | **0.71** | **0.70** | **539** |
| **Weighted average** | **0.71** | **0.71** | **0.71** | **539** |

### Interpretation

- Class 7 had the strongest test F1 (0.90), followed by class 2 (0.82) and class 6 (0.81).
- Class 5 was the most difficult class by F1 (0.49), with recall of 0.43; fewer than half of its test examples were recovered.
- The final test macro F1 of approximately 0.7035 is higher than the selected model's validation macro F1 of approximately 0.6463. This difference can occur because the finite test set may be easier or harder than the validation set; it should not be interpreted as proof that the model improved after test evaluation.
- These metrics describe performance on this dataset and this cluster-aware split. They do not prove experimental enzyme activity, full EC-number correctness, or performance on every unseen protein family.

## D. Saved outputs and version control

The notebook and selected results were copied into the Git repository and pushed to GitHub.

- Notebook: `protein_project/notebooks/esm2-enzyme-function-prediction.ipynb`
- Results directory: `protein_project/results/esm2_final/`
- Reported result files include:
  - `final_test_metrics.json`
  - `final_test_classification_report.csv`
  - `final_test_confusion_matrix.csv`
  - `final_test_confusion_matrix.png`
- The initial ESM-2 results commit was `458fb38` (`Add final ESM-2 enzyme classification results`).
- A later commit, `bcbf311` (`Add k-mer baseline and SVM experiments`), added k-mer baseline and SVM experiments.

The full embedding arrays were generated under `/kaggle/working/esm2_embeddings/`. Kaggle working-directory files are temporary unless downloaded or copied to persistent storage. The repository results listed above are the documented outputs; do not assume that all large `.npz` embedding arrays were committed.

## E. Current experiment status

| Component | Status |
|---|---|
| Dataset collection and quality control | Complete |
| MMseqs2 clustering and cluster-aware splits | Complete; integrity checks passed |
| Amino-acid composition + Logistic Regression | Complete |
| 3-mer / 4-mer / 5-mer baselines and Linear SVM comparisons | Complete, according to recorded validation results |
| ESM-2 embedding extraction | Complete for train, validation, and test |
| Validation-based SVM tuning | Complete; selected `C = 2.0` |
| Final held-out test evaluation | Complete once after model selection |
| Duplicate-flag audit | Still needs investigation |
| Reproducibility audit and final research write-up | Remaining work |

## F. Recommended next steps

1. **Freeze the current test result.** Do not keep tuning against the same test scores. Use validation data for any further development and clearly label any new experiment.
2. **Audit the eight inherited `is_exact_duplicate=True` flags** noted earlier, and record whether recalculation is needed.
3. **Check repository reproducibility.** Ensure the notebook records all required paths, model name, windowing parameters, pooling logic, classifier settings, random seed, and package versions.
4. **Preserve model artifacts intentionally.** Decide whether to save the final fitted classifier and any required scaler/configuration. Avoid committing large cached model weights or embedding arrays unless there is a clear reason.
5. **Perform error analysis.** Review the final confusion matrix and per-class precision/recall, particularly class 5, and distinguish observed confusions from hypotheses about their biological cause.
6. **Write the final methodology and limitations.** Explain the reviewed-UniProt data source, class-label definition, sequence filtering, MMseqs2 thresholds, cluster-aware split, ESM-2 model, pooling/windowing, validation-only model selection, and one-time test evaluation.
7. **For a stronger low-data claim, run a planned learning-curve experiment** using progressively smaller, cluster-aware subsets of the training data while keeping validation fixed. This is a new experiment; do not infer low-data performance from the full-training result alone.
8. **Consider an external or temporally separated test set** for a future study if the project needs a further unbiased final evaluation after additional experimentation.

## G. Updated handover summary

The current project contains 2,673 quality-checked protein sequences across seven top-level EC classes. MMseqs2 created 1,427 clusters and the cluster-aware partitions contain 1,757 training, 377 validation, and 539 test records with zero cluster overlap between splits.

The amino-acid-composition baseline and k-mer baselines provide progressively stronger conventional-feature references. ESM-2 (`facebook/esm2_t12_35M_UR50D`) produces 480-dimensional sequence embeddings. A Linear SVM with `C = 2.0`, selected using validation macro F1, was finally trained on train plus validation data and evaluated on the held-out test set. It achieved 71.43% accuracy, 71.12% balanced accuracy, 0.7035 macro F1, and 0.7093 weighted F1 on 539 test proteins.

**Next action:** finish the duplicate-flag audit, check reproducibility and artifact retention, complete error analysis, and write the research report. If continuing model development, use validation-based experiments and reserve a fresh external test set for a future final evaluation.


---

# Addendum — 11 October 2026: Low-data experiments and follow-up audits

This addendum supersedes earlier status statements in the chronological log that described the k-mer/ESM-2 experiments as planned or the held-out test set as unevaluated. Historical notes are retained as a record of the project timeline; the current status is summarized here and in the updated README.

## 1. Protein-level low-data learning curve

Using the existing 480-dimensional ESM-2 embeddings, the fixed 377-protein validation split, fractions of 10%, 25%, 50%, and 100% of the 1,757-protein training split, and seeds 42, 123, and 2026, both a Linear SVM (`C=2.0`, balanced class weights) and logistic regression were evaluated. Sampling was stratified by top-level EC class.

| Target fraction | Training proteins | SVM macro F1 (mean) | SVM macro F1 (SD) | Logistic regression macro F1 (mean) | Logistic regression macro F1 (SD) |
|---:|---:|---:|---:|---:|---:|
| 10% | 173 | 0.5281 | 0.0297 | 0.5336 | 0.0245 |
| 25% | 436 | 0.5642 | 0.0278 | 0.5912 | 0.0426 |
| 50% | 877 | 0.5681 | 0.0117 | 0.5878 | 0.0191 |
| 100% | 1,757 | 0.6463 | 0.0000 | 0.6333 | 0.0000 |

Other mean validation metrics from the recorded summary:

| Fraction | Model | Balanced accuracy | Accuracy | Weighted F1 |
|---:|---|---:|---:|---:|
| 10% | Linear SVM | 0.5369 | 0.5402 | 0.5313 |
| 10% | Logistic regression | 0.5319 | 0.5367 | 0.5370 |
| 25% | Linear SVM | 0.5651 | 0.5694 | 0.5672 |
| 25% | Logistic regression | 0.5910 | 0.5968 | 0.5965 |
| 50% | Linear SVM | 0.5681 | 0.5774 | 0.5741 |
| 50% | Logistic regression | 0.5879 | 0.5959 | 0.5931 |
| 100% | Linear SVM | 0.6522 | 0.6605 | 0.6506 |
| 100% | Logistic regression | 0.6342 | 0.6419 | 0.6373 |

Files were created in Kaggle under `/kaggle/working/esm2_results/low_data/`: `low_data_runs.csv`, `low_data_summary.csv`, and `low_data_learning_curve.png`. These outputs need to be downloaded or otherwise persisted before they can be relied on as repository artifacts.

## 2. Cluster-aware low-data learning curve

The training CSV contained a `cluster_id` column and matched embedding entry IDs exactly: 1,757 training entries and 377 validation entries were aligned with zero ID mismatches. The training split contains 914 distinct clusters, no missing cluster IDs, a largest cluster of 37 proteins, and one cluster containing more than one label.

For each subset fraction, complete clusters were sampled, retaining all proteins within each selected cluster. Because cluster sizes vary, actual protein counts differ from the target fraction. Results below are means across seeds 42, 123, and 2026.

| Target fraction | Mean proteins selected | Mean clusters selected | SVM macro F1 (mean) | SVM macro F1 (SD) | Logistic regression macro F1 (mean) | Logistic regression macro F1 (SD) |
|---:|---:|---:|---:|---:|---:|---:|
| 10% | 187.3 | 103.0 | 0.4675 | 0.0268 | 0.4336 | 0.0678 |
| 25% | 440.0 | 234.3 | 0.5753 | 0.0446 | 0.5743 | 0.0604 |
| 50% | 879.7 | 466.0 | 0.5825 | 0.0303 | 0.5982 | 0.0180 |
| 100% | 1,757 | 914 | 0.6463 | 0.0000 | 0.6333 | 0.0000 |

| Target fraction | Model | Balanced accuracy | Accuracy | Weighted F1 |
|---:|---|---:|---:|---:|
| 10% | Linear SVM | 0.4735 | 0.4783 | 0.4716 |
| 10% | Logistic regression | 0.4518 | 0.4607 | 0.4415 |
| 25% | Linear SVM | 0.5780 | 0.5862 | 0.5802 |
| 25% | Logistic regression | 0.5765 | 0.5818 | 0.5779 |
| 50% | Linear SVM | 0.5820 | 0.5880 | 0.5843 |
| 50% | Logistic regression | 0.5967 | 0.6021 | 0.6000 |
| 100% | Linear SVM | 0.6522 | 0.6605 | 0.6506 |
| 100% | Logistic regression | 0.6342 | 0.6419 | 0.6373 |

Files were created in Kaggle under `/kaggle/working/esm2_results/cluster_low_data/`: `cluster_low_data_runs.csv`, `cluster_low_data_summary.csv`, and `cluster_low_data_learning_curve.png`.

Interpretation: validation macro F1 improves overall as the training fraction grows, but the curve is not strictly monotonic at every intermediate fraction. Cluster-level subsets are the more relevant low-data setting for this project, because related proteins are kept together during sampling. However, these curves use the same fixed validation split and do not establish generalization to entirely unseen protein families. The 100% result repeats across seeds because all training examples are used.

## 3. Duplicate-flag audit

A read-only audit of `protein_project/processed/uniprot_expanded_model_ready.csv` reported 2,673 rows, 2,673 unique exact-sequence hashes, 2,673 unique raw `Sequence` values, and 2,673 unique `sequence_clean` values. Eight rows nevertheless had `is_exact_duplicate=True`; recomputing exact duplicates from the current cleaned sequences found zero duplicated rows. The eight accession IDs shown by the audit were P0D176, Q15185, P68106, P63038, Q9DBJ1, Q6AYT3, P06686, and P06687.

This points to an inherited/stale flag or a different historical definition of the flag, but its origin has not been proven. The original model-ready CSV was not changed, and the eight records were not removed. Audit outputs were saved locally at `protein_project/results/data_audit/previously_flagged_duplicates.csv` and `protein_project/results/data_audit/duplicate_flag_audit.csv`.

## 4. External UniProt pool attempt

UniProt access was confirmed, but the external data-collection attempt did not yield a usable single-label evaluation set. After excluding the original 2,673 accessions, one bounded collection produced 179 unique proteins. Applying the parser used for the primary single-label external evaluation left only 2 proteins with one top-level EC class; 177 had annotations spanning multiple top-level classes. The two remaining proteins were all class 1, which is far too small and imbalanced for external performance metrics. No external accuracy, macro F1, or balanced accuracy was calculated.

The intermediate pool and audit files were created in Kaggle `/kaggle/working/enzyme_candidates/`, but these are temporary unless downloaded or persisted. Do not describe this attempt as a successful external benchmark. Multi-EC annotations should not be silently collapsed to the class used in the retrieval query.

## 5. Current completion status and next actions

**Completed:** dataset preparation, MMseqs2 clustering and split checks; amino-acid-composition and k-mer baselines; ESM-2 embedding generation; validation-based SVM tuning; final held-out test evaluation; protein-level low-data experiment; cluster-aware low-data experiment; duplicate-flag audit.

**Not completed:** a valid external performance evaluation; a validated novel-candidate discovery system; final consolidated report and reproducibility package.

Recommended order to finish efficiently:
1. Download/persist the two low-data CSV summaries and plots if they are to be included in the repository.
2. Add the low-data tables and figures to the README/report and clearly label them as validation results.
3. Preserve the final test metrics as frozen; do not tune against the already-evaluated test set.
4. Decide whether a minimal candidate-ranking prototype is feasible using a separate suitable candidate set. Label confidence and sequence similarity as computational prioritization signals, not proof of novel function.
5. Complete error analysis, methodology, limitations, and final GitHub cleanup.

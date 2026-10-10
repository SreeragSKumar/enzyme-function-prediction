# Complete Project Progress and Experiment Log
## Low-Data Protein Language Model for Enzyme Function Prediction and Novel Candidate Discovery

**Last updated:** 9 October 2026  
**Status:** Dataset preparation, MMseqs2 clustering, cluster-aware splitting, and the first amino-acid-composition baseline are complete. The k-mer and ESM-2 experiments are planned but have not been run. The test split has not been evaluated.

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

Steps 1–4 are complete; the first baseline's validation results are recorded below.

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

**Next action:** build and evaluate the k-mer baseline in a new script, preserving the current baseline and all split files.

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

The held-out test set has not been used. No conclusions about test performance should be made yet.

### Next step

Run a 4-mer TF-IDF baseline using the same training and validation splits. Compare it against the amino-acid composition and 3-mer models before selecting a model for final test evaluation.
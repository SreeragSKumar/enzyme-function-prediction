
import pandas as pd
from pathlib import Path

file_path = (
    Path.cwd()
    / "protein_project"
    / "uniprot_model_ready_pilot.csv"
)

df = pd.read_csv(file_path)

print("Total model-ready proteins:", len(df))
print("Unique sequences:", df["sequence"].nunique())

print("\nExact class counts:")
print(df["label"].value_counts().sort_index())

print("\nClass proportions (%):")
print(
    (df["label"].value_counts(normalize=True).sort_index() * 100)
    .round(2)
)

print("\nSequence-length statistics:")
print(df["sequence_length"].describe())

print("\nMissing labels:", df["label"].isna().sum())
print("Missing sequences:", df["sequence"].isna().sum())

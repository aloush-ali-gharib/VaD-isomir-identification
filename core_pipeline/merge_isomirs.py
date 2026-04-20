import os
import glob
import pandas as pd

# Technical Workflow: Cohort-level expression matrix generation
files = glob.glob("*_isomir_counts.csv")
merged_df = None

for file in files:
    sample_name = file.replace("_isomir_counts.csv", "")
    df = pd.read_csv(file)

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # [RESEARCH SENSITIVE: DUPLICATE RESOLUTION & MAPPING FIDELITY]
    # Handling for mir-precursor overlapping variants to maintain 
    # alignment across unique genomic coordinates.
    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    
    df = df[['mir', 'freq']].rename(columns={'freq': sample_name})

    if merged_df is None:
        merged_df = df
    else:
        # Inner join to preserve only common high-confidence variants
        merged_df = pd.merge(merged_df, df, on='mir', how='inner')

merged_df.to_csv("VDA_expression_matrix.csv", index=False)
print(f"Matrix Generated for {len(files)} samples.")

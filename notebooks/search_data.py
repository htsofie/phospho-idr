#%% [markdown]
## Notebooks setup and data config

import pandas as pd
import sys
import yaml
from pathlib import Path
from pandas._config import config

# set project root (phospho_root)
PROJECT_ROOT = Path(__file__).resolve().parents[1] if "__file__" in globals() else Path.cwd().parents[0]

#Pick config manually here: (currently rat)
config_path = PROJECT_ROOT / "configs/rat.yaml"

with open(config_path) as f:
    config = yaml.safe_load(f)

print(config)

#%%
def search_uniprot_id(file_path, target_id):
    """
    Read Excel file and search for a specific UniProt ID
    
    Args:
        file_path (str): Path to the Excel file
        target_id (str): UniProt ID to search for
    
    Returns:
        bool: True if found, False otherwise
    """
    try:
        # Read the Excel file
        print(f"Reading Excel file: {file_path}")
        df = pd.read_excel(file_path)
        
        # Display basic info about the data
        print(f"\nData shape: {df.shape}")
        print(f"Columns: {list(df.columns)}")
        print(f"\nFirst few rows:")
        print(df.head())
        
        # Search for the target UniProt ID in all columns
        print(f"\nSearching for UniProt ID: {target_id}")
        found = False
        matches = []
        
        for col in df.columns:
            # Check if any cell in this column contains the target ID
            mask = df[col].astype(str).str.contains(target_id, case=False, na=False)
            if mask.any():
                found = True
                matching_rows = df[mask]
                for idx, row in matching_rows.iterrows():
                    matches.append({
                        'column': col,
                        'row_index': idx,
                        'value': row[col]
                    })
        
        if found:
            print(f"\n✓ Found {len(matches)} match(es) for UniProt ID {target_id}:")
            for match in matches:
                print(f"  - Column '{match['column']}', Row {match['row_index']}: {match['value']}")
        else:
            print(f"\n✗ No matches found for UniProt ID {target_id}")
            
        return found
        
    except Exception as e:
        print(f"Error reading file: {e}")
        return False

if __name__ == "__main__":
    # File path and target UniProt ID
    file_path = config["raw_data"]
    target_uniprot_id = "Q9R1N3"
    
    print("=" * 50)
    print("UniProt ID Search Tool")
    print("=" * 50)
    
    # Debug information
    import os
    print(f"Python executable: {sys.executable}")
    print(f"Python version: {sys.version}")
    print(f"Current working directory: {os.getcwd()}")
    print(f"File exists: {os.path.exists(file_path)}")
    print(f"File path (absolute): {os.path.abspath(file_path)}")
    print("=" * 50)
    
    # Search for the UniProt ID
    result = search_uniprot_id(file_path, target_uniprot_id)
    
    print("\n" + "=" * 50)
    if result:
        print(f"RESULT: UniProt ID {target_uniprot_id} was FOUND in the data")
    else:
        print(f"RESULT: UniProt ID {target_uniprot_id} was NOT FOUND in the data")
    print("=" * 50)

# %%
# check # of missing value in both Uniprot and ENSEMBL
df = pd.read_excel(config["raw_data"])
missing_count = df[['Uniprot', 'ENSEMBL']].isna().all(axis=1).sum()
print(f"Missing count: {missing_count}")
print(f"Total count: {len(df)}")
print(f"Percentage of missing values: {missing_count / len(df) * 100:.2f}%")

# %% [markdown]
# To look at range of values of localized prob to determine if all sites in dataset are localized or if there are some that need to be removed

df = pd.read_excel(config["raw_data"])
col = "Localization Prob"
min_val = df[col].min()
max_val = df[col].max()
print("Minimum:", min_val)
print("Maximum:", max_val)
# Min is 0.745606 and Max is 1.0
# Some scores are less than 0.75 so should be removed
# To determine number of scores less than 0.75
num_rows = (df[col] < 0.75).sum()
print("Number of rows with localization_prob < 0.75:", num_rows)
# %% [markdown]
# To look at the values of these 11 rows w localization prob less than 0.75:
df.loc[df[col] < 0.75, col]

# %%
# Count rows missing BOTH UniProt AND ENSEMBL values
df = pd.read_excel(config["raw_data"])
missing_both = df[['Uniprot', 'ENSEMBL']].isna().all(axis=1).sum()
total_rows = len(df)

print(f"Rows missing BOTH UniProt AND ENSEMBL: {missing_both}")
print(f"Total rows: {total_rows}")
print(f"Percentage missing both: {missing_both / total_rows * 100:.2f}%")

# %%
# Count unique phosphoproteins in the Protein column
df = pd.read_excel(config["raw_data"])
unique_proteins = df['Protein'].nunique()
total_rows = len(df)

print(f"Number of unique phosphoproteins: {unique_proteins}")
print(f"Total phosphorylation sites: {total_rows}")
print(f"Average sites per protein: {total_rows / unique_proteins:.2f}")

# %%
# Calculate percentage of dataset with aligned phosphosites
# Adjust the file path below as needed
import pandas as pd

aligned_data_path = "data/processed/rat/full_total_blast_aligned.csv"

aligned_df = pd.read_csv(aligned_data_path)

# Count total number of phosphosites
total_phosphosites = len(aligned_df)

# Count phosphosites with successful alignment
aligned_phosphosites = aligned_df['alignment_success'].sum()

# Calculate percentage
alignment_percentage = (aligned_phosphosites / total_phosphosites) * 100

print(f"Total phosphosites in dataset: {total_phosphosites:,}")
print(f"Successfully aligned phosphosites: {aligned_phosphosites:,}")
print(f"Percentage of dataset with aligned phosphosites: {alignment_percentage:.2f}%")
print(f"Percentage of dataset with failed alignment: {100 - alignment_percentage:.2f}%")

# %%

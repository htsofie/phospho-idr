# %% [markdown]
# # Clean Data - Extract Specific Columns
# %%
# read the parquet file
import pandas as pd
import re

#%%
# Read the parquet file
df = pd.read_parquet("Phosphomouse_phosphorylation_sites.parquet")

#%%
print(df.columns)
# %% [markdown]
# # Extract only the required columns and clean the Protein IPI column
# %%

# Function to extract SWISS-PROT identifier from Protein IPI column
def extract_swiss_prot(protein_ipi):
    if pd.isna(protein_ipi):
        return None
    # Look for pattern SWISS-PROT:XXXXX| or SWISS-PROT:XXXXX
    match = re.search(r'SWISS-PROT:([^|]+)', str(protein_ipi))
    if match:
        return match.group(1)
    return None

# Select only the required columns (including Localized for filtering)
required_columns = [
    'Protein IPI (Assigned by Sequest)',  # Will be cleaned
    'Gene Symbol',
    'Protein Description', 
    'Sequence',
    'Residue',
    'Site',
    'Localized'  # Needed for filtering
]

# Create new dataframe with only required columns
df_clean = df[required_columns].copy()

# Filter out rows where Localized is 'Amb' (ambiguous)
print(f"Original data shape: {df_clean.shape}")
df_clean = df_clean[df_clean['Localized'] != 'Amb']
print(f"After removing 'Amb' rows: {df_clean.shape}")

# Remove the Localized column as it's no longer needed
df_clean = df_clean.drop(columns=['Localized'])

# Clean the Protein IPI column to extract only SWISS-PROT identifier
df_clean['Protein Identifier (SWISS-PROT)'] = df_clean['Protein IPI (Assigned by Sequest)'].apply(extract_swiss_prot)

# Drop the original Protein IPI column and rename the new one
df_clean = df_clean.drop(columns=['Protein IPI (Assigned by Sequest)'])

# Reorder columns to match your requested order
final_columns = [
    'Protein Identifier (SWISS-PROT)',
    'Gene Symbol',
    'Protein Description',
    'Sequence', 
    'Residue',
    'Site'
]

df_clean = df_clean[final_columns]

#%%
# Print the cleaned data info
print("Cleaned dataframe shape:", df_clean.shape)
print("Cleaned columns:", df_clean.columns.tolist())
print("\nFirst few rows of cleaned data:")
print(df_clean.head())

#%%
# Save the cleaned data to a new parquet file
output_file = "cleaned_specific_columns.parquet"
df_clean.to_parquet(output_file, index=False)
print(f"\nCleaned data saved to: {output_file}")

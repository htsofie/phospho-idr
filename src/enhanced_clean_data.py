# Enhanced Data Cleaning Script
# Extracts all identifier types from Protein IPI column and filters for localized sites only

import pandas as pd
import re

def extract_identifiers(protein_ipi):
    """
    Extract all identifier types from the Protein IPI column.
    Returns a dictionary with all found identifiers.
    """
    if pd.isna(protein_ipi):
        return {
            'swiss_prot': None,
            'ensembl': None, 
            'refseq': None,
            'trembl': None
        }
    
    protein_str = str(protein_ipi)
    identifiers = {
        'swiss_prot': None,
        'ensembl': None,
        'refseq': None, 
        'trembl': None
    }
    
    # Extract SWISS-PROT identifier
    swiss_match = re.search(r'SWISS-PROT:([^|]+)', protein_str)
    if swiss_match:
        identifiers['swiss_prot'] = swiss_match.group(1)
    
    # Extract Ensembl identifier (note: it's ENSEMBL in uppercase)
    ensembl_match = re.search(r'ENSEMBL:([^|]+)', protein_str)
    if ensembl_match:
        # Take the first Ensembl ID if there are multiple (separated by semicolon)
        ensembl_id = ensembl_match.group(1).split(';')[0]
        identifiers['ensembl'] = ensembl_id
    
    # Extract RefSeq identifier (if present)
    refseq_match = re.search(r'RefSeq:([^|]+)', protein_str)
    if refseq_match:
        identifiers['refseq'] = refseq_match.group(1)
    
    # Extract TREMBL identifier
    trembl_match = re.search(r'TREMBL:([^|]+)', protein_str)
    if trembl_match:
        identifiers['trembl'] = trembl_match.group(1)
    
    return identifiers

# Read the original parquet file
print("Reading original data...")
df = pd.read_parquet("data/processed/Phosphomouse_phosphorylation_sites.parquet")
print(f"Original data shape: {df.shape}")

# Select required columns
required_columns = [
    'Protein IPI (Assigned by Sequest)',
    'Gene Symbol',
    'Protein Description', 
    'Sequence',
    'Residue',
    'Site',
    'Localized'
]

df_clean = df[required_columns].copy()

# Check Localized column distribution
print("\nLocalized column distribution:")
print(df_clean['Localized'].value_counts())

# Filter to keep only localized sites (remove 'Amb' ambiguous sites)
print(f"\nOriginal data shape: {df_clean.shape}")
df_clean = df_clean[df_clean['Localized'] != 'Amb']
print(f"After removing 'Amb' rows: {df_clean.shape}")

# Extract all identifier types
print("\nExtracting identifier types...")
identifier_data = df_clean['Protein IPI (Assigned by Sequest)'].apply(extract_identifiers)

# Add identifier columns
df_clean['Protein Identifier (SWISS-PROT)'] = [x['swiss_prot'] for x in identifier_data]
df_clean['Protein Identifier (Ensembl)'] = [x['ensembl'] for x in identifier_data]
df_clean['Protein Identifier (RefSeq)'] = [x['refseq'] for x in identifier_data]
df_clean['Protein Identifier (TREMBL)'] = [x['trembl'] for x in identifier_data]

# Remove the original Protein IPI column and Localized column
df_clean = df_clean.drop(columns=['Protein IPI (Assigned by Sequest)', 'Localized'])

# Reorder columns
final_columns = [
    'Gene Symbol',
    'Protein Description',
    'Sequence', 
    'Residue',
    'Site',
    'Protein Identifier (SWISS-PROT)',
    'Protein Identifier (Ensembl)',
    'Protein Identifier (RefSeq)',
    'Protein Identifier (TREMBL)'
]

df_clean = df_clean[final_columns]

# Print statistics
print(f"\nFinal cleaned data shape: {df_clean.shape}")
print("\nIdentifier coverage:")
for col in ['Protein Identifier (SWISS-PROT)', 'Protein Identifier (Ensembl)', 
            'Protein Identifier (RefSeq)', 'Protein Identifier (TREMBL)']:
    count = df_clean[col].notna().sum()
    percentage = (count / len(df_clean)) * 100
    print(f"  {col}: {count:,} sites ({percentage:.1f}%)")

# Save the enhanced cleaned data
output_file = "data/processed/cleaned_phosphomouse_site_data.parquet"
df_clean.to_parquet(output_file, index=False)
print(f"\nEnhanced cleaned data saved to: {output_file}")

print("\nFirst few rows of cleaned data:")
print(df_clean.head())

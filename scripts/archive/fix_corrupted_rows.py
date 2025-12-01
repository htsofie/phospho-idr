#!/usr/bin/env python3
"""
Fix corrupted rows in aligned_interfaces.csv where protein_id contains sequence data.

This script identifies rows where protein_id contains sequence data instead of an ID,
and attempts to fix them by reconstructing the correct data structure.
"""

import pandas as pd
import argparse
import logging
from pathlib import Path

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def identify_corrupted_rows(df):
    """Identify rows where protein_id looks like sequence data or IRES data."""
    # Protein IDs are typically 6-10 characters, alphanumeric
    # Anything longer than 20 chars is suspicious
    # And if it's all amino acids, it's definitely sequence data
    # If it starts with [, it's IRES data
    amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    
    def is_corrupted(protein_id):
        if pd.isna(protein_id):
            return False  # Missing protein_id is handled separately
        protein_id_str = str(protein_id).strip()
        if protein_id_str == '' or protein_id_str == 'nan':
            return False  # Missing protein_id is handled separately
        
        # Check if it looks like IRES data (starts with [)
        if protein_id_str.startswith('[') and ']' in protein_id_str:
            return True
        
        if len(protein_id_str) <= 20:
            return False
        # Check if it looks like sequence data (all amino acids in first 20 chars)
        if all(c in amino_acids for c in protein_id_str[:20]):
            return True
        return False
    
    corrupted = df[df['protein_id'].astype(str).apply(is_corrupted)]
    return corrupted.index.tolist()


def try_to_fix_row(corrupted_row, df_all, source_df=None):
    """
    Try to reconstruct the correct data for a corrupted row.
    
    Strategy:
    1. The corrupted protein_id is actually sequence data (1590 chars)
    2. Search source file (human_ortholog_mapped.csv) for this sequence
    3. Use that row's protein_id and data as the correct values
    4. Preserve any alignment data from the corrupted row
    """
    corrupted_protein_id = str(corrupted_row['protein_id'])
    corrupted_seq_start = corrupted_protein_id[:100]
    
    # First try to find in the current dataframe
    matching_rows = df_all[df_all['full_sequence'].astype(str).str.startswith(
        corrupted_seq_start, na=False
    )]
    
    # If not found, try source file
    if len(matching_rows) == 0 and source_df is not None:
        matching_rows = source_df[source_df['full_sequence'].astype(str).str.startswith(
            corrupted_seq_start, na=False
        )]
        if len(matching_rows) > 0:
            logger.info(f"Found match in source file")
            # Convert source row to match current dataframe structure
            source_row = matching_rows.iloc[0]
            fixed_row = pd.Series(dtype=object)
            
            # Copy base columns from source
            base_cols = ['protein_id', 'full_sequence', 'human_ortholog_id', 
                        'human_ortholog_sequence']
            for col in base_cols:
                if col in source_row.index:
                    fixed_row[col] = source_row[col]
            
            # Initialize other columns as None
            for col in df_all.columns:
                if col not in fixed_row.index:
                    fixed_row[col] = None
            
            # Preserve alignment data from corrupted row if it exists
            alignment_cols = ['aligned_IRES', 'alignment_percent', 'matched_IRES', 
                            'mismatched_IRES', 'alignment_error', 'total_alignment',
                            'interface_alignment_percent', 'matched_IRES_count',
                            'human_phosphoprotein_IRES_count']
            for col in alignment_cols:
                if col in corrupted_row.index:
                    corrupted_val = corrupted_row[col]
                    if pd.notna(corrupted_val) and corrupted_val != '':
                        fixed_row[col] = corrupted_val
            
            logger.info(f"Reconstructed row with protein_id: {fixed_row.get('protein_id', 'N/A')}")
            return fixed_row
    
    if len(matching_rows) == 0:
        logger.warning(f"Could not find matching row for corrupted sequence")
        return None
    
    # Use the first match from current dataframe
    correct_row = matching_rows.iloc[0]
    correct_protein_id = correct_row['protein_id']
    
    logger.info(f"Found match: corrupted sequence matches protein_id '{correct_protein_id}'")
    
    # Create fixed row by copying the correct row's structure
    fixed_row = correct_row.copy()
    
    # Preserve alignment columns from corrupted row if they exist and are not None
    alignment_cols = ['aligned_IRES', 'alignment_percent', 'matched_IRES', 
                      'mismatched_IRES', 'alignment_error', 'total_alignment',
                      'interface_alignment_percent', 'matched_IRES_count',
                      'human_phosphoprotein_IRES_count']
    for col in alignment_cols:
        if col in corrupted_row.index and col in fixed_row.index:
            corrupted_val = corrupted_row[col]
            if pd.notna(corrupted_val) and corrupted_val != '':
                fixed_row[col] = corrupted_val
    
    return fixed_row


def fix_corrupted_csv(input_file, output_file=None, source_file=None):
    """Fix corrupted rows in the CSV file."""
    if output_file is None:
        output_file = input_file.replace('.csv', '_fixed.csv')
    
    logger.info(f"Loading file: {input_file}")
    df = pd.read_csv(input_file)
    logger.info(f"Loaded {len(df)} rows")
    
    # Try to load source file if provided
    source_df = None
    if source_file and Path(source_file).exists():
        logger.info(f"Loading source file: {source_file}")
        source_df = pd.read_csv(source_file)
        logger.info(f"Loaded {len(source_df)} rows from source")
    elif not source_file:
        # Try to guess source file location
        if 'aligned_interfaces' in input_file:
            # Try human_ortholog_mapped.csv in same directory
            source_file = input_file.replace('aligned_interfaces.csv', 'human_ortholog_mapped.csv')
            if Path(source_file).exists():
                logger.info(f"Auto-detected source file: {source_file}")
                source_df = pd.read_csv(source_file)
                logger.info(f"Loaded {len(source_df)} rows from source")
    
    # Identify corrupted rows
    corrupted_indices = identify_corrupted_rows(df)
    logger.info(f"Found {len(corrupted_indices)} corrupted rows: {corrupted_indices}")
    
    if len(corrupted_indices) == 0:
        logger.info("No corrupted rows found. File is clean!")
        return
    
    # Try to fix each corrupted row
    fixed_count = 0
    rows_to_remove = []
    
    for idx in corrupted_indices:
        logger.info(f"\nFixing row {idx}...")
        corrupted_row = df.iloc[idx]
        
        fixed_row = try_to_fix_row(corrupted_row, df, source_df)
        
        if fixed_row is not None:
            # Update the dataframe
            for col in df.columns:
                if col in fixed_row.index:
                    df.at[idx, col] = fixed_row[col]
            fixed_count += 1
            logger.info(f"  ✓ Fixed row {idx}")
        else:
            logger.warning(f"  ✗ Could not fix row {idx} - will be removed")
            rows_to_remove.append(idx)
    
    # Remove rows that couldn't be fixed
    if rows_to_remove:
        logger.info(f"\nRemoving {len(rows_to_remove)} corrupted rows that couldn't be fixed: {rows_to_remove}")
        df = df.drop(index=rows_to_remove).reset_index(drop=True)
        logger.info(f"After removal: {len(df)} rows remaining")
    
    # Also check for missing protein_ids
    missing_protein_id_indices = []
    for idx, row in df.iterrows():
        protein_id = str(row.get('protein_id', '')).strip()
        if pd.isna(row.get('protein_id')) or protein_id == '' or protein_id.lower() == 'nan':
            missing_protein_id_indices.append(idx)
    
    if missing_protein_id_indices:
        logger.info(f"\nFound {len(missing_protein_id_indices)} rows with missing protein_id: {missing_protein_id_indices}")
        df = df.drop(index=missing_protein_id_indices).reset_index(drop=True)
        logger.info(f"After removal: {len(df)} rows remaining")
    
    logger.info(f"\nFixed {fixed_count} out of {len(corrupted_indices)} corrupted rows")
    logger.info(f"Removed {len(rows_to_remove)} rows that couldn't be fixed")
    if missing_protein_id_indices:
        logger.info(f"Removed {len(missing_protein_id_indices)} rows with missing protein_id")
    
    # Save fixed file
    logger.info(f"Saving fixed file to: {output_file}")
    df.to_csv(output_file, index=False)
    logger.info(f"✓ Saved {len(df)} rows to {output_file}")
    
    return output_file


def main():
    parser = argparse.ArgumentParser(
        description='Fix corrupted rows in aligned_interfaces.csv'
    )
    parser.add_argument(
        '--input',
        type=str,
        default='data/processed/mouse/aligned_interfaces.csv',
        help='Input CSV file to fix'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Output CSV file (default: input_fixed.csv)'
    )
    parser.add_argument(
        '--source',
        type=str,
        default=None,
        help='Source file to use for reconstruction (e.g., human_ortholog_mapped.csv)'
    )
    
    args = parser.parse_args()
    
    fix_corrupted_csv(args.input, args.output, args.source)


if __name__ == '__main__':
    main()


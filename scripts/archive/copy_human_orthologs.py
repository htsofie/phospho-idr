#!/usr/bin/env python3
"""
Copy human ortholog data from existing mapped files to a new full_disorder.csv file.

This script searches for protein_id matches in existing human_ortholog_mapped.csv files
and copies the human ortholog data to avoid unnecessary API calls.

Usage:
    python scripts/archive/copy_human_orthologs.py \
        --input data/processed/mouse/full_disorder.csv \
        --mapped-file1 data/processed/mouse/human_ortholog_mapped.csv \
        --mapped-file2 data/processed/mouse/human_ortholog_mapped1.csv \
        --output data/processed/mouse/full_disorder_with_orthologs.csv
"""

import pandas as pd
import argparse
import logging
import os
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Columns to copy from mapped files
ORTHOLOG_COLUMNS = [
    'human_ortholog_id',
    'human_ortholog_description',
    'human_ortholog_sequence',
    'api_error'
]


def load_mapped_file(filepath: str) -> pd.DataFrame:
    """
    Load a human ortholog mapped file.
    
    Args:
        filepath: Path to mapped CSV file
        
    Returns:
        DataFrame with protein_id and ortholog columns, or None if error
    """
    if not os.path.exists(filepath):
        logger.warning(f"File not found: {filepath}")
        return None
    
    try:
        df = pd.read_csv(filepath)
        logger.info(f"Loaded {len(df)} rows from {filepath}")
        
        # Check if protein_id column exists
        if 'protein_id' not in df.columns:
            logger.warning(f"  Warning: 'protein_id' column not found in {filepath}")
            return None
        
        # Check which ortholog columns are available
        available_cols = [col for col in ORTHOLOG_COLUMNS if col in df.columns]
        logger.info(f"  Available ortholog columns: {available_cols}")
        
        return df
    except Exception as e:
        logger.error(f"Error loading {filepath}: {e}")
        return None


def copy_ortholog_data(input_file: str, mapped_files: list, output_file: str) -> None:
    """
    Copy human ortholog data from mapped files to input file.
    
    Args:
        input_file: Path to new full_disorder.csv file
        mapped_files: List of paths to existing human_ortholog_mapped.csv files
        output_file: Path to output CSV file
    """
    logger.info("=" * 80)
    logger.info("Copy Human Ortholog Data")
    logger.info("=" * 80)
    logger.info(f"Input file: {input_file}")
    logger.info(f"Mapped files: {mapped_files}")
    logger.info(f"Output file: {output_file}")
    logger.info("")
    
    # Load input file
    logger.info(f"Loading input file: {input_file}")
    try:
        df_input = pd.read_csv(input_file)
        logger.info(f"Loaded {len(df_input)} rows")
    except Exception as e:
        logger.error(f"Error reading input file: {e}")
        return
    
    # Check required column
    if 'protein_id' not in df_input.columns:
        logger.error("Error: 'protein_id' column not found in input file")
        logger.error(f"Available columns: {list(df_input.columns)}")
        return
    
    # Initialize ortholog columns in input dataframe
    for col in ORTHOLOG_COLUMNS:
        if col not in df_input.columns:
            df_input[col] = None
    
    # Load all mapped files and combine them
    logger.info("")
    logger.info("Loading mapped files...")
    all_mapped_data = []
    
    for mapped_file in mapped_files:
        df_mapped = load_mapped_file(mapped_file)
        if df_mapped is not None:
            # Select only protein_id and ortholog columns
            cols_to_use = ['protein_id'] + [col for col in ORTHOLOG_COLUMNS if col in df_mapped.columns]
            df_subset = df_mapped[cols_to_use].copy()
            
            # Group by protein_id and take first non-null value for each column
            # This handles cases where same protein_id appears multiple times
            grouped = df_subset.groupby('protein_id').first().reset_index()
            all_mapped_data.append(grouped)
            logger.info(f"  Extracted {len(grouped)} unique protein_ids from {mapped_file}")
    
    if not all_mapped_data:
        logger.error("No valid mapped files could be loaded")
        return
    
    # Combine all mapped data (later files take precedence if same protein_id exists)
    logger.info("")
    logger.info("Combining mapped data...")
    df_combined = pd.concat(all_mapped_data, ignore_index=True)
    
    # If same protein_id appears in multiple files, prefer non-null values
    # Group by protein_id and fill missing values
    df_combined = df_combined.groupby('protein_id').agg({
        col: lambda x: x.dropna().iloc[0] if x.dropna().any() else None 
        for col in ORTHOLOG_COLUMNS if col in df_combined.columns
    }).reset_index()
    
    logger.info(f"Combined {len(df_combined)} unique protein_ids from all mapped files")
    
    # Create lookup dictionary for faster matching
    logger.info("")
    logger.info("Creating lookup dictionary...")
    ortholog_lookup = {}
    for _, row in df_combined.iterrows():
        protein_id = row['protein_id']
        if pd.notna(protein_id):
            ortholog_lookup[protein_id] = {
                col: row.get(col) for col in ORTHOLOG_COLUMNS if col in df_combined.columns
            }
    
    logger.info(f"Created lookup for {len(ortholog_lookup)} protein_ids")
    
    # Match and copy data
    logger.info("")
    logger.info("Matching and copying ortholog data...")
    matches_found = 0
    
    for idx, row in df_input.iterrows():
        protein_id = row['protein_id']
        
        if pd.notna(protein_id) and protein_id in ortholog_lookup:
            ortholog_data = ortholog_lookup[protein_id]
            
            # Copy each column if it has a value
            for col in ORTHOLOG_COLUMNS:
                if col in ortholog_data and pd.notna(ortholog_data[col]):
                    # Only update if current value is empty/None
                    if pd.isna(df_input.at[idx, col]) or df_input.at[idx, col] is None:
                        df_input.at[idx, col] = ortholog_data[col]
                        matches_found += 1
    
    logger.info(f"Copied ortholog data for {matches_found} entries")
    
    # Count how many rows got complete ortholog data
    complete_data = (
        df_input['human_ortholog_id'].notna() &
        df_input['human_ortholog_description'].notna() &
        df_input['human_ortholog_sequence'].notna()
    ).sum()
    
    logger.info(f"Rows with complete ortholog data: {complete_data} / {len(df_input)}")
    
    # Save output file
    logger.info("")
    logger.info(f"Saving results to {output_file}...")
    try:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        df_input.to_csv(output_file, index=False)
        logger.info(f"✓ Successfully saved {len(df_input)} rows")
    except Exception as e:
        logger.error(f"Error saving CSV file: {e}")
        return
    
    # Print summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    logger.info(f"Input rows: {len(df_input)}")
    logger.info(f"Unique protein_ids in input: {df_input['protein_id'].nunique()}")
    logger.info(f"Unique protein_ids in mapped files: {len(ortholog_lookup)}")
    logger.info(f"Rows with complete ortholog data: {complete_data}")
    logger.info(f"Output file: {output_file}")
    logger.info("")


def main():
    """Main workflow function."""
    parser = argparse.ArgumentParser(
        description='Copy human ortholog data from existing mapped files to a new full_disorder.csv'
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Path to input full_disorder.csv file'
    )
    parser.add_argument(
        '--mapped-file1',
        type=str,
        required=True,
        help='Path to first human_ortholog_mapped.csv file'
    )
    parser.add_argument(
        '--mapped-file2',
        type=str,
        default=None,
        help='Path to second human_ortholog_mapped.csv file (optional)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Path to output CSV file (default: input file with _with_orthologs suffix)'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        return
    
    # Determine output file path
    if args.output:
        output_path = args.output
    else:
        base, ext = os.path.splitext(args.input)
        output_path = f"{base}_with_orthologs{ext}"
    
    # Collect mapped files
    mapped_files = [args.mapped_file1]
    if args.mapped_file2:
        mapped_files.append(args.mapped_file2)
    
    # Run the copy operation
    copy_ortholog_data(args.input, mapped_files, output_path)


if __name__ == '__main__':
    main()


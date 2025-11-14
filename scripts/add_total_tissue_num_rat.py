#!/usr/bin/env python3
"""
Add total_tissue_num column to rat full_disorder.csv

This script reads the rat full_disorder.csv file and adds a total_tissue_num column
that sums up the tissue presence values across the specified tissue columns.
The column is inserted before the 'all_brain' column.

Usage:
    python add_total_tissue_num_rat.py [input_file] [output_file]
    
    If output_file is not specified, the input file will be overwritten.
"""

import pandas as pd
import sys
import os
from pathlib import Path

def add_total_tissue_num(input_file, output_file=None):
    """
    Add total_tissue_num column to rat full_disorder.csv
    
    Args:
        input_file (str): Path to input CSV file
        output_file (str, optional): Path to output CSV file. If None, overwrites input.
    """
    # Tissue columns to sum (excluding all_brain from the sum list, but we'll include it)
    tissue_columns = [
        'all_brain',
        'testicle',
        'pancreas',
        'stomach',
        'liver',
        'fat',
        'intestine',
        'kidney',
        'spleen',
        'thymus',
        'lung',
        'muscle',
        'heart',
        'blood'
    ]
    
    print(f"Reading {input_file}...")
    df = pd.read_csv(input_file)
    
    print(f"Loaded {len(df)} rows")
    print(f"Columns: {list(df.columns)}")
    
    # Check which tissue columns exist in the dataframe
    existing_tissue_cols = [col for col in tissue_columns if col in df.columns]
    missing_cols = [col for col in tissue_columns if col not in df.columns]
    
    if missing_cols:
        print(f"Warning: The following tissue columns are missing: {missing_cols}")
        print(f"Will only sum across existing columns: {existing_tissue_cols}")
    
    if not existing_tissue_cols:
        print("Error: None of the specified tissue columns were found in the dataframe!")
        return
    
    # Fill NaN values with 0 and convert to numeric, then sum across tissue columns
    print(f"Calculating total_tissue_num by summing: {existing_tissue_cols}")
    df['total_tissue_num'] = df[existing_tissue_cols].fillna(0).astype(float).sum(axis=1).astype(int)
    
    # Find the position of 'all_brain' column
    if 'all_brain' not in df.columns:
        print("Warning: 'all_brain' column not found. Inserting total_tissue_num at the end.")
        # If all_brain doesn't exist, just add at the end
        cols = df.columns.tolist()
        cols.remove('total_tissue_num')
        cols.append('total_tissue_num')
        df = df[cols]
    else:
        # Get column order
        cols = df.columns.tolist()
        # Remove total_tissue_num from its current position
        cols.remove('total_tissue_num')
        # Find index of all_brain
        all_brain_idx = cols.index('all_brain')
        # Insert total_tissue_num before all_brain
        cols.insert(all_brain_idx, 'total_tissue_num')
        # Reorder dataframe
        df = df[cols]
    
    # Determine output file
    if output_file is None:
        output_file = input_file
        print(f"Overwriting {input_file}...")
    else:
        print(f"Writing to {output_file}...")
    
    # Save the dataframe
    df.to_csv(output_file, index=False)
    
    print(f"Successfully added total_tissue_num column!")
    print(f"Total tissue number statistics:")
    print(df['total_tissue_num'].describe())
    print(f"\nDistribution:")
    print(df['total_tissue_num'].value_counts().sort_index())


def main():
    """Main function"""
    if len(sys.argv) < 2:
        # Default to rat full_disorder.csv
        script_dir = Path(__file__).parent
        project_root = script_dir.parent
        default_input = project_root / "data" / "processed" / "rat" / "full_disorder.csv"
        
        if not default_input.exists():
            print(f"Error: Default input file not found: {default_input}")
            print("\nUsage: python add_total_tissue_num_rat.py [input_file] [output_file]")
            print("  input_file: Path to rat full_disorder.csv")
            print("  output_file: (optional) Path to output file. If not specified, input will be overwritten.")
            sys.exit(1)
        
        input_file = str(default_input)
        output_file = None
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.exists(input_file):
        print(f"Error: Input file not found: {input_file}")
        sys.exit(1)
    
    add_total_tissue_num(input_file, output_file)


if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""
Compare PUNCH2 Results with Input Data

This script compares CSV files from a PUNCH2 results directory with the ID_matches 
column in an input file to identify which protein IDs did not get PUNCH2 predictions.

Usage:
    python compare_punch2_results.py input_file.csv punch2_directory [output_file]
    
Example:
    python compare_punch2_results.py data/processed/mouse/full_total_blast_aligned.csv data/processed/mouse/Punch2_results_mouse missing_predictions.txt
"""

import pandas as pd
import os
import sys
import glob
from pathlib import Path

def extract_protein_id_from_filename(filename):
    """
    Extract protein ID from PUNCH2 result filename.
    
    Args:
        filename (str): Filename like 'sp|A0A087WPF7|A0A087WPF7_MOUSE.csv'
        
    Returns:
        str: Protein ID like 'A0A087WPF7'
    """
    # Remove .csv extension
    name = filename.replace('.csv', '')
    
    # Split by | and get the middle part (ID_matches)
    parts = name.split('|')
    if len(parts) >= 2:
        return parts[1]  # Return the middle part (ID_matches)
    else:
        return None

def get_punch2_protein_ids(punch2_directory):
    """
    Get all protein IDs from PUNCH2 result files.
    
    Args:
        punch2_directory (str): Path to PUNCH2 results directory
        
    Returns:
        set: Set of protein IDs that have PUNCH2 results
    """
    punch2_ids = set()
    
    # Find all CSV files in the directory
    pattern = os.path.join(punch2_directory, "*.csv")
    csv_files = glob.glob(pattern)
    
    print(f"Found {len(csv_files)} PUNCH2 result files")
    
    for csv_file in csv_files:
        filename = os.path.basename(csv_file)
        protein_id = extract_protein_id_from_filename(filename)
        
        if protein_id:
            punch2_ids.add(protein_id)
        else:
            print(f"Warning: Could not extract protein ID from {filename}")
    
    return punch2_ids

def get_input_protein_ids(input_file):
    """
    Get all protein IDs from the input file's ID_matches column.
    
    Args:
        input_file (str): Path to input CSV file
        
    Returns:
        set: Set of protein IDs from input file
    """
    print(f"Reading input file: {input_file}")
    
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading input file: {e}")
        return set()
    
    # Check if ID_matches column exists
    if 'ID_matches' not in df.columns:
        print(f"Error: 'ID_matches' column not found in input file")
        print(f"Available columns: {list(df.columns)}")
        return set()
    
    # Get unique protein IDs from ID_matches column
    input_ids = set(df['ID_matches'].dropna().astype(str))
    
    print(f"Found {len(input_ids)} unique protein IDs in input file")
    
    return input_ids

def find_missing_predictions(input_file, punch2_directory, output_file=None):
    """
    Find protein IDs that don't have PUNCH2 predictions.
    
    Args:
        input_file (str): Path to input CSV file
        punch2_directory (str): Path to PUNCH2 results directory
        output_file (str, optional): Path to output file for missing IDs
    """
    
    print("=" * 60)
    print("PUNCH2 Results Comparison")
    print("=" * 60)
    
    # Get protein IDs from input file
    input_ids = get_input_protein_ids(input_file)
    if not input_ids:
        return
    
    # Get protein IDs from PUNCH2 results
    punch2_ids = get_punch2_protein_ids(punch2_directory)
    
    # Find missing predictions
    missing_ids = input_ids - punch2_ids
    successful_ids = input_ids & punch2_ids
    
    # Print summary
    print(f"\nSummary:")
    print(f"Total proteins in input: {len(input_ids)}")
    print(f"Proteins with PUNCH2 results: {len(successful_ids)}")
    print(f"Proteins missing PUNCH2 results: {len(missing_ids)}")
    print(f"Success rate: {len(successful_ids)/len(input_ids)*100:.1f}%")
    
    # Show some examples of missing proteins
    if missing_ids:
        print(f"\nFirst 10 missing protein IDs:")
        for i, protein_id in enumerate(sorted(missing_ids)):
            if i < 10:
                print(f"  {protein_id}")
            else:
                break
        
        if len(missing_ids) > 10:
            print(f"  ... and {len(missing_ids) - 10} more")
    
    # Save missing IDs to file if requested
    if output_file and missing_ids:
        with open(output_file, 'w') as f:
            f.write("Missing PUNCH2 Predictions\n")
            f.write("=" * 30 + "\n")
            f.write(f"Total missing: {len(missing_ids)}\n")
            f.write(f"Success rate: {len(successful_ids)/len(input_ids)*100:.1f}%\n\n")
            f.write("Missing Protein IDs:\n")
            for protein_id in sorted(missing_ids):
                f.write(f"{protein_id}\n")
        
        print(f"\nMissing protein IDs saved to: {output_file}")
    
    return missing_ids, successful_ids

def main():
    """Main function to handle command line arguments."""
    
    if len(sys.argv) < 3:
        print("Usage: python compare_punch2_results.py input_file.csv punch2_directory [output_file]")
        print("\nExample:")
        print("  python compare_punch2_results.py data/processed/mouse/full_total_blast_aligned.csv data/processed/mouse/Punch2_results_mouse")
        print("  python compare_punch2_results.py data/processed/mouse/full_total_blast_aligned.csv data/processed/mouse/Punch2_results_mouse missing_predictions.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    punch2_directory = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    # Check if PUNCH2 directory exists
    if not os.path.exists(punch2_directory):
        print(f"Error: PUNCH2 directory '{punch2_directory}' not found")
        sys.exit(1)
    
    # Run the comparison
    find_missing_predictions(input_file, punch2_directory, output_file)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Simple Motif Position Corrector

This script finds rows where position < 7 and updates the motif_position column
to match the position value.

Usage:
    python correct_aligned_positions.py input_csv [output_csv]
    
Example:
    python correct_aligned_positions.py data/processed/mouse/full_total_blast_aligned.csv data/processed/mouse/full_total_blast_aligned_corrected.csv
"""

import pandas as pd
import os
import sys

def process_csv(input_file, output_file):
    """
    Process CSV file to update motif_position for rows where position < 7.
    
    Args:
        input_file (str): Path to input CSV file
        output_file (str): Path to output CSV file
    """
    print("=" * 60)
    print("Motif Position Corrector")
    print("=" * 60)
    
    # Read input CSV
    print(f"Reading input CSV: {input_file}")
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    print(f"Found {len(df)} total rows")
    
    # Check required columns
    required_columns = ['position', 'motif_position']
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        print(f"Error: Missing required columns: {missing_columns}")
        print(f"Available columns: {list(df.columns)}")
        return
    
    # Find rows where position < 7
    position_less_than_7 = df['position'] < 7
    rows_to_update = df[position_less_than_7]
    
    print(f"Found {len(rows_to_update)} rows with position < 7")
    
    if len(rows_to_update) == 0:
        print("No rows to update. All positions are >= 7.")
        # Still write the output file
        df.to_csv(output_file, index=False)
        print(f"Output file written: {output_file}")
        return
    
    # Update motif_position for these rows
    updates_made = 0
    update_log = []
    
    for idx, row in rows_to_update.iterrows():
        old_motif_position = row['motif_position']
        new_motif_position = row['position']
        
        # Update the motif_position
        df.at[idx, 'motif_position'] = new_motif_position
        updates_made += 1
        
        # Log the change
        update_log.append(f"Row {idx+1}: motif_position {old_motif_position} → {new_motif_position}")
    
    # Print summary
    print(f"\nProcessing Summary:")
    print(f"  Total rows processed: {len(df)}")
    print(f"  Rows with position < 7: {len(rows_to_update)}")
    print(f"  Updates made: {updates_made}")
    print(f"  Success rate: {updates_made/len(rows_to_update)*100:.1f}%")
    
    # Show first 10 updates
    if update_log:
        print(f"\nFirst 10 updates:")
        for log_entry in update_log[:10]:
            print(f"  {log_entry}")
        if len(update_log) > 10:
            print(f"  ... and {len(update_log) - 10} more")
    
    # Write output file
    print(f"\nWriting corrected CSV to: {output_file}")
    df.to_csv(output_file, index=False)
    
    print(f"Successfully updated {updates_made} rows")
    print(f"Output file: {output_file}")

def main():
    """Main function to handle command line arguments."""
    
    if len(sys.argv) < 2:
        print("Usage: python correct_aligned_positions.py input_csv [output_csv]")
        print("\nExample:")
        print("  python correct_aligned_positions.py data/processed/mouse/full_total_blast_aligned.csv")
        print("  python correct_aligned_positions.py data/processed/mouse/full_total_blast_aligned.csv data/processed/mouse/full_total_blast_aligned_corrected.csv")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    # Set default output filename if not provided
    if not output_file:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_corrected.csv"
    
    # Process the CSV file
    process_csv(input_file, output_file)

if __name__ == "__main__":
    main()

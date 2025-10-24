#!/usr/bin/env python3
"""
Integrate Disorder Data with Phosphorylation Data

This script processes phosphorylation data CSV files, keeping only relevant columns
and integrating PUNCH2 disorder predictions based on amino acid positions.

Usage:
    python integrate_disorder_data.py input_csv punch2_directory [output_file]
    
Example:
    python integrate_disorder_data.py data/processed/mouse/full_total_blast_aligned.csv outputs/punch2/Punch2_results_mouse
    python integrate_disorder_data.py data/processed/rat/full_total_blast_aligned.csv outputs/punch2/Punch2_results_rat missing_rat_disorder.csv
"""

import pandas as pd
import os
import sys
import csv
from pathlib import Path

def detect_species(input_file):
    """
    Detect species from the input file path.
    
    Args:
        input_file (str): Path to input CSV file
        
    Returns:
        str: 'mouse' or 'rat'
    """
    if 'data/processed/mouse' in input_file:
        return 'mouse'
    elif 'data/processed/rat' in input_file:
        return 'rat'
    else:
        # Default to mouse if path doesn't match
        return 'mouse'

def get_columns_for_species(species):
    """
    Get ordered list of columns to keep based on species.
    
    Args:
        species (str): 'mouse' or 'rat'
        
    Returns:
        list: Ordered list of column names to keep
    """
    # Common columns for both species
    common_columns = [
        'ID_matches',  # Will be renamed to protein_id
        'ID_types',
        'protein_description',
        'full_sequence',
        'length',
        'amino_acid',
        'aligned_position',
        'position_difference',
        'site_motif',
        'cleaned_site_motif'
    ]
    
    if species == 'mouse':
        tissue_columns = [
            'total_tissue_num',
            'brain',
            'brownfat',
            'heart',
            'kidney',
            'liver',
            'lung',
            'pancreas',
            'spleen',
            'testis'
        ]
    elif species == 'rat':
        tissue_columns = [
            'all_brain',
            'cortex',
            'brainstem',
            'cerebellum',
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
    else:
        raise ValueError(f"Unknown species: {species}")
    
    return common_columns + tissue_columns

def build_punch2_filename(id_types, id_matches, species):
    """
    Build PUNCH2 filename from protein information.
    
    Args:
        id_types (str): ID type (e.g., 'sp', 'tr')
        id_matches (str): Protein ID
        species (str): 'mouse' or 'rat'
        
    Returns:
        str: PUNCH2 filename
    """
    species_suffix = 'MOUSE' if species == 'mouse' else 'RAT'
    return f"{id_types}|{id_matches}|{id_matches}_{species_suffix}.csv"

def get_disorder_data(punch2_dir, filename, position, expected_aa):
    """
    Get disorder data from PUNCH2 CSV file by searching for matching position.
    
    Args:
        punch2_dir (str): Directory containing PUNCH2 results
        filename (str): PUNCH2 filename
        position (int or float): Position to search for in first column
        expected_aa (str): Expected amino acid
        
    Returns:
        tuple: (disorder_score, disordered, error_message)
    """
    filepath = os.path.join(punch2_dir, filename)
    
    if not os.path.exists(filepath):
        return None, None, f"PUNCH2 file not found: {filename}"
    
    try:
        # Convert position to integer, handling float values
        if pd.isna(position):
            return None, None, f"Position is NaN"
        
        position_int = int(float(position))
        
        with open(filepath, 'r') as f:
            reader = csv.reader(f)
            rows = list(reader)
        
        # Search for the position in the first column
        matching_row = None
        
        # Try exact position match
        for row in rows:
            if len(row) >= 1:
                try:
                    row_position = int(float(row[0]))  # First column is position
                    if row_position == position_int:
                        matching_row = row
                        break
                except (ValueError, IndexError):
                    continue
        
        if matching_row is None:
            return None, None, f"Position {position_int} not found in PUNCH2 file"
        
        if len(matching_row) < 4:
            return None, None, f"Invalid row format: {matching_row}"
        
        # Extract data
        row_aa = matching_row[1]  # Column 2 (amino acid)
        disorder_score = float(matching_row[2])  # Column 3 (disorder score)
        disordered = int(matching_row[3])  # Column 4 (disorder prediction: 1=disordered, 0=ordered)
        
        # Verify amino acid match
        if row_aa != expected_aa:
            return None, None, f"Amino acid mismatch: expected {expected_aa}, found {row_aa}"
        
        return disorder_score, disordered, None
        
    except Exception as e:
        return None, None, f"Error reading PUNCH2 file: {str(e)}"

def process_data(input_csv, punch2_dir, output_file):
    """
    Main processing function.
    
    Args:
        input_csv (str): Path to input CSV file
        punch2_dir (str): Path to PUNCH2 results directory
        output_file (str): Path to output CSV file
    """
    print("=" * 60)
    print("Disorder Data Integration")
    print("=" * 60)
    
    # Detect species
    species = detect_species(input_csv)
    print(f"Detected species: {species}")
    
    # Read input CSV
    print(f"Reading input CSV: {input_csv}")
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    print(f"Found {len(df)} total rows")
    
    # Get columns to keep
    columns_to_keep = get_columns_for_species(species)
    print(f"Keeping {len(columns_to_keep)} columns for {species} data")
    
    # Check if all required columns exist
    missing_columns = [col for col in columns_to_keep if col not in df.columns]
    if missing_columns:
        print(f"Error: Missing required columns: {missing_columns}")
        print(f"Available columns: {list(df.columns)}")
        return
    
    # Select and rename columns
    df_selected = df[columns_to_keep].copy()
    df_selected = df_selected.rename(columns={'ID_matches': 'protein_id'})
    
    # Add disorder columns
    df_selected['disorder_score'] = None
    df_selected['disordered?'] = None
    
    print(f"Processing {len(df_selected)} phosphorylation sites...")
    
    # Process each row
    successful_count = 0
    skipped_count = 0
    skipped_reasons = {}
    amino_acid_mismatches = []  # Track amino acid mismatches
    
    for idx, row in df_selected.iterrows():
        protein_id = row['protein_id']
        id_types = row['ID_types']
        amino_acid = row['amino_acid']
        aligned_position = row['aligned_position']
        
        # Skip rows with invalid positions
        if pd.isna(aligned_position):
            skipped_count += 1
            skipped_reasons['Position is NaN'] = skipped_reasons.get('Position is NaN', 0) + 1
            if skipped_count <= 10:
                print(f"  Skipped row {idx+1}: Position is NaN")
            continue
        
        # Build PUNCH2 filename
        punch2_filename = build_punch2_filename(id_types, protein_id, species)
        
        # Get disorder data
        disorder_score, disordered, error_msg = get_disorder_data(
            punch2_dir, punch2_filename, aligned_position, amino_acid
        )
        
        if error_msg:
            skipped_count += 1
            reason = error_msg.split(':')[0]  # Get first part of error message
            skipped_reasons[reason] = skipped_reasons.get(reason, 0) + 1
            
            # Track amino acid mismatches
            if "Amino acid mismatch" in error_msg:
                amino_acid_mismatches.append({
                    'row': idx + 1,
                    'protein_id': protein_id,
                    'position': aligned_position,
                    'expected_aa': amino_acid,
                    'error': error_msg
                })
            
            if skipped_count <= 10:  # Show first 10 errors
                print(f"  Skipped row {idx+1}: {error_msg}")
                # Debug info for amino acid mismatches
                if "Amino acid mismatch" in error_msg:
                    print(f"    Debug - Position: {aligned_position}, Expected: {amino_acid}, Protein: {protein_id}")
        else:
            successful_count += 1
            df_selected.at[idx, 'disorder_score'] = disorder_score
            df_selected.at[idx, 'disordered?'] = disordered
    
    # Print summary
    print(f"\nProcessing Summary:")
    print(f"  Total rows processed: {len(df_selected)}")
    print(f"  Successfully integrated: {successful_count}")
    print(f"  Skipped: {skipped_count}")
    print(f"  Success rate: {successful_count/len(df_selected)*100:.1f}%")
    
    if skipped_reasons:
        print(f"\nSkipped reasons:")
        for reason, count in skipped_reasons.items():
            print(f"  {reason}: {count} rows")
    
    # List all amino acid mismatches
    if amino_acid_mismatches:
        print(f"\nAmino Acid Mismatches ({len(amino_acid_mismatches)} total):")
        print("=" * 80)
        for mismatch in amino_acid_mismatches:
            print(f"Row {mismatch['row']}: Protein {mismatch['protein_id']}, Position {mismatch['position']}, Expected {mismatch['expected_aa']}")
            print(f"  Error: {mismatch['error']}")
        print("=" * 80)
    
    # Write output file
    print(f"\nWriting output to: {output_file}")
    df_selected.to_csv(output_file, index=False)
    
    print(f"Successfully created disorder-integrated dataset with {len(df_selected)} rows")
    print(f"Output file: {output_file}")

def main():
    """Main function to handle command line arguments."""
    
    if len(sys.argv) < 3:
        print("Usage: python integrate_disorder_data.py input_csv punch2_directory [output_file]")
        print("\nExample:")
        print("  python integrate_disorder_data.py data/processed/mouse/full_total_blast_aligned.csv outputs/punch2/Punch2_results_mouse")
        print("  python integrate_disorder_data.py data/processed/rat/full_total_blast_aligned.csv outputs/punch2/Punch2_results_rat missing_rat_disorder.csv")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    punch2_dir = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    # Check if input file exists
    if not os.path.exists(input_csv):
        print(f"Error: Input file '{input_csv}' not found")
        sys.exit(1)
    
    # Check if PUNCH2 directory exists
    if not os.path.exists(punch2_dir):
        print(f"Error: PUNCH2 directory '{punch2_dir}' not found")
        sys.exit(1)
    
    # Set default output filename if not provided
    if not output_file:
        input_dir = os.path.dirname(input_csv)
        output_file = os.path.join(input_dir, "full_disorder.csv")
    
    # Process the data
    process_data(input_csv, punch2_dir, output_file)

if __name__ == "__main__":
    main()

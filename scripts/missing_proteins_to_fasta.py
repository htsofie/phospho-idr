#!/usr/bin/env python3
"""
Generate FASTA file from missing PUNCH2 proteins

This script takes a list of missing protein IDs and generates a FASTA file
by looking up their sequences in the original CSV file.

Usage:
    python missing_proteins_to_fasta.py missing_ids_file.csv input_csv_file [output_fasta_file]
    
Example:
    python missing_proteins_to_fasta.py outputs/punch2/missed_proteins data/processed/mouse/full_total_blast_aligned.csv missing_proteins.fasta
"""

import pandas as pd
import sys
import os
from pathlib import Path

def detect_species_from_path(csv_file):
    """
    Detect species from the CSV file path.
    
    Args:
        csv_file (str): Path to CSV file
        
    Returns:
        tuple: (species_id, species_name, tax_id)
    """
    if 'data/processed/mouse' in csv_file:
        return 'MOUSE', 'Mus musculus', '10090'
    elif 'data/processed/rat' in csv_file:
        return 'RAT', 'Rattus norvegicus', '10116'
    else:
        # Default to mouse if path doesn't match
        return 'MOUSE', 'Mus musculus', '10090'

def read_missing_ids(missing_ids_file):
    """
    Read missing protein IDs from file.
    
    Args:
        missing_ids_file (str): Path to file containing missing protein IDs
        
    Returns:
        set: Set of missing protein IDs
    """
    missing_ids = set()
    
    try:
        with open(missing_ids_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and not line.startswith('Missing'):
                    missing_ids.add(line)
    except Exception as e:
        print(f"Error reading missing IDs file: {e}")
        return set()
    
    print(f"Read {len(missing_ids)} missing protein IDs")
    return missing_ids

def create_fasta_from_missing_ids(missing_ids_file, input_csv_file, output_fasta_file=None):
    """
    Create FASTA file from missing protein IDs.
    
    Args:
        missing_ids_file (str): Path to file containing missing protein IDs
        input_csv_file (str): Path to original CSV file with sequences
        output_fasta_file (str, optional): Path to output FASTA file
    """
    
    print("=" * 60)
    print("Missing Proteins to FASTA Converter")
    print("=" * 60)
    
    # Read missing protein IDs
    missing_ids = read_missing_ids(missing_ids_file)
    if not missing_ids:
        print("No missing protein IDs found")
        return
    
    # Detect species from input CSV path
    species_id, species_name, tax_id = detect_species_from_path(input_csv_file)
    print(f"Detected species: {species_name} ({species_id})")
    
    # Read input CSV file
    print(f"Reading input CSV file: {input_csv_file}")
    try:
        df = pd.read_csv(input_csv_file)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    # Check required columns
    required_columns = ['ID_matches', 'ID_types', 'protein_description', 'full_sequence']
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        print(f"Error: Missing required columns: {missing_columns}")
        print(f"Available columns: {list(df.columns)}")
        return
    
    # Filter for missing protein IDs
    missing_df = df[df['ID_matches'].isin(missing_ids)]
    print(f"Found {len(missing_df)} entries for missing protein IDs")
    
    if len(missing_df) == 0:
        print("No matching entries found for missing protein IDs")
        return
    
    # Group by ID_matches to get unique proteins
    unique_missing = missing_df.groupby('ID_matches').first().reset_index()
    print(f"Found {len(unique_missing)} unique missing proteins")
    
    # Filter out proteins without sequences
    proteins_with_sequences = []
    skipped_count = 0
    
    for _, protein in unique_missing.iterrows():
        sequence = str(protein['full_sequence'])
        id_types = str(protein['ID_types'])
        id_matches = str(protein['ID_matches'])
        description = str(protein['protein_description'])
        
        # Skip if sequence is empty, NaN, or 'nan'
        if pd.isna(sequence) or sequence.strip() == '' or sequence.lower() == 'nan':
            skipped_count += 1
            continue
        
        # Skip if ID_types or ID_matches are None, NaN, or 'nan'
        if (pd.isna(id_types) or id_types.lower() in ['none', 'nan', ''] or 
            pd.isna(id_matches) or id_matches.lower() in ['none', 'nan', '']):
            skipped_count += 1
            continue
        
        proteins_with_sequences.append(protein)
    
    print(f"Skipped {skipped_count} proteins (missing sequences or ID information)")
    print(f"Processing {len(proteins_with_sequences)} proteins with complete data")
    
    if len(proteins_with_sequences) == 0:
        print("No proteins with complete data found")
        return
    
    # Set default output filename if not provided
    if not output_fasta_file:
        output_fasta_file = f"missing_proteins_{species_id.lower()}.fasta"
    
    print(f"Writing {len(proteins_with_sequences)} proteins to {output_fasta_file}")
    
    # Write FASTA file
    with open(output_fasta_file, 'w') as f:
        for protein in proteins_with_sequences:
            # Create FASTA header in the format: >tr|ID|ID_SPECIES description OS=species OX=taxid GN=gene PE=1 SV=1
            protein_id = protein['ID_matches']
            id_types = str(protein['ID_types'])
            id_matches = str(protein['ID_matches'])
            description = str(protein['protein_description'])
            
            # Create header with detected species
            header = f'>{id_types}|{id_matches}|{id_matches}_{species_id} {description} OS={species_name} OX={tax_id} GN={id_matches} PE=1 SV=1'
            
            # Get sequence (we already know it exists)
            sequence = str(protein['full_sequence']).strip()
            
            # Write header
            f.write(f"{header}\n")
            
            # Write sequence with line wrapping (60 characters per line)
            for i in range(0, len(sequence), 60):
                f.write(f"{sequence[i:i+60]}\n")
    
    print(f"Successfully generated FASTA file: {output_fasta_file}")
    print(f"Total proteins written: {len(proteins_with_sequences)}")
    
    return output_fasta_file

def main():
    """Main function to handle command line arguments."""
    
    if len(sys.argv) < 3:
        print("Usage: python missing_proteins_to_fasta.py missing_ids_file input_csv_file [output_fasta_file]")
        print("\nExample:")
        print("  python missing_proteins_to_fasta.py outputs/punch2/missed_proteins data/processed/mouse/full_total_blast_aligned.csv")
        print("  python missing_proteins_to_fasta.py outputs/punch2/missed_proteins data/processed/mouse/full_total_blast_aligned.csv missing_mouse_proteins.fasta")
        sys.exit(1)
    
    missing_ids_file = sys.argv[1]
    input_csv_file = sys.argv[2]
    output_fasta_file = sys.argv[3] if len(sys.argv) > 3 else None
    
    # Check if missing IDs file exists
    if not os.path.exists(missing_ids_file):
        print(f"Error: Missing IDs file '{missing_ids_file}' not found")
        sys.exit(1)
    
    # Check if input CSV file exists
    if not os.path.exists(input_csv_file):
        print(f"Error: Input CSV file '{input_csv_file}' not found")
        sys.exit(1)
    
    # Create FASTA file from missing proteins
    create_fasta_from_missing_ids(missing_ids_file, input_csv_file, output_fasta_file)

if __name__ == "__main__":
    main()

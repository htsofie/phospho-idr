#!/usr/bin/env python3
"""
CSV to FASTA Converter

This script takes a CSV file with protein data and generates FASTA files.
- Groups proteins by IPI ID (Protein column)
- Generates one FASTA entry per unique protein
- Skips proteins without sequences in "full_sequence" column
- Configurable number of proteins per file (default: 65, use 0 for unlimited)
- Auto-detects species from input path (mouse/rat)
- FASTA header format: >tr|ID|ID_SPECIES description OS=species OX=taxid GN=gene PE=1 SV=1
- Sequences are wrapped to 60 characters per line for readability
- Sequence from "full_sequence" column

Species Detection:
- Paths containing 'data/processed/mouse' → ID_MOUSE, Mus musculus, OX=10090
- Paths containing 'data/processed/rat' → ID_RAT, Rattus norvegicus, OX=10116

Usage:
    python csv_to_fasta.py input.csv [output_prefix] [max_proteins_per_file]
    
Examples:
    python csv_to_fasta.py data/processed/mouse/proteins.csv                    # Mouse proteins, 65 per file
    python csv_to_fasta.py data/processed/rat/proteins.csv rat_proteins 100     # Rat proteins, 100 per file
    python csv_to_fasta.py data/processed/mouse/proteins.csv mouse_proteins 0   # Mouse proteins, unlimited
"""

import pandas as pd
import sys
import os
from collections import defaultdict

def process_csv_to_fasta(csv_file, output_prefix="proteins", max_proteins_per_file=65):
    """
    Process CSV file and generate FASTA files.
    
    Args:
        csv_file (str): Path to input CSV file
        output_prefix (str): Prefix for output FASTA files
        max_proteins_per_file (int): Maximum proteins per FASTA file (default: 65, use 0 for unlimited)
    """
    
    # Determine species based on input path
    if 'data/processed/mouse' in csv_file:
        species_id = 'MOUSE'
        species_name = 'Mus musculus'
        tax_id = '10090'
    elif 'data/processed/rat' in csv_file:
        species_id = 'RAT'
        species_name = 'Rattus norvegicus'
        tax_id = '10116'
    else:
        # Default to rat if path doesn't match
        species_id = 'RAT'
        species_name = 'Rattus norvegicus'
        tax_id = '10116'
    
    print(f"Detected species: {species_name} ({species_id})")
    
    print(f"Reading CSV file: {csv_file}")
    
    # Read CSV file
    try:
        df = pd.read_csv(csv_file)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return
    
    # Check required columns
    required_columns = ['Protein', 'ID_types', 'ID_matches', 'protein_description', 'full_sequence']
    missing_columns = [col for col in required_columns if col not in df.columns]
    
    if missing_columns:
        print(f"Error: Missing required columns: {missing_columns}")
        print(f"Available columns: {list(df.columns)}")
        return
    
    print(f"Found {len(df)} total entries")
    
    # Group by Protein (IPI ID) and get unique proteins
    unique_proteins = df.groupby('Protein').first().reset_index()
    print(f"Found {len(unique_proteins)} unique proteins")
    
    # Filter out proteins without sequences or with missing ID information BEFORE processing
    proteins_with_sequences = []
    skipped_count = 0
    
    for _, protein in unique_proteins.iterrows():
        sequence = str(protein['full_sequence'])
        id_types = str(protein['ID_types'])
        id_matches = str(protein['ID_matches'])
        
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
    
    # Handle unlimited proteins (single file)
    if max_proteins_per_file == 0:
        num_files = 1
        print(f"Will generate 1 FASTA file with all {len(proteins_with_sequences)} proteins")
    else:
        # Calculate number of FASTA files needed
        num_files = (len(proteins_with_sequences) + max_proteins_per_file - 1) // max_proteins_per_file
        print(f"Will generate {num_files} FASTA file(s) with max {max_proteins_per_file} proteins each")
    
    # Generate FASTA files
    for file_num in range(num_files):
        if max_proteins_per_file == 0:
            # Unlimited - use all proteins
            start_idx = 0
            end_idx = len(proteins_with_sequences)
        else:
            # Limited - calculate subset
            start_idx = file_num * max_proteins_per_file
            end_idx = min((file_num + 1) * max_proteins_per_file, len(proteins_with_sequences))
        
        # Get subset of proteins for this file
        proteins_subset = proteins_with_sequences[start_idx:end_idx]
        
        # Generate output filename
        if num_files == 1:
            output_file = f"{output_prefix}.fasta"
        else:
            output_file = f"{output_prefix}_part{file_num + 1:03d}.fasta"
        
        print(f"Writing {len(proteins_subset)} proteins to {output_file}")
        
        # Write FASTA file
        with open(output_file, 'w') as f:
            for protein in proteins_subset:
                # Create FASTA header in the format: >tr|ID|ID_SPECIES description OS=species OX=taxid GN=gene PE=1 SV=1
                protein_id = protein['Protein']
                id_types = str(protein['ID_types'])
                id_matches = str(protein['ID_matches'])
                description = str(protein['protein_description'])
                
                # Create a more readable header format with detected species
                header = f'>{id_types}|{id_matches}|{id_matches}_{species_id} {description} OS={species_name} OX={tax_id} GN={id_matches} PE=1 SV=1'
                
                # Get sequence (we already know it exists)
                sequence = str(protein['full_sequence']).strip()
                
                # Write header
                f.write(f"{header}\n")
                
                # Write sequence with line wrapping (60 characters per line)
                for i in range(0, len(sequence), 60):
                    f.write(f"{sequence[i:i+60]}\n")
        
        print(f"Completed {output_file}")
    
    print(f"Successfully generated {num_files} FASTA file(s) with {len(proteins_with_sequences)} total proteins")

def main():
    """Main function to handle command line arguments."""
    
    if len(sys.argv) < 2:
        print("Usage: python csv_to_fasta.py input.csv [output_prefix] [max_proteins_per_file]")
        print("Example: python csv_to_fasta.py proteins.csv mouse_proteins 65")
        print("Example: python csv_to_fasta.py proteins.csv mouse_proteins 0  # unlimited")
        print("Default: max_proteins_per_file = 65")
        sys.exit(1)
    
    csv_file = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else "proteins"
    max_proteins = int(sys.argv[3]) if len(sys.argv) > 3 else 65
    
    # Check if input file exists
    if not os.path.exists(csv_file):
        print(f"Error: Input file '{csv_file}' not found")
        sys.exit(1)
    
    # Process the CSV file
    process_csv_to_fasta(csv_file, output_prefix, max_proteins)

if __name__ == "__main__":
    main()

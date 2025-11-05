#!/usr/bin/env python3
"""
Interactome BLAST Homology Search

This script BLASTs sequences from full_disorder.csv against an interactome BLAST database
and adds homologous protein information to the output CSV file.

Usage:
    python interactome_blast_search.py <full_disorder_file> <blast_database_path> <output_csv>

Arguments:
    full_disorder_file    Path to full_disorder.csv file
    blast_database_path   Path to BLAST database (without extension, e.g., output/interactome_proteins_blast)
    output_csv           Path to output CSV file with BLAST results

Example:
    python interactome_blast_search.py data/processed/mouse/full_disorder.csv \
        output/interactome_proteins_blast \
        output/full_disorder_with_interactome.csv
"""

from Bio import SeqIO
import pandas as pd
import sys
import os
import subprocess
import tempfile

def blast_sequence_against_database(sequence, db_path, evalue=1e-5, max_target_seqs=50):
    """
    BLAST a single sequence against the interactome database.
    
    Args:
        sequence (str): Protein sequence to BLAST
        db_path (str): Path to BLAST database (without extension)
        evalue (float): E-value threshold (default: 1e-5)
        max_target_seqs (int): Maximum number of target sequences (default: 50)
    
    Returns:
        list: List of dictionaries with hit information
        Format: [{'subject_id': '...', 'identity': 95.5, 'coverage': 80.0, 'evalue': 1e-10}, ...]
    """
    if not sequence or pd.isna(sequence):
        return []
    
    try:
        # Create temporary input FASTA file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_input:
            temp_input.write(f">query\n{sequence}\n")
            temp_input_path = temp_input.name
        
        # Create temporary output file for tabular results
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_output:
            temp_output_path = temp_output.name
        
        # Run BLAST with tabular output format
        # Format: qseqid sseqid pident length qstart qend evalue bitscore
        cmd = [
            'blastp', '-query', temp_input_path, '-db', db_path,
            '-out', temp_output_path, '-outfmt', '6 qseqid sseqid pident length qstart qend evalue bitscore',
            '-evalue', str(evalue), '-max_target_seqs', str(max_target_seqs),
            '-word_size', '3', '-matrix', 'BLOSUM62'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        # Clean up temporary input file
        os.unlink(temp_input_path)
        
        if result.returncode != 0:
            # Clean up output file if it exists
            if os.path.exists(temp_output_path):
                os.unlink(temp_output_path)
            return []
        
        # Parse tabular results
        # Use dictionary to track best hit per subject (in case of multiple HSPs)
        hits_dict = {}
        hits = []  # Initialize hits list
        
        if os.path.exists(temp_output_path) and os.path.getsize(temp_output_path) > 0:
            with open(temp_output_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) >= 8:
                        subject_id = fields[1]
                        identity = float(fields[2])
                        alignment_length = int(fields[3])
                        qstart = int(fields[4])
                        qend = int(fields[5])
                        evalue_hit = float(fields[6])
                        
                        # Calculate coverage as percentage of query sequence covered
                        query_length = len(sequence)
                        if query_length > 0:
                            # Handle both forward and reverse alignments
                            coverage_length = abs(qend - qstart) + 1
                            coverage = (coverage_length / query_length) * 100
                        else:
                            coverage = 0.0
                        
                        # Keep best hit per subject (lowest e-value, or highest identity if e-value is same)
                        if subject_id not in hits_dict:
                            hits_dict[subject_id] = {
                                'subject_id': subject_id,
                                'identity': identity,
                                'coverage': coverage,
                                'evalue': evalue_hit
                            }
                        else:
                            # Update if this hit is better (lower e-value or higher identity)
                            existing = hits_dict[subject_id]
                            if evalue_hit < existing['evalue'] or \
                               (evalue_hit == existing['evalue'] and identity > existing['identity']):
                                hits_dict[subject_id] = {
                                    'subject_id': subject_id,
                                    'identity': identity,
                                    'coverage': coverage,
                                    'evalue': evalue_hit
                                }
            
            # Convert dictionary to list
            hits = list(hits_dict.values())
            
            # Clean up temporary output file
            os.unlink(temp_output_path)
        else:
            # Clean up temporary output file if it exists but is empty
            if os.path.exists(temp_output_path):
                os.unlink(temp_output_path)
        
        return hits
        
    except Exception as e:
        print(f"Error running BLAST: {e}")
        return []

def blast_full_disorder_against_database(full_disorder_file, db_path, output_csv):
    """
    BLAST sequences from full_disorder.csv against the interactome database.
    
    Args:
        full_disorder_file (str): Path to full_disorder.csv file
        db_path (str): Path to BLAST database (without extension)
        output_csv (str): Path to output CSV file with BLAST results
    
    Returns:
        bool: True if successful, False otherwise
    """
    print("=" * 80)
    print("BLASTING FULL_DISORDER SEQUENCES AGAINST INTERACTOME DATABASE")
    print("=" * 80)
    
    if not os.path.exists(full_disorder_file):
        print(f"Error: full_disorder file not found: {full_disorder_file}")
        return False
    
    if not os.path.exists(f"{db_path}.phr"):
        print(f"Error: BLAST database not found: {db_path}")
        print(f"Expected database file: {db_path}.phr")
        return False
    
    print(f"Loading full_disorder file: {full_disorder_file}")
    try:
        df = pd.read_csv(full_disorder_file)
        print(f"Loaded {len(df)} rows from full_disorder file")
    except Exception as e:
        print(f"Error loading full_disorder file: {e}")
        return False
    
    if 'full_sequence' not in df.columns:
        print("Error: 'full_sequence' column not found in full_disorder file")
        print(f"Available columns: {list(df.columns)}")
        return False
    
    # Initialize new columns
    df['homologous_ids'] = ''
    df['sequence_identity_percentages'] = ''
    df['alignment_coverage'] = ''
    df['e_values'] = ''
    
    # Group by full_sequence to avoid redundant BLAST searches
    print(f"\nGrouping by unique sequences to optimize BLAST searches...")
    unique_sequences = df['full_sequence'].dropna().unique()
    print(f"Found {len(unique_sequences)} unique sequences to BLAST")
    
    # Create a mapping from sequence to BLAST results
    sequence_to_hits = {}
    
    for idx, sequence in enumerate(unique_sequences):
        if pd.isna(sequence) or not isinstance(sequence, str) or len(sequence.strip()) == 0:
            continue
        
        print(f"  BLASTing sequence {idx + 1}/{len(unique_sequences)} (length: {len(sequence)})...")
        hits = blast_sequence_against_database(sequence, db_path)
        
        if hits:
            # Store hits for this sequence
            sequence_to_hits[sequence] = hits
            print(f"    Found {len(hits)} hits")
        else:
            sequence_to_hits[sequence] = []
    
    # Apply results to all rows with matching sequences
    print(f"\nApplying BLAST results to dataframe...")
    print(f"Total rows in dataframe: {len(df)}")
    
    for idx, row in df.iterrows():
        sequence = row['full_sequence']
        
        if pd.isna(sequence) or sequence not in sequence_to_hits:
            continue
        
        hits = sequence_to_hits[sequence]
        
        if hits:
            # Create semicolon-separated lists
            homologous_ids = ';'.join([hit['subject_id'] for hit in hits])
            identities = ';'.join([f"{hit['identity']:.2f}" for hit in hits])
            coverages = ';'.join([f"{hit['coverage']:.2f}" for hit in hits])
            evalues = ';'.join([f"{hit['evalue']:.2e}" for hit in hits])
            
            df.at[idx, 'homologous_ids'] = homologous_ids
            df.at[idx, 'sequence_identity_percentages'] = identities
            df.at[idx, 'alignment_coverage'] = coverages
            df.at[idx, 'e_values'] = evalues
    
    # Count sequences with multiple rows
    sequence_counts = df['full_sequence'].value_counts()
    multi_row_sequences = sequence_counts[sequence_counts > 1]
    if len(multi_row_sequences) > 0:
        print(f"Found {len(multi_row_sequences)} sequences appearing in multiple rows")
        print(f"  Example: {sequence_counts.max()} rows share the same sequence")
    
    # Count rows with hits
    rows_with_hits = df[df['homologous_ids'] != ''].shape[0]
    print(f"Rows with BLAST hits: {rows_with_hits} / {len(df)}")
    
    # Save results
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    try:
        df.to_csv(output_csv, index=False)
        print(f"Successfully saved results to: {output_csv}")
        return True
    except Exception as e:
        print(f"Error saving output CSV: {e}")
        return False

def main():
    """Main function to handle command line arguments."""
    
    if len(sys.argv) != 4:
        print("Usage: python interactome_blast_search.py <full_disorder_file> <blast_database_path> <output_csv>")
        print("\nArguments:")
        print("  full_disorder_file    Path to full_disorder.csv file")
        print("  blast_database_path   Path to BLAST database (without extension)")
        print("  output_csv           Path to output CSV file with BLAST results")
        print("\nExample:")
        print("  python interactome_blast_search.py data/processed/mouse/full_disorder.csv \\")
        print("      output/interactome_proteins_blast \\")
        print("      output/full_disorder_with_interactome.csv")
        sys.exit(1)
    
    full_disorder_file = sys.argv[1]
    db_path = sys.argv[2]
    output_csv = sys.argv[3]
    
    # Run BLAST search
    success = blast_full_disorder_against_database(full_disorder_file, db_path, output_csv)
    
    if success:
        print("\n" + "=" * 80)
        print("✅ SUCCESS! BLAST search completed successfully")
        print("=" * 80)
        print(f"Output CSV: {output_csv}")
        print(f"Database path: {db_path}")
        print(f"\nAdded columns to output CSV:")
        print(f"  - homologous_ids (semicolon-separated)")
        print(f"  - sequence_identity_percentages (semicolon-separated)")
        print(f"  - alignment_coverage (semicolon-separated)")
        print(f"  - e_values (semicolon-separated)")
    else:
        print("\n❌ Failed to complete BLAST search")
        sys.exit(1)

if __name__ == "__main__":
    main()


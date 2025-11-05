#!/usr/bin/env python3
"""
Interactome BLAST Database Creator

This script always creates a BLAST database from a FASTA file.
Optionally, it can also create a FASTA file from an interactome CSV file
by extracting unique protein IDs and searching for them in a reference FASTA file.

Usage:
    # Create BLAST database from existing FASTA file:
    python interactome_blast_dbs.py <fasta_file> --db-path <db_path>
    
    # Create FASTA file from CSV and then create BLAST database:
    python interactome_blast_dbs.py <fasta_file> --db-path <db_path> \\
        --create-fasta <interactome_file> <search_fasta> <output_fasta>

Options:
    --create-fasta      Create FASTA file from interactome CSV file (optional)
    --db-path PATH     Path for BLAST database (required, without extension)

Example:
    # Create BLAST database from existing FASTA:
    python interactome_blast_dbs.py output/interactome_proteins.fasta \\
        --db-path output/interactome_blast
    
    # Create FASTA from CSV and then create BLAST database:
    python interactome_blast_dbs.py output/interactome_proteins.fasta \\
        --db-path output/interactome_blast \\
        --create-fasta data/interactome_insider/M_musculus_interfacesALL.csv \\
        data/blast_dbs/uniprotkb_taxonomy_id_10090_2025_10_06.fasta \\
        output/interactome_proteins.fasta
"""

from Bio import SeqIO
import pandas as pd
import sys
import os
import re
import subprocess
import argparse

def extract_unique_protein_ids(interactome_file):
    """
    Extract unique protein IDs from P1 and P2 columns of interactome file.
    
    Args:
        interactome_file (str): Path to interactome CSV file
        
    Returns:
        set: Set of unique protein IDs
    """
    print("=" * 80)
    print("EXTRACTING UNIQUE PROTEIN IDs FROM INTERACTOME FILE")
    print("=" * 80)
    
    if not os.path.exists(interactome_file):
        print(f"Error: Interactome file not found: {interactome_file}")
        return None
    
    print(f"Reading interactome file: {interactome_file}")
    try:
        # Read the CSV file - it may be tab-separated or comma-separated
        # Try comma first, then tab
        try:
            pairs = pd.read_csv(interactome_file, sep=',')
        except:
            pairs = pd.read_csv(interactome_file, sep='\t')
        
        print(f"Loaded {len(pairs)} interaction pairs")
        
        # Check if P1 and P2 columns exist
        if 'P1' not in pairs.columns or 'P2' not in pairs.columns:
            print(f"Error: 'P1' or 'P2' columns not found in interactome file")
            print(f"Available columns: {list(pairs.columns)}")
            return None
        
        # Extract unique IDs from P1 and P2 columns (matching user's code style)
        unique_ids = set(pairs['P1']).union(set(pairs['P2']))
        
        print(f"Found {len(unique_ids)} unique protein IDs")
        print(f"Sample IDs: {list(unique_ids)[:10]}")
        
        return unique_ids
        
    except Exception as e:
        print(f"Error reading interactome file: {e}")
        return None

def search_fasta_for_ids(search_fasta, unique_ids):
    """
    Search FASTA file for protein IDs and extract matching sequences.
    Optimized to use set lookups instead of nested loops.
    
    Args:
        search_fasta (str): Path to FASTA file to search
        unique_ids (set): Set of protein IDs to search for
        
    Returns:
        list: List of SeqRecord objects for matched sequences
    """
    print("\n" + "=" * 80)
    print("SEARCHING FASTA FILE FOR MATCHING IDs")
    print("=" * 80)
    
    if not os.path.exists(search_fasta):
        print(f"Error: Search FASTA file not found: {search_fasta}")
        return []
    
    print(f"Searching in: {search_fasta}")
    print(f"Looking for {len(unique_ids)} unique protein IDs")
    
    # Pre-process unique_ids for faster lookup
    # Create sets for different ID formats and mapping dictionaries
    unique_ids_set = set(unique_ids)
    unique_ids_clean = set()  # IDs without prefixes
    clean_to_original = {}  # Map clean ID to original ID(s)
    
    for query_id in unique_ids:
        clean_query = query_id.split('|')[-1] if '|' in query_id else query_id
        unique_ids_clean.add(clean_query)
        # Map clean ID to original (if multiple, take first)
        if clean_query not in clean_to_original:
            clean_to_original[clean_query] = query_id
    
    matched_records = []
    found_ids = set()
    processed_count = 0
    
    try:
        # Parse FASTA file
        for record in SeqIO.parse(search_fasta, "fasta"):
            processed_count += 1
            if processed_count % 100000 == 0:
                print(f"  Processed {processed_count:,} FASTA records, found {len(found_ids)} matches so far...")
            
            # Extract ID from FASTA header
            header_id = record.id
            description = record.description
            
            # Check for matches using set lookups (O(1) instead of O(n))
            matched = False
            matched_id = None
            
            # Direct match with full ID
            if header_id in unique_ids_set:
                matched_id = header_id
                matched = True
            # Check if clean header ID matches
            elif '|' in header_id:
                clean_header = header_id.split('|')[-1]
                if clean_header in unique_ids_clean:
                    matched_id = clean_to_original.get(clean_header)
                    matched = True
            else:
                # Check if header_id matches any clean query ID
                if header_id in unique_ids_clean:
                    matched_id = clean_to_original.get(header_id)
                    matched = True
            
            # Also check description for matches (but only if not already matched)
            # This is slower but necessary for IDs that appear in description
            if not matched and description:
                # Check if any query ID appears in description
                for query_id in unique_ids:
                    if query_id in description:
                        matched_id = query_id
                        matched = True
                        break
                    # Check clean version
                    clean_query = query_id.split('|')[-1] if '|' in query_id else query_id
                    if clean_query in description:
                        matched_id = query_id
                        matched = True
                        break
            
            if matched and matched_id not in found_ids:
                matched_records.append(record)
                found_ids.add(matched_id)
        
        print(f"Processed {processed_count:,} FASTA records")
        print(f"Found {len(matched_records)} matching sequences")
        print(f"Matched IDs: {len(found_ids)}")
        
        if len(found_ids) < len(unique_ids):
            missing_count = len(unique_ids) - len(found_ids)
            print(f"Warning: {missing_count} IDs from interactome file were not found in FASTA")
        
        return matched_records
        
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return []

def write_matched_fasta(matched_records, output_fasta):
    """
    Write matched FASTA records to output file.
    
    Args:
        matched_records (list): List of SeqRecord objects
        output_fasta (str): Path to output FASTA file
    
    Returns:
        bool: True if successful, False otherwise
    """
    print("\n" + "=" * 80)
    print("WRITING MATCHED SEQUENCES TO OUTPUT FASTA")
    print("=" * 80)
    
    if not matched_records:
        print("No matched records to write")
        return False
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_fasta)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")
    
    try:
        SeqIO.write(matched_records, output_fasta, "fasta")
        print(f"Successfully wrote {len(matched_records)} sequences to: {output_fasta}")
        return True
    except Exception as e:
        print(f"Error writing output FASTA file: {e}")
        return False

def create_blast_database(fasta_path, db_path):
    """
    Create BLAST database from FASTA file using makeblastdb.
    
    Args:
        fasta_path (str): Path to input FASTA file
        db_path (str): Path for output BLAST database (without extension)
    
    Returns:
        bool: True if successful, False otherwise
    """
    print("\n" + "=" * 80)
    print("CREATING BLAST DATABASE")
    print("=" * 80)
    
    if not os.path.exists(fasta_path):
        print(f"Error: FASTA file not found: {fasta_path}")
        return False
    
    print(f"Creating BLAST database from: {fasta_path}")
    print(f"Database path: {db_path}")
    
    try:
        # Run makeblastdb command
        cmd = ['makeblastdb', '-in', fasta_path, '-dbtype', 'prot', '-out', db_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error creating BLAST database:")
            print(f"STDOUT: {result.stdout}")
            print(f"STDERR: {result.stderr}")
            return False
        
        print(f"BLAST database created successfully!")
        print(f"Database files: {db_path}.phr, {db_path}.pin, {db_path}.psq")
        return True
        
    except FileNotFoundError:
        print("Error: makeblastdb command not found. Please ensure BLAST+ is installed and in your PATH.")
        return False
    except Exception as e:
        print(f"Error running makeblastdb: {e}")
        return False

def main():
    """Main function to handle command line arguments."""
    
    parser = argparse.ArgumentParser(
        description='Create BLAST database from FASTA file. Optionally create FASTA file from interactome CSV.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Create BLAST database from existing FASTA file:
  python interactome_blast_dbs.py output/interactome_proteins.fasta \\
      --db-path output/interactome_blast
  
  # Create FASTA file from CSV and then create BLAST database:
  python interactome_blast_dbs.py output/interactome_proteins.fasta \\
      --db-path output/interactome_blast \\
      --create-fasta data/interactome_insider/M_musculus_interfacesALL.csv \\
      data/blast_dbs/uniprotkb_taxonomy_id_10090_2025_10_06.fasta \\
      output/interactome_proteins.fasta
        """
    )
    
    parser.add_argument('fasta_file', help='Path to FASTA file (will be created if --create-fasta is used)')
    parser.add_argument('--db-path', type=str, required=True,
                       help='Path for BLAST database (without extension)')
    parser.add_argument('--create-fasta', nargs=3, metavar=('INTERACTOME', 'SEARCH_FASTA', 'OUTPUT_FASTA'),
                       help='Create FASTA file from interactome CSV. Requires: interactome_file search_fasta output_fasta')
    
    args = parser.parse_args()
    
    fasta_file = args.fasta_file
    db_path = args.db_path
    create_fasta = args.create_fasta
    
    # Optional: Create FASTA file from CSV
    if create_fasta:
        interactome_file, search_fasta, output_fasta = create_fasta
        
        # Step 1: Extract unique protein IDs from interactome file
        unique_ids = extract_unique_protein_ids(interactome_file)
        if unique_ids is None:
            print("\n❌ Failed to extract protein IDs from interactome file")
            sys.exit(1)
        
        # Step 2: Search FASTA file for matching IDs
        matched_records = search_fasta_for_ids(search_fasta, unique_ids)
        if not matched_records:
            print("\n❌ No matching sequences found in FASTA file")
            sys.exit(1)
        
        # Step 3: Write matched sequences to output FASTA file
        success = write_matched_fasta(matched_records, output_fasta)
        
        if not success:
            print("\n❌ Failed to write output FASTA file")
            sys.exit(1)
        
        print("\n" + "=" * 80)
        print("✅ SUCCESS! Output FASTA file created successfully")
        print("=" * 80)
        print(f"Output file: {output_fasta}")
        print(f"Total sequences: {len(matched_records)}")
        print(f"Matched IDs: {len(matched_records)} / {len(unique_ids)}")
        
        # Use the created FASTA file for BLAST database
        fasta_file = output_fasta
    
    # Always create BLAST database from FASTA file
    db_success = create_blast_database(fasta_file, db_path)
    
    if db_success:
        print("\n" + "=" * 80)
        print("✅ SUCCESS! BLAST database created successfully")
        print("=" * 80)
        print(f"Database path: {db_path}")
        print(f"FASTA file: {fasta_file}")
        print(f"You can now use this database for BLAST searches")
    else:
        print("\n❌ Failed to create BLAST database")
        sys.exit(1)

if __name__ == "__main__":
    main()

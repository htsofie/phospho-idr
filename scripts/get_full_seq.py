#!/usr/bin/env python3
"""
Extract full protein sequences from UniProt using uniprot_mapped_id.
"""

import pandas as pd
import argparse
import os
from typing import Optional, Tuple
from Bio import ExPASy
from Bio import SwissProt
import time


def get_uniprot_sequence(uniprot_id: str) -> Optional[Tuple[str, str]]:
    """Get full protein sequence from UniProt, handling isoforms."""
    if pd.isna(uniprot_id) or not uniprot_id:
        return None
    
    try:
        clean_id = str(uniprot_id).strip()
        
        # Handle multiple IDs separated by semicolons
        if ';' in clean_id:
            ids = [id.strip() for id in clean_id.split(';') if id.strip()]
            for single_id in ids:
                result = _get_single_uniprot_sequence(single_id)
                if result:
                    return result
            return None
        else:
            return _get_single_uniprot_sequence(clean_id)
    except Exception as e:
        print(f"  ✗ UniProt sequence fetch failed for {uniprot_id}: {e}")
        return None


def _get_single_uniprot_sequence(uniprot_id: str) -> Optional[Tuple[str, str]]:
    """Get sequence for a single UniProt ID, using specific isoform if available."""
    try:
        # First, try to get the sequence directly
        handle = ExPASy.get_sprot_raw(uniprot_id)
        record = SwissProt.read(handle)
        handle.close()
        
        sequence = record.sequence
        
        # Check if this is an isoform by looking at the original uniprot_id
        if '-' in uniprot_id and uniprot_id.split('-')[-1].isdigit():
            # This is a specific isoform - use it directly
            print(f"    Using specific isoform: {uniprot_id}")
            return sequence, "isoform"
        else:
            # This is canonical - use it directly
            print(f"    Using canonical sequence: {uniprot_id}")
            return sequence, "canonical"
            
    except Exception as e:
        print(f"  ✗ Single UniProt sequence fetch failed for {uniprot_id}: {e}")
        return None


def process_dataset(input_path: str, species: str, output_path: str) -> None:
    """Process the dataset to extract full protein sequences."""
    print(f"Loading dataset: {input_path}")
    df = pd.read_csv(input_path)
    print(f"Loaded {len(df)} rows")
    
    # Check if uniprot_mapped_id column exists
    if 'uniprot_mapped_id' not in df.columns:
        print("Error: 'uniprot_mapped_id' column not found in dataset")
        return
    
    # Add new columns for sequence data
    df['full_sequence'] = None
    df['sequence_length'] = None
    df['sequence_type'] = None
    df['sequence_fetch_success'] = False
    df['sequence_fetch_error'] = None
    
    # Process each row
    for idx, row in df.iterrows():
        print(f"\nRow {idx + 1}/{len(df)}")
        print(f"  UniProt ID: {row.get('uniprot_mapped_id', 'N/A')}")
        
        # Get full sequence from UniProt
        sequence_data = get_uniprot_sequence(row['uniprot_mapped_id'])
        
        if sequence_data:
            full_sequence, sequence_type = sequence_data
            df.at[idx, 'full_sequence'] = full_sequence
            df.at[idx, 'sequence_length'] = len(full_sequence)
            df.at[idx, 'sequence_type'] = sequence_type
            df.at[idx, 'sequence_fetch_success'] = True
            print(f"  ✓ Success: {sequence_type}, {len(full_sequence)} amino acids")
        else:
            df.at[idx, 'sequence_fetch_success'] = False
            df.at[idx, 'sequence_fetch_error'] = 'Could not fetch sequence from UniProt'
            print(f"  ✗ Failed to fetch sequence")
        
        # Add a small delay to be respectful to the APIs
        time.sleep(0.1)
    
    # Reorder columns to put sequence data after mapping results
    mapping_cols = ['uniprot_mapped_id', 'mapping_source', 'mapping_success']
    sequence_cols = ['full_sequence', 'sequence_length', 'sequence_type', 'sequence_fetch_success', 'sequence_fetch_error']
    other_cols = [col for col in df.columns if col not in mapping_cols + sequence_cols]
    df_reordered = df[mapping_cols + sequence_cols + other_cols]
    
    # Save the results
    print(f"\nSaving results to: {output_path}")
    df_reordered.to_csv(output_path, index=False)
    
    # Print summary statistics
    successful_fetches = df['sequence_fetch_success'].sum()
    print(f"\nProcessing Summary:")
    print(f"  Total rows: {len(df)}")
    print(f"  Successful sequence fetches: {successful_fetches}")
    print(f"  Failed sequence fetches: {len(df) - successful_fetches}")
    print(f"  Success rate: {successful_fetches/len(df)*100:.1f}%")
    
    # Print error breakdown
    if len(df) - successful_fetches > 0:
        print(f"\nError breakdown:")
        error_counts = df[~df['sequence_fetch_success']]['sequence_fetch_error'].value_counts()
        for error, count in error_counts.items():
            print(f"  {error}: {count}")


def main():
    parser = argparse.ArgumentParser(description='Extract full protein sequences from UniProt')
    parser.add_argument('--input', '-i', required=True, help='Path to input CSV file (mapped data)')
    parser.add_argument('--species', '-s', required=True, choices=['mouse', 'rat'], help='Species name')
    parser.add_argument('--output', '-o', help='Path to output CSV file (default: data/processed/{species}/input_filename_sequences.csv)')
    
    args = parser.parse_args()
    
    # Determine output path
    if args.output:
        output_path = args.output
    else:
        # Create output directory
        output_dir = os.path.join('data', 'processed', args.species)
        os.makedirs(output_dir, exist_ok=True)
        
        # Get input filename and create sequences filename
        input_filename = os.path.basename(args.input)
        base_name = os.path.splitext(input_filename)[0]
        output_filename = f"{base_name}_sequences.csv"
        output_path = os.path.join(output_dir, output_filename)
    
    # Process the dataset
    process_dataset(args.input, args.species, output_path)


if __name__ == "__main__":
    main()

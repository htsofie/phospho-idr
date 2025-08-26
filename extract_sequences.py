#!/usr/bin/env python3
"""
Script to extract sequences from test_af.tsv and create separate files for each sequence
with their corresponding uniprot IDs.
"""

import csv
import os
import re
import json
from pathlib import Path

def extract_sequences_from_tsv(tsv_file_path, output_dir="sequence_files"):
    """
    Extract sequences from TSV file and create a single JSON file with all jobs
    in AlphaFold server format with their corresponding uniprot IDs.
    
    Args:
        tsv_file_path (str): Path to the input TSV file
        output_dir (str): Directory to save the sequence files
    """
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(exist_ok=True)
    
    # List to store all jobs
    all_jobs = []
    
    # Read the TSV file
    with open(tsv_file_path, 'r', encoding='utf-8') as file:
        reader = csv.DictReader(file, delimiter='\t')
        
        for row_num, row in enumerate(reader, start=1):
            # Extract the required fields
            uniprot_1 = row['uniprot_1']
            uniprot_2 = row['uniprot_2']
            sequence = row['sequence']
            job_name = row['job_name']
            
            # Split sequence at colon and separate into two parts
            sequence_parts = sequence.split(':')
            sequence_part_1 = sequence_parts[0] if len(sequence_parts) > 0 else ""
            sequence_part_2 = sequence_parts[1] if len(sequence_parts) > 1 else ""
            
            # Create AlphaFold job format
            af3_job = {
                "name": job_name,
                "modelSeeds": [],
                "sequences": [
                    {
                        "proteinChain": {
                            "sequence": sequence_part_1
                        }
                    },
                    {
                        "proteinChain": {
                            "sequence": sequence_part_2
                        }
                    }
                ],
                "dialect": "alphafoldserver",
                "version": 1
            }
            
            # Add job to the list
            all_jobs.append(af3_job)
            
            print(f"Added job: {job_name}")
            print(f"  UniProt 1: {uniprot_1}")
            print(f"  UniProt 2: {uniprot_2}")
            print(f"  Sequence length: {len(sequence)} characters")
            print()
    
    # Write all jobs to a single JSON file
    output_file = os.path.join(output_dir, "alphafold_jobs.json")
    with open(output_file, 'w', encoding='utf-8') as json_file:
        json.dump(all_jobs, json_file, indent=2)
    
    print(f"Created file: {output_file}")
    print(f"Total jobs: {len(all_jobs)}")

def main():
    """Main function to run the sequence extraction."""
    
    # Path to the TSV file
    tsv_file = "test_af.tsv"
    
    # Check if the file exists
    if not os.path.exists(tsv_file):
        print(f"Error: File '{tsv_file}' not found!")
        print("Please make sure the script is run from the same directory as the TSV file.")
        return
    
    print(f"Processing file: {tsv_file}")
    print("=" * 60)
    
    # Extract sequences
    extract_sequences_from_tsv(tsv_file)
    
    print("=" * 60)
    print("Sequence extraction completed!")
    print("Check the 'sequence_files' directory for the output files.")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Grouped integrated pipeline for protein ID mapping after BLAST.
Works on rows where manual_review == TRUE.
Searches all IDs for each protein group and extracts sequences.
Maintains ID_matches and full_sequence as comma-separated lists in same order.
"""

import pandas as pd
import argparse
import os
import time
from typing import Dict, Optional, Tuple, List, Any
from Bio import ExPASy
from Bio import SwissProt
from Bio.Align import PairwiseAligner
import requests

# Excel cell limit (32,767 characters)
EXCEL_CELL_LIMIT = 32767


def get_uniprot_organism(uniprot_id: str) -> Optional[str]:
    """Get organism information for a UniProt ID using REST API."""
    try:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            organism = data.get('organism', {})
            scientific_name = organism.get('scientificName', '')
            return scientific_name.lower()
        return None
    except Exception as e:
        print(f"  ✗ Failed to get organism for {uniprot_id}: {e}")
        return None


def verify_species_match(uniprot_id: str, expected_species: str) -> bool:
    """Verify if the mapped UniProt ID belongs to the expected species."""
    if not uniprot_id:
        return False
    
    organism = get_uniprot_organism(uniprot_id)
    if not organism:
        return False
    
    # Map expected species to organism names
    species_mapping = {
        'mouse': ['mus musculus', 'mus'],
        'rat': ['rattus norvegicus', 'rattus']
    }
    
    expected_organisms = species_mapping.get(expected_species.lower(), [])
    return any(exp_org in organism for exp_org in expected_organisms)


def map_uniprot_id_to_uniprot(uniprot_id: str) -> Optional[str]:
    """Map UniProt ID to verified UniProt ID."""
    if pd.isna(uniprot_id) or not uniprot_id:
        return None
    
    try:
        clean_id = str(uniprot_id).strip()
        if ';' in clean_id:
            ids = [id.strip() for id in clean_id.split(';') if id.strip()]
            for single_id in ids:
                result = _verify_single_uniprot_id(single_id)
                if result:
                    return result
            return None
        else:
            return _verify_single_uniprot_id(clean_id)
    except Exception:
        return None


def _verify_single_uniprot_id(uniprot_id: str) -> Optional[str]:
    """Verify a single UniProt ID exists."""
    try:
        handle = ExPASy.get_sprot_raw(uniprot_id)
        record = SwissProt.read(handle)
        handle.close()
        return uniprot_id
    except Exception:
        return None


def map_swissprot_to_uniprot(swissprot_id: str) -> Optional[str]:
    """Map SwissProt ID to UniProt ID."""
    if pd.isna(swissprot_id) or not swissprot_id:
        return None
    
    try:
        clean_id = str(swissprot_id).strip()
        if ';' in clean_id:
            ids = [id.strip() for id in clean_id.split(';') if id.strip()]
            for single_id in ids:
                result = _verify_single_uniprot_id(single_id)
                if result:
                    return result
            return None
        else:
            return _verify_single_uniprot_id(clean_id)
    except Exception:
        return None


def map_trembl_to_uniprot(trembl_id: str) -> Optional[str]:
    """Map TrEMBL ID to UniProt ID."""
    if pd.isna(trembl_id) or not trembl_id:
        return None
    
    try:
        clean_id = str(trembl_id).strip()
        if ';' in clean_id:
            ids = [id.strip() for id in clean_id.split(';') if id.strip()]
            for single_id in ids:
                result = _verify_single_uniprot_id(single_id)
                if result:
                    return result
            return None
        else:
            return _verify_single_uniprot_id(clean_id)
    except Exception:
        return None


def get_uniprot_sequence(uniprot_id: str) -> Optional[Tuple[str, str, str]]:
    """Get full protein sequence from UniProt and determine db_type."""
    if pd.isna(uniprot_id) or not uniprot_id:
        return None
    
    try:
        clean_id = str(uniprot_id).strip()
        handle = ExPASy.get_sprot_raw(clean_id)
        record = SwissProt.read(handle)
        handle.close()
        
        sequence = record.sequence
        
        # Determine db_type: SwissProt entries have reviewed status
        # For simplicity, check entry name format or use data_class
        db_type = 'swissprot' if record.data_class == 'STANDARD' or record.data_class == 'Reviewed' else 'trembl'
        
        if '-' in clean_id and clean_id.split('-')[-1].isdigit():
            return sequence, "isoform", db_type
        else:
            return sequence, "canonical", db_type
    except Exception:
        return None


def collect_all_ids_from_group(group_df: pd.DataFrame) -> List[str]:
    """Collect all unique IDs from all rows in a protein group."""
    all_ids = []
    id_columns = ['uniprot', 'uniprot_id', 'swissprot_id', 'trembl_id', 'uniprot_mapped_id']
    
    for idx, row in group_df.iterrows():
        for col in id_columns:
            if col in row and pd.notna(row[col]) and row[col] != '':
                ids = str(row[col]).split(';')
                all_ids.extend([id.strip() for id in ids if id.strip()])
    
    # Remove duplicates while preserving order
    unique_ids = list(dict.fromkeys(all_ids))
    return unique_ids


def validate_sequence_length(sequence: str, phospho_position: int) -> bool:
    """Check if sequence length is greater than phosphorylation position."""
    if not sequence or pd.isna(phospho_position):
        return False
    return len(sequence) > phospho_position




def check_cell_limit_and_sort(ids: List[str], sequences: List[str], lengths: List[int], 
                              db_types: List[str]) -> Tuple[List[str], List[str], List[int]]:
    """Sort IDs (SwissProt first, then TrEMBL) and ensure results fit within Excel cell limits."""
    if not ids:
        return [], [], []
    
    # Create tuples for sorting
    combined = list(zip(ids, sequences, lengths, db_types))
    
    # Sort: SwissProt first, then TrEMBL
    combined.sort(key=lambda x: (0 if x[3] == 'swissprot' else 1, x[0]))
    
    # Unpack sorted data
    sorted_ids, sorted_seqs, sorted_lengths, _ = zip(*combined)
    sorted_ids = list(sorted_ids)
    sorted_seqs = list(sorted_seqs)
    sorted_lengths = list(sorted_lengths)
    
    # Check if all entries fit within cell limit
    id_str = ';'.join(sorted_ids)
    seq_str = ';'.join(sorted_seqs)
    length_str = ';'.join(map(str, sorted_lengths))
    
    max_len = max(len(id_str), len(seq_str), len(length_str))
    
    if max_len <= EXCEL_CELL_LIMIT:
        print(f"    Using all {len(sorted_ids)} matches (fits within cell limit)")
        return sorted_ids, sorted_seqs, sorted_lengths
    
    # If doesn't fit, reduce from end one by one
    print(f"    Warning: {len(sorted_ids)} matches exceed cell limit, reducing...")
    for i in range(len(sorted_ids) - 1, 0, -1):
        test_ids = sorted_ids[:i]
        test_seqs = sorted_seqs[:i]
        test_lengths = sorted_lengths[:i]
        
        id_str = ';'.join(test_ids)
        seq_str = ';'.join(test_seqs)
        length_str = ';'.join(map(str, test_lengths))
        
        max_len = max(len(id_str), len(seq_str), len(length_str))
        
        if max_len <= EXCEL_CELL_LIMIT:
            print(f"    Using {i} matches to fit within cell limit")
            return test_ids, test_seqs, test_lengths
    
    # If even 1 doesn't fit, return empty
    print(f"    Cannot fit even 1 match within cell limit")
    return [], [], []


def process_protein_group(group_df: pd.DataFrame, group_idx: int, total_groups: int, 
                         species: str, sequence_cache: Dict[str, Tuple[str, str, str]]) -> Dict[str, Any]:
    """Process a protein group by searching all IDs and extracting sequences."""
    protein_id = group_df.iloc[0]['Protein']
    group_size = len(group_df)
    
    print(f"\nGroup {group_idx + 1}/{total_groups} (Protein: {protein_id}, {group_size} phosphosites)")
    
    # Initialize result for all rows in group
    group_result = {
        'ID_matches': '',
        'full_sequence': '',
        'length': '',
        'manual_review': True,  # Will be set to False if we find valid results
        'match_method': '',
        'processing_notes': []
    }
    
    # Collect all IDs from all rows in the group
    all_ids = collect_all_ids_from_group(group_df)
    
    if not all_ids:
        group_result['processing_notes'].append("No IDs found in group")
        print(f"  ✗ No IDs found in group")
        return group_result
    
    print(f"  Found {len(all_ids)} unique IDs to search: {all_ids[:5]}{'...' if len(all_ids) > 5 else ''}")
    
    # Try to verify and fetch sequences for all IDs
    matched_ids = []
    matched_sequences = []
    matched_lengths = []
    matched_db_types = []
    
    for uniprot_id in all_ids:
        # Verify ID exists
        verified_id = _verify_single_uniprot_id(uniprot_id)
        if not verified_id:
            print(f"    ✗ {uniprot_id}: Not found in UniProt")
            continue
        
        # Verify species match
        if not verify_species_match(verified_id, species):
            print(f"    ✗ {verified_id}: Species mismatch")
            continue
        
        # Check cache first
        if verified_id in sequence_cache:
            sequence, seq_type, db_type = sequence_cache[verified_id]
            print(f"    ✓ {verified_id}: Cached ({len(sequence)} aa, {seq_type}, {db_type})")
        else:
            # Get sequence from UniProt
            sequence_data = get_uniprot_sequence(verified_id)
            if not sequence_data:
                print(f"    ✗ {verified_id}: Failed to fetch sequence")
                continue
            
            sequence, seq_type, db_type = sequence_data
            sequence_cache[verified_id] = (sequence, seq_type, db_type)
            print(f"    ✓ {verified_id}: Fetched ({len(sequence)} aa, {seq_type}, {db_type})")
        
        # Add to matched lists
        matched_ids.append(verified_id)
        matched_sequences.append(sequence)
        matched_lengths.append(len(sequence))
        matched_db_types.append(db_type)
        
        # Small delay to be respectful to API
        time.sleep(0.05)
    
    if matched_ids:
        # Sort by db_type and check cell limits
        limited_ids, limited_seqs, limited_lengths = check_cell_limit_and_sort(
            matched_ids, matched_sequences, matched_lengths, matched_db_types
        )
        
        if limited_ids:
            # Update results with semicolon-separated lists
            group_result['ID_matches'] = ';'.join(limited_ids)
            group_result['full_sequence'] = ';'.join(limited_seqs)
            group_result['length'] = ';'.join(map(str, limited_lengths))
            group_result['manual_review'] = False  # Found valid results, clear the flag
            group_result['match_method'] = 'uniprot_search'
            group_result['processing_notes'].append(f"Found {len(limited_ids)} matching IDs with sequences")
            print(f"  ✓ Success: Found {len(limited_ids)} matching IDs with sequences")
        else:
            group_result['processing_notes'].append("Results exceed Excel cell limit")
            print(f"  ✗ Results exceed Excel cell limit")
    else:
        group_result['processing_notes'].append("No valid IDs found with species match and sequences")
        print(f"  ✗ No valid IDs found with species match and sequences")
    
    return group_result


def apply_group_result_to_all_rows(df: pd.DataFrame, group_index: pd.Index, group_result: Dict[str, Any]) -> None:
    """Apply the group result to all rows in the main dataframe for the given group index."""
    for key, value in group_result.items():
        # For each row in the group, assign the same value
        for idx in group_index:
            df.at[idx, key] = value


def process_dataset(input_path: str, species: str, output_path: str) -> None:
    """Process the entire dataset - only groups with manual_review == TRUE."""
    print(f"Loading dataset: {input_path}")
    df = pd.read_csv(input_path)
    print(f"Loaded {len(df)} rows")
    
    # Check if required columns exist
    if 'Protein' not in df.columns:
        print("Error: 'Protein' column not found. This pipeline requires protein grouping.")
        return
    
    if 'manual_review' not in df.columns:
        print("Error: 'manual_review' column not found. This script should run after blast.py.")
        return
    
    # Filter to only manual_review == TRUE
    manual_review_mask = df['manual_review'] == True
    rows_to_process = df[manual_review_mask]
    print(f"Found {len(rows_to_process)} rows with manual_review == TRUE")
    
    if len(rows_to_process) == 0:
        print("No rows to process. Exiting.")
        return
    
    # Group by Protein column (only the filtered rows)
    protein_groups = rows_to_process.groupby('Protein')
    total_groups = len(protein_groups)
    print(f"Found {total_groups} unique proteins to process")
    
    # Process each protein group
    sequence_cache = {}
    successful_groups = 0
    failed_groups = 0
    
    for group_idx, (protein_id, group_df) in enumerate(protein_groups):
        # Process the group
        group_result = process_protein_group(group_df, group_idx, total_groups, species, sequence_cache)
        
        # Apply result to all rows in the main dataframe for this group's indices
        # Update existing columns
        for idx in group_df.index:
            if group_result['ID_matches']:
                df.at[idx, 'ID_matches'] = group_result['ID_matches']
            if group_result['full_sequence']:
                df.at[idx, 'full_sequence'] = group_result['full_sequence']
            if group_result['length']:
                df.at[idx, 'length'] = group_result['length']
            # Update manual_review and match_method if we found results
            if group_result['ID_matches']:
                df.at[idx, 'manual_review'] = False
                df.at[idx, 'match_method'] = 'uniprot_search'
        
        # Update statistics
        if group_result['ID_matches']:
            successful_groups += 1
        else:
            failed_groups += 1
        
        # Add delay to be respectful to APIs
        time.sleep(0.1)
    
    # Remove error_reason column before saving
    if 'error_reason' in df.columns:
        df = df.drop(columns=['error_reason'])
        print(f"Removed 'error_reason' column from output")
    
    # Save results
    print(f"\nSaving results to: {output_path}")
    df.to_csv(output_path, index=False)
    
    # Print summary
    print(f"\nProcessing Summary:")
    print(f"  Total rows in dataset: {len(df)}")
    print(f"  Rows with manual_review == TRUE: {len(rows_to_process)}")
    print(f"  Total protein groups processed: {total_groups}")
    print(f"  Successful groups (found IDs & sequences): {successful_groups}")
    print(f"  Failed groups: {failed_groups}")
    print(f"  Sequence cache size: {len(sequence_cache)}")
    



def main():
    parser = argparse.ArgumentParser(description='Grouped integrated pipeline for protein processing')
    parser.add_argument('--input', '-i', required=True, help='Path to input CSV file')
    parser.add_argument('--species', '-s', required=True, choices=['mouse', 'rat'], help='Species name')
    parser.add_argument('--output', '-o', help='Path to output CSV file (default: data/processed/{species}/input_filename_grouped_processed.csv)')
    
    args = parser.parse_args()
    
    # Determine output path
    if args.output:
        output_path = args.output
    else:
        output_dir = os.path.join('data', 'processed', args.species)
        os.makedirs(output_dir, exist_ok=True)
        input_filename = os.path.basename(args.input)
        base_name = os.path.splitext(input_filename)[0]
        output_filename = f"{base_name}_grouped_processed.csv"
        output_path = os.path.join(output_dir, output_filename)
    
    process_dataset(args.input, args.species, output_path)


if __name__ == "__main__":
    main()

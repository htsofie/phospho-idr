#!/usr/bin/env python3
"""
Integrated pipeline for protein ID mapping, sequence extraction, and validation.
Processes each row individually with retry logic and cross-row sequence sharing.
"""

import pandas as pd
import argparse
import os
import time
from typing import Dict, Optional, Tuple, List
from Bio import ExPASy
from Bio import SwissProt
from Bio.Align import PairwiseAligner
import requests


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


def map_ensembl_to_uniprot(ensembl_id: str) -> Optional[str]:
    """Map Ensembl ID to UniProt ID using REST API."""
    if pd.isna(ensembl_id) or not ensembl_id:
        return None
    
    try:
        clean_id = str(ensembl_id).strip()
        if ';' in clean_id:
            ids = [id.strip() for id in clean_id.split(';') if id.strip()]
            for single_id in ids:
                result = _map_single_ensembl_to_uniprot(single_id)
                if result:
                    return result
            return None
        else:
            return _map_single_ensembl_to_uniprot(clean_id)
    except Exception:
        return None


def _map_single_ensembl_to_uniprot(ensembl_id: str) -> Optional[str]:
    """Map a single Ensembl ID to UniProt ID."""
    try:
        url = f"https://rest.ensembl.org/xrefs/id/{ensembl_id}"
        headers = {"Content-Type": "application/json"}
        response = requests.get(url, headers=headers)
        
        if response.status_code == 200:
            data = response.json()
            for entry in data:
                if entry.get('dbname') in ['Uniprot/SWISSPROT', 'Uniprot/SPTREMBL', 'UniProt']:
                    uniprot_id = entry.get('primary_id')
                    if uniprot_id:
                        return uniprot_id
        return None
    except Exception:
        return None


def get_uniprot_sequence(uniprot_id: str) -> Optional[Tuple[str, str]]:
    """Get full protein sequence from UniProt."""
    if pd.isna(uniprot_id) or not uniprot_id:
        return None
    
    try:
        clean_id = str(uniprot_id).strip()
        handle = ExPASy.get_sprot_raw(clean_id)
        record = SwissProt.read(handle)
        handle.close()
        
        sequence = record.sequence
        if '-' in clean_id and clean_id.split('-')[-1].isdigit():
            return sequence, "isoform"
        else:
            return sequence, "canonical"
    except Exception:
        return None


def try_id_mapping(row: pd.Series, mapping_priority: List[str]) -> Optional[Tuple[str, str]]:
    """Try ID mapping in priority order."""
    for mapping_type in mapping_priority:
        if mapping_type == 'uniprot_id' and not pd.isna(row.get('uniprot_id')):
            mapped_id = map_uniprot_id_to_uniprot(row['uniprot_id'])
            if mapped_id:
                return mapped_id, 'uniprot_id'
        
        elif mapping_type == 'swissprot_id' and not pd.isna(row.get('swissprot_id')):
            mapped_id = map_swissprot_to_uniprot(row['swissprot_id'])
            if mapped_id:
                return mapped_id, 'swissprot_id'
        
        elif mapping_type == 'trembl_id' and not pd.isna(row.get('trembl_id')):
            mapped_id = map_trembl_to_uniprot(row['trembl_id'])
            if mapped_id:
                return mapped_id, 'trembl_id'
        
        elif mapping_type == 'ensembl_id' and not pd.isna(row.get('ensembl_id')):
            mapped_id = map_ensembl_to_uniprot(row['ensembl_id'])
            if mapped_id:
                return mapped_id, 'ensembl_id'
    
    return None, None


def validate_sequence_length(sequence: str, phospho_position: int) -> bool:
    """Check if sequence length is greater than phosphorylation position."""
    if not sequence or pd.isna(phospho_position):
        return False
    return len(sequence) > phospho_position


def process_single_row(row: pd.Series, row_idx: int, total_rows: int, 
                      sequence_cache: Dict[str, Tuple[str, str, str]]) -> Dict[str, any]:
    """Process a single row with retry logic."""
    print(f"\nRow {row_idx + 1}/{total_rows}")
    print(f"  Position: {row.get('position', 'N/A')}")
    
    result = {
        'uniprot_mapped_id': None,
        'mapping_source': None,
        'mapping_success': False,
        'full_sequence': None,
        'sequence_length': None,
        'sequence_type': None,
        'sequence_fetch_success': False,
        'sequence_length_warning': None,
        'manual_review_flag': False,
        'processing_notes': []
    }
    
    # Get phosphorylation position
    phospho_position = None
    if not pd.isna(row.get('position')) and row.get('position'):
        try:
            phospho_position = int(row['position'])
        except (ValueError, TypeError):
            pass
    
    if not phospho_position:
        result['manual_review_flag'] = True
        result['processing_notes'].append("No valid phosphorylation position")
        print(f"  ✗ No valid phosphorylation position")
        return result
    
    # Define mapping priority
    mapping_priority = ['uniprot_id', 'swissprot_id', 'trembl_id', 'ensembl_id']
    
    # Try mapping and sequence extraction with retry logic
    for attempt in range(len(mapping_priority)):
        print(f"  Attempt {attempt + 1}: Trying {mapping_priority[attempt]}")
        
        # Try ID mapping
        mapped_id, mapping_source = try_id_mapping(row, mapping_priority[attempt:])
        
        if not mapped_id:
            result['processing_notes'].append(f"Failed to map {mapping_priority[attempt]}")
            continue
        
        print(f"    Mapped to: {mapped_id} (via {mapping_source})")
        
        # Check if we already have this sequence in cache
        if mapped_id in sequence_cache:
            sequence, seq_type, source = sequence_cache[mapped_id]
            print(f"    Using cached sequence: {len(sequence)} amino acids")
        else:
            # Get sequence from UniProt
            sequence_data = get_uniprot_sequence(mapped_id)
            if not sequence_data:
                result['processing_notes'].append(f"Failed to fetch sequence for {mapped_id}")
                continue
            
            sequence, seq_type = sequence_data
            sequence_cache[mapped_id] = (sequence, seq_type, mapping_source)
            print(f"    Fetched sequence: {len(sequence)} amino acids ({seq_type})")
        
        # Validate sequence length
        if validate_sequence_length(sequence, phospho_position):
            # Success!
            result['uniprot_mapped_id'] = mapped_id
            result['mapping_source'] = mapping_source
            result['mapping_success'] = True
            result['full_sequence'] = sequence
            result['sequence_length'] = len(sequence)
            result['sequence_type'] = seq_type
            result['sequence_fetch_success'] = True
            result['processing_notes'].append(f"Success with {mapping_source}")
            print(f"  ✓ Success: {mapping_source} -> {mapped_id}, {len(sequence)} amino acids")
            return result
        else:
            warning = f"Sequence length ({len(sequence)}) <= phosphorylation position ({phospho_position})"
            result['processing_notes'].append(f"{mapping_source}: {warning}")
            print(f"    ⚠ {warning}")
    
    # All attempts failed
    result['manual_review_flag'] = True
    result['processing_notes'].append("All mapping attempts failed or sequence too short")
    print(f"  ✗ All attempts failed - flagged for manual review")
    return result


def apply_sequence_to_matching_rows(df: pd.DataFrame, successful_row_idx: int, 
                                  sequence_cache: Dict[str, Tuple[str, str, str]]) -> None:
    """Apply successful sequence to other rows with matching identifiers."""
    successful_row = df.iloc[successful_row_idx]
    successful_id = successful_row['uniprot_mapped_id']
    successful_sequence, successful_type, successful_source = sequence_cache[successful_id]
    
    # Find rows with matching identifiers that don't have sequences yet
    for idx, row in df.iterrows():
        if idx == successful_row_idx or row['sequence_fetch_success']:
            continue
        
        # Check if this row has any matching identifiers
        matching = False
        for col in ['uniprot_id', 'swissprot_id', 'trembl_id', 'ensembl_id']:
            if not pd.isna(row.get(col)) and not pd.isna(successful_row.get(col)):
                if str(row[col]).strip() == str(successful_row[col]).strip():
                    matching = True
                    break
        
        if matching:
            # Apply the successful sequence
            df.at[idx, 'uniprot_mapped_id'] = successful_id
            df.at[idx, 'mapping_source'] = f"{successful_source}_shared"
            df.at[idx, 'mapping_success'] = True
            df.at[idx, 'full_sequence'] = successful_sequence
            df.at[idx, 'sequence_length'] = len(successful_sequence)
            df.at[idx, 'sequence_type'] = successful_type
            df.at[idx, 'sequence_fetch_success'] = True
            df.at[idx, 'processing_notes'] = f"Applied sequence from row {successful_row_idx + 1}"
            print(f"    Applied sequence to row {idx + 1} (matching identifiers)")


def process_dataset(input_path: str, species: str, output_path: str) -> None:
    """Process the entire dataset with integrated pipeline."""
    print(f"Loading dataset: {input_path}")
    df = pd.read_csv(input_path)
    print(f"Loaded {len(df)} rows")
    
    # Initialize result columns
    result_columns = [
        'uniprot_mapped_id', 'mapping_source', 'mapping_success',
        'full_sequence', 'sequence_length', 'sequence_type', 'sequence_fetch_success',
        'sequence_length_warning', 'manual_review_flag', 'processing_notes'
    ]
    
    for col in result_columns:
        df[col] = None
    
    # Process each row
    sequence_cache = {}
    successful_rows = []
    
    for idx, row in df.iterrows():
        result = process_single_row(row, idx, len(df), sequence_cache)
        
        # Update dataframe with results
        for key, value in result.items():
            df.at[idx, key] = value
        
        # If successful, apply to matching rows
        if result['sequence_fetch_success']:
            successful_rows.append(idx)
            apply_sequence_to_matching_rows(df, idx, sequence_cache)
        
        # Add delay to be respectful to APIs
        time.sleep(0.1)
    
    # Reorder columns
    original_cols = [col for col in df.columns if col not in result_columns]
    df_reordered = df[original_cols + result_columns]
    
    # Save results
    print(f"\nSaving results to: {output_path}")
    df_reordered.to_csv(output_path, index=False)
    
    # Print summary
    successful = df['sequence_fetch_success'].sum()
    manual_review = df['manual_review_flag'].sum()
    
    print(f"\nProcessing Summary:")
    print(f"  Total rows: {len(df)}")
    print(f"  Successful: {successful} ({successful/len(df)*100:.1f}%)")
    print(f"  Manual review needed: {manual_review} ({manual_review/len(df)*100:.1f}%)")
    print(f"  Sequence cache size: {len(sequence_cache)}")
    
    # Show manual review rows
    if manual_review > 0:
        print(f"\nRows flagged for manual review:")
        manual_rows = df[df['manual_review_flag'] == True]
        for idx, row in manual_rows.iterrows():
            print(f"  Row {idx + 1}: {row.get('processing_notes', 'No notes')}")


def main():
    parser = argparse.ArgumentParser(description='Integrated pipeline for protein processing')
    parser.add_argument('--input', '-i', required=True, help='Path to input CSV file')
    parser.add_argument('--species', '-s', required=True, choices=['mouse', 'rat'], help='Species name')
    parser.add_argument('--output', '-o', help='Path to output CSV file (default: data/processed/{species}/input_filename_processed.csv)')
    
    args = parser.parse_args()
    
    # Determine output path
    if args.output:
        output_path = args.output
    else:
        output_dir = os.path.join('data', 'processed', args.species)
        os.makedirs(output_dir, exist_ok=True)
        input_filename = os.path.basename(args.input)
        base_name = os.path.splitext(input_filename)[0]
        output_filename = f"{base_name}_processed.csv"
        output_path = os.path.join(output_dir, output_filename)
    
    process_dataset(args.input, args.species, output_path)


if __name__ == "__main__":
    main()

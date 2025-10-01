#!/usr/bin/env python3
"""
Grouped integrated pipeline for protein ID mapping, sequence extraction, and validation.
Groups phosphosites by protein and processes each group as a unit.
If species mismatch occurs for the first row of a group, flags all rows in that group.
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
        if mapping_type == 'uniprot' and not pd.isna(row.get('uniprot')):
            mapped_id = map_uniprot_id_to_uniprot(row['uniprot'])
            if mapped_id:
                return mapped_id, 'uniprot'
        
        elif mapping_type == 'swissprot_id' and not pd.isna(row.get('swissprot_id')):
            mapped_id = map_swissprot_to_uniprot(row['swissprot_id'])
            if mapped_id:
                return mapped_id, 'swissprot_id'
        
        elif mapping_type == 'trembl_id' and not pd.isna(row.get('trembl_id')):
            mapped_id = map_trembl_to_uniprot(row['trembl_id'])
            if mapped_id:
                return mapped_id, 'trembl_id'
    
    return None, None


def validate_sequence_length(sequence: str, phospho_position: int) -> bool:
    """Check if sequence length is greater than phosphorylation position."""
    if not sequence or pd.isna(phospho_position):
        return False
    return len(sequence) > phospho_position




def process_protein_group(group_df: pd.DataFrame, group_idx: int, total_groups: int, 
                         species: str, sequence_cache: Dict[str, Tuple[str, str, str]]) -> Dict[str, any]:
    """Process a protein group by trying to map the first row."""
    protein_id = group_df.iloc[0]['Protein']
    group_size = len(group_df)
    
    print(f"\nGroup {group_idx + 1}/{total_groups} (Protein: {protein_id}, {group_size} phosphosites)")
    
    # Initialize result for all rows in group
    group_result = {
        'uniprot_mapped_id': None,
        'mapping_source': None,
        'full_sequence': None,
        'sequence_length': None,
        'sequence_type': None,
        'manual_review_flag': False,
        'processing_notes': []
    }
    
    # Get first row for mapping attempt
    first_row = group_df.iloc[0]
    
    # Get phosphorylation position from first row
    phospho_position = None
    if not pd.isna(first_row.get('position')) and first_row.get('position'):
        try:
            phospho_position = int(first_row['position'])
        except (ValueError, TypeError):
            pass
    
    if not phospho_position:
        group_result['manual_review_flag'] = True
        group_result['processing_notes'].append("No valid phosphorylation position in first row")
        print(f"  ✗ No valid phosphorylation position in first row")
        return group_result
    
    # Define mapping priority (only UniProt, SWISS-PROT, TREMBL)
    mapping_priority = ['uniprot', 'swissprot_id', 'trembl_id']
    
    # Try mapping and sequence extraction
    for attempt in range(len(mapping_priority)):
        print(f"  Attempt {attempt + 1}: Trying {mapping_priority[attempt]}")
        
        # Try ID mapping
        mapped_id, mapping_source = try_id_mapping(first_row, mapping_priority[attempt:])
        
        if not mapped_id:
            group_result['processing_notes'].append(f"Failed to map {mapping_priority[attempt]}")
            continue
        
        print(f"    Mapped to: {mapped_id} (via {mapping_source})")
        
        # Verify species match
        if not verify_species_match(mapped_id, species):
            group_result['processing_notes'].append(f"Species mismatch for {mapped_id}")
            print(f"    ⚠ Species mismatch for {mapped_id}, trying next source...")
            continue
        
        print(f"    ✓ Species verified for {mapped_id}")
        
        # Check if we already have this sequence in cache
        if mapped_id in sequence_cache:
            sequence, seq_type, source = sequence_cache[mapped_id]
            print(f"    Using cached sequence: {len(sequence)} amino acids")
        else:
            # Get sequence from UniProt
            sequence_data = get_uniprot_sequence(mapped_id)
            if not sequence_data:
                group_result['processing_notes'].append(f"Failed to fetch sequence for {mapped_id}")
                continue
            
            sequence, seq_type = sequence_data
            sequence_cache[mapped_id] = (sequence, seq_type, mapping_source)
            print(f"    Fetched sequence: {len(sequence)} amino acids ({seq_type})")
        
        # Validate sequence length
        if validate_sequence_length(sequence, phospho_position):
            # Success! Set basic mapping results
            group_result['uniprot_mapped_id'] = mapped_id
            group_result['mapping_source'] = mapping_source
            group_result['full_sequence'] = sequence
            group_result['sequence_length'] = len(sequence)
            group_result['sequence_type'] = seq_type
            group_result['processing_notes'].append(f"Success with {mapping_source}")
            print(f"  ✓ Success: {mapping_source} -> {mapped_id}, {len(sequence)} amino acids")
            return group_result
        else:
            warning = f"Sequence length ({len(sequence)}) <= phosphorylation position ({phospho_position})"
            group_result['processing_notes'].append(f"{mapping_source}: {warning}")
            print(f"    ⚠ {warning}")
    
    # All attempts failed
    group_result['manual_review_flag'] = True
    group_result['processing_notes'].append("All mapping attempts failed or had species mismatches")
    print(f"  ⚠ All attempts failed or had species mismatches - flagged for manual review")
    
    return group_result


def apply_group_result_to_all_rows(df: pd.DataFrame, group_index: pd.Index, group_result: Dict[str, Any]) -> None:
    """Apply the group result to all rows in the main dataframe for the given group index."""
    for key, value in group_result.items():
        # For each row in the group, assign the same value
        for idx in group_index:
            df.at[idx, key] = value


def process_dataset(input_path: str, species: str, output_path: str) -> None:
    """Process the entire dataset with grouped pipeline."""
    print(f"Loading dataset: {input_path}")
    df = pd.read_csv(input_path)
    print(f"Loaded {len(df)} rows")
    
    # Check if Protein column exists
    if 'Protein' not in df.columns:
        print("Error: 'Protein' column not found. This pipeline requires protein grouping.")
        return
    
    # Initialize result columns
    result_columns = [
        'uniprot_mapped_id', 'mapping_source', 'full_sequence', 'sequence_length', 'sequence_type',
        'manual_review_flag', 'processing_notes'
    ]
    
    for col in result_columns:
        df[col] = None
    
    # Group by Protein column
    protein_groups = df.groupby('Protein')
    total_groups = len(protein_groups)
    print(f"Found {total_groups} unique proteins")
    
    # Process each protein group
    sequence_cache = {}
    successful_groups = 0
    failed_groups = 0
    
    for group_idx, (protein_id, group_df) in enumerate(protein_groups):
        # Process the group
        group_result = process_protein_group(group_df, group_idx, total_groups, species, sequence_cache)
        
        # Apply result to all rows in the main dataframe for this group's indices
        apply_group_result_to_all_rows(df, group_df.index, group_result)
        
        # Update statistics
        if group_result['uniprot_mapped_id'] is not None:
            successful_groups += 1
        else:
            failed_groups += 1
        
        # Add delay to be respectful to APIs
        time.sleep(0.1)
    
    # Reorder columns
    original_cols = [col for col in df.columns if col not in result_columns]
    df_reordered = df[original_cols + result_columns]
    
    # Save results
    print(f"\nSaving results to: {output_path}")
    df_reordered.to_csv(output_path, index=False)
    
    # Print summary
    successful_rows = df['uniprot_mapped_id'].notna().sum()
    manual_review_rows = df['manual_review_flag'].sum()
    
    print(f"\nProcessing Summary:")
    print(f"  Total rows: {len(df)}")
    print(f"  Total protein groups: {total_groups}")
    print(f"  Successful groups: {successful_groups}")
    print(f"  Failed groups: {failed_groups}")
    print(f"  Successful rows: {successful_rows} ({successful_rows/len(df)*100:.1f}%)")
    print(f"  Manual review rows: {manual_review_rows} ({manual_review_rows/len(df)*100:.1f}%)")
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

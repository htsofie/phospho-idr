#!/usr/bin/env python3
"""
Match MoRF residues to proteins in human_ortholog_mapped.csv.

This script:
1. Parses CV1-CV4 FASTA files to extract MoRF data
2. Maps Disprot/Ideal IDs to UniProt using UniProt ID Mapping API
3. Matches MoRF residues to proteins in human_ortholog_mapped.csv
4. Generates matched_morfs.csv with matched data

Usage:
    python scripts/match_morf_interfaces.py --input data/processed/mouse/human_ortholog_mapped.csv --species mouse
    python scripts/match_morf_interfaces.py --input data/processed/rat/human_ortholog_mapped.csv --species rat
"""

import pandas as pd
import requests
import time
import argparse
import logging
import os
from pathlib import Path
from typing import Optional, List, Dict
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# API rate limiting delay (seconds)
API_DELAY = 0.15


def extract_morf_positions(annotation_line: str, sequence: str) -> str:
    """
    Parse annotation line and convert '1' positions to range format.
    
    Args:
        annotation_line: String of 0, 1, - characters
        sequence: Protein sequence (for validation)
        
    Returns:
        Formatted string: [20-30,40,50-60] (1-indexed positions)
    """
    positions = []
    
    # Find all positions where character is '1' (1-indexed)
    for i, char in enumerate(annotation_line, start=1):
        if char == '1':
            positions.append(i)
    
    if not positions:
        return '[]'
    
    # Group consecutive positions into ranges
    ranges = []
    start = positions[0]
    end = positions[0]
    
    for i in range(1, len(positions)):
        if positions[i] == end + 1:
            # Consecutive, extend range
            end = positions[i]
        else:
            # Gap found, save current range
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{end}")
            start = positions[i]
            end = positions[i]
    
    # Add last range
    if start == end:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}-{end}")
    
    return '[' + ','.join(ranges) + ']'


def classify_database_type(protein_id: str) -> str:
    """
    Classify ID as uniprot, disprot, or ideal.
    
    Args:
        protein_id: Protein ID string
        
    Returns:
        'disprot', 'ideal', or 'uniprot'
    """
    if protein_id.startswith('DP'):
        return 'disprot'
    elif protein_id.startswith('IID'):
        return 'ideal'
    else:
        return 'uniprot'


def parse_morf_fasta(fasta_file: str) -> pd.DataFrame:
    """
    Parse FASTA file and extract MoRF data.
    
    Format:
    >ID
    sequence
    MoRFTest annotation (skip)
    MoRFTrain annotation (use this - third line after ID)
    
    Args:
        fasta_file: Path to FASTA file
        
    Returns:
        DataFrame with columns: ID, sequence, morf_residues, database_type
    """
    logger.info(f"Parsing {fasta_file}...")
    
    data = []
    
    with open(fasta_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Check for ID line
        if line.startswith('>'):
            protein_id = line[1:].strip()  # Remove '>'
            
            # Get sequence (next line)
            if i + 1 < len(lines):
                sequence = lines[i + 1].strip()
            else:
                logger.warning(f"Missing sequence for {protein_id}")
                i += 1
                continue
            
            # Skip MoRFTest annotation (second line after ID)
            # Get MoRFTrain annotation (third line after ID)
            if i + 3 < len(lines):
                morftrain_annotation = lines[i + 3].strip()
            else:
                logger.warning(f"Missing MoRFTrain annotation for {protein_id}")
                i += 1
                continue
            
            # Extract MoRF positions
            morf_residues = extract_morf_positions(morftrain_annotation, sequence)
            
            # Classify database type
            database_type = classify_database_type(protein_id)
            
            data.append({
                'ID': protein_id,
                'sequence': sequence,
                'morf_residues': morf_residues,
                'database_type': database_type
            })
            
            i += 4  # Move to next entry
        else:
            i += 1
    
    df = pd.DataFrame(data)
    logger.info(f"  Extracted {len(df)} entries from {fasta_file}")
    return df


def compile_morf_data(fasta_files: List[str]) -> pd.DataFrame:
    """
    Process CV1-CV4 files and combine into single DataFrame.
    
    Args:
        fasta_files: List of paths to FASTA files
        
    Returns:
        Combined DataFrame
    """
    logger.info("Compiling MoRF data from FASTA files...")
    
    all_dfs = []
    for fasta_file in fasta_files:
        if os.path.exists(fasta_file):
            df = parse_morf_fasta(fasta_file)
            all_dfs.append(df)
        else:
            logger.warning(f"File not found: {fasta_file}")
    
    if not all_dfs:
        logger.error("No FASTA files could be parsed")
        return pd.DataFrame(columns=['ID', 'sequence', 'morf_residues', 'database_type'])
    
    combined_df = pd.concat(all_dfs, ignore_index=True)
    logger.info(f"Total entries compiled: {len(combined_df)}")
    
    return combined_df


def map_id_to_uniprot(protein_id: str, from_db: str, cache: dict) -> Optional[str]:
    """
    Map Disprot/Ideal ID to UniProt using UniProt ID Mapping API.
    
    Args:
        protein_id: Disprot or Ideal ID
        from_db: 'DisProt' or 'IDEAL'
        cache: Dictionary to cache results
        
    Returns:
        UniProt ID or None if not found
    """
    # Check cache first
    cache_key = f"{from_db}:{protein_id}"
    if cache_key in cache:
        return cache[cache_key]
    
    # UniProt ID Mapping API
    url = "https://rest.uniprot.org/idmapping/run"
    
    # Map from DisProt/IDEAL to UniProtKB
    data = {
        'from': from_db,
        'to': 'UniProtKB',
        'ids': protein_id
    }
    
    try:
        response = requests.post(url, data=data, timeout=30)
        response.raise_for_status()
        
        job_id = response.json().get('jobId')
        if not job_id:
            logger.error(f"No job ID returned for {protein_id}")
            cache[cache_key] = None
            return None
        
        # Poll for results
        status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
        max_attempts = 30
        for attempt in range(max_attempts):
            time.sleep(API_DELAY)
            status_response = requests.get(status_url, timeout=30)
            status_response.raise_for_status()
            status_data = status_response.json()
            
            job_status = status_data.get('jobStatus')
            
            # Check if results are directly in the status response (API sometimes returns results here)
            if 'results' in status_data and len(status_data['results']) > 0:
                # Extract UniProt ID from results
                result = status_data['results'][0]
                if 'to' in result:
                    to_data = result['to']
                    if isinstance(to_data, dict) and 'primaryAccession' in to_data:
                        uniprot_id = to_data['primaryAccession']
                        cache[cache_key] = uniprot_id
                        return uniprot_id
                    elif isinstance(to_data, str):
                        # Sometimes 'to' is just a string (UniProt ID)
                        uniprot_id = to_data
                        cache[cache_key] = uniprot_id
                        return uniprot_id
            
            if job_status == 'FINISHED':
                # Get results from stream endpoint
                results_url = f"https://rest.uniprot.org/idmapping/stream/{job_id}?format=tsv"
                results_response = requests.get(results_url, timeout=30)
                results_response.raise_for_status()
                
                # Parse TSV results
                lines = results_response.text.strip().split('\n')
                if len(lines) > 1:  # Header + data
                    # First result line (skip header)
                    parts = lines[1].split('\t')
                    if len(parts) > 1:
                        uniprot_id = parts[1].strip()
                        if uniprot_id:
                            cache[cache_key] = uniprot_id
                            return uniprot_id
                
                logger.warning(f"No mapping found for {protein_id}")
                cache[cache_key] = None
                return None
            elif job_status == 'ERROR':
                error_msg = status_data.get('error', status_data.get('messages', ['Unknown error']))
                logger.error(f"Job error for {protein_id}: {error_msg}")
                cache[cache_key] = None
                return None
            elif job_status in ['RUNNING', 'NEW']:
                continue
            elif job_status is None:
                # Job status is None but no results yet - wait and retry
                if attempt < max_attempts - 1:
                    time.sleep(1)
                    continue
                else:
                    logger.warning(f"Job status is None for {protein_id} after {max_attempts} attempts")
                    cache[cache_key] = None
                    return None
            else:
                logger.warning(f"Unexpected job status for {protein_id}: {job_status}")
                cache[cache_key] = None
                return None
        
        logger.warning(f"Timeout waiting for mapping result for {protein_id}")
        cache[cache_key] = None
        return None
        
    except requests.exceptions.RequestException as e:
        logger.error(f"Error mapping {protein_id} from {from_db}: {e}")
        cache[cache_key] = None
        return None


def map_all_ids_to_uniprot(morf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Map all Disprot/Ideal IDs to UniProt and add uniprot_id column.
    
    Args:
        morf_df: DataFrame with MoRF data
        
    Returns:
        DataFrame with uniprot_id and mapping_error columns added
    """
    logger.info("Mapping Disprot/Ideal IDs to UniProt...")
    
    morf_df = morf_df.copy()
    morf_df['uniprot_id'] = None
    morf_df['mapping_error'] = None
    
    cache = {}
    
    # Filter rows that need mapping
    to_map = morf_df[morf_df['database_type'].isin(['disprot', 'ideal'])]
    logger.info(f"  Found {len(to_map)} IDs to map")
    
    for idx, row in to_map.iterrows():
        protein_id = row['ID']
        db_type = row['database_type']
        
        # Determine source database (UniProt API uses 'DisProt' and 'IDEAL')
        from_db = 'DisProt' if db_type == 'disprot' else 'IDEAL'
        
        logger.info(f"  Mapping {protein_id} ({from_db})...")
        uniprot_id = map_id_to_uniprot(protein_id, from_db, cache)
        time.sleep(API_DELAY)
        
        if uniprot_id:
            morf_df.at[idx, 'uniprot_id'] = uniprot_id
            logger.info(f"    -> {uniprot_id}")
        else:
            morf_df.at[idx, 'mapping_error'] = f"Failed to map {protein_id} from {from_db}"
            logger.warning(f"    -> Mapping failed")
    
    # For uniprot type, use original ID
    uniprot_mask = morf_df['database_type'] == 'uniprot'
    morf_df.loc[uniprot_mask, 'uniprot_id'] = morf_df.loc[uniprot_mask, 'ID']
    
    logger.info(f"  Mapping complete. {len(morf_df[morf_df['uniprot_id'].notna()])} IDs have UniProt mappings")
    
    return morf_df


def load_and_group_ortholog_data(input_file: str) -> pd.DataFrame:
    """
    Load CSV and group by unique (protein_id, human_ortholog_id) pairs.
    
    Args:
        input_file: Path to human_ortholog_mapped.csv
        
    Returns:
        DataFrame with one row per unique pair
    """
    logger.info(f"Loading ortholog data from {input_file}...")
    
    df = pd.read_csv(input_file)
    logger.info(f"  Loaded {len(df)} rows")
    
    # Group by unique (protein_id, human_ortholog_id) pairs
    # Keep: protein_id, sequence (from full_sequence), length, human_ortholog_id, human_ortholog_sequence
    grouped = df.groupby(['protein_id', 'human_ortholog_id']).first().reset_index()
    
    # Extract required columns
    result_df = pd.DataFrame({
        'protein_id': grouped['protein_id'],
        'sequence': grouped.get('full_sequence', grouped.get('sequence', '')),
        'length': grouped.get('length', None),
        'human_ortholog_id': grouped['human_ortholog_id'],
        'human_ortholog_sequence': grouped.get('human_ortholog_sequence', '')
    })
    
    logger.info(f"  Grouped to {len(result_df)} unique (protein_id, human_ortholog_id) pairs")
    
    return result_df


def match_morf_to_orthologs(morf_df: pd.DataFrame, ortholog_df: pd.DataFrame) -> pd.DataFrame:
    """
    Match MoRF uniprot_ids to protein_id or human_ortholog_id columns.
    
    Args:
        morf_df: DataFrame with MoRF data (must have uniprot_id column)
        ortholog_df: DataFrame with ortholog data
        
    Returns:
        DataFrame with matched data
    """
    logger.info("Matching MoRF data to ortholog data...")
    
    result_rows = []
    
    # Create index for faster lookup
    ortholog_by_protein = ortholog_df.set_index('protein_id')
    ortholog_by_human = ortholog_df.set_index('human_ortholog_id')
    
    for idx, morf_row in morf_df.iterrows():
        uniprot_id = morf_row.get('uniprot_id')
        
        if pd.isna(uniprot_id) or not uniprot_id:
            continue
        
        morf_residues = morf_row.get('morf_residues', '[]')
        
        # Check if matches protein_id
        if uniprot_id in ortholog_by_protein.index:
            ortholog_row = ortholog_by_protein.loc[uniprot_id]
            
            # Handle both Series (single row) and DataFrame (multiple rows) cases
            if isinstance(ortholog_row, pd.DataFrame):
                ortholog_row = ortholog_row.iloc[0]
            
            # Convert to dict to access all values, including index
            row_dict = ortholog_row.to_dict()
            # Add the index value (protein_id) to the dict
            row_dict['protein_id'] = uniprot_id
            
            result_rows.append({
                'protein_id': row_dict['protein_id'],
                'sequence': row_dict['sequence'],
                'length': row_dict.get('length', None),
                'morf_residues': morf_residues,
                'human_ortholog_id': row_dict['human_ortholog_id'],
                'human_ortholog_sequence': row_dict['human_ortholog_sequence'],
                'human_morf_residues': None
            })
            logger.info(f"  Matched {uniprot_id} to protein_id")
        
        # Check if matches human_ortholog_id
        if uniprot_id in ortholog_by_human.index:
            ortholog_row = ortholog_by_human.loc[uniprot_id]
            
            # Handle both Series (single row) and DataFrame (multiple rows) cases
            if isinstance(ortholog_row, pd.DataFrame):
                ortholog_row = ortholog_row.iloc[0]
            
            # Convert to dict to access all values, including index
            row_dict = ortholog_row.to_dict()
            # Add the index value (human_ortholog_id) to the dict
            row_dict['human_ortholog_id'] = uniprot_id
            
            result_rows.append({
                'protein_id': row_dict['protein_id'],
                'sequence': row_dict['sequence'],
                'length': row_dict.get('length', None),
                'morf_residues': None,
                'human_ortholog_id': row_dict['human_ortholog_id'],
                'human_ortholog_sequence': row_dict['human_ortholog_sequence'],
                'human_morf_residues': morf_residues
            })
            logger.info(f"  Matched {uniprot_id} to human_ortholog_id")
    
    result_df = pd.DataFrame(result_rows)
    logger.info(f"  Total matches: {len(result_df)}")
    
    return result_df


def main():
    """Main workflow function."""
    parser = argparse.ArgumentParser(
        description='Match MoRF residues to proteins in human_ortholog_mapped.csv'
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Path to input human_ortholog_mapped.csv file'
    )
    parser.add_argument(
        '--species',
        type=str,
        required=True,
        choices=['mouse', 'rat'],
        help='Species: mouse or rat'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Path to output matched_morfs.csv file (default: same directory as input)'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        return
    
    # Determine output file path
    if args.output:
        output_path = args.output
    else:
        input_dir = os.path.dirname(args.input)
        output_path = os.path.join(input_dir, 'matched_morfs.csv')
    
    # Determine FASTA file paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    morf_dir = os.path.join(project_root, 'data', 'MoRF2')
    
    fasta_files = [
        os.path.join(morf_dir, 'CV1.af'),
        os.path.join(morf_dir, 'CV2.af'),
        os.path.join(morf_dir, 'CV3.af'),
        os.path.join(morf_dir, 'CV4.af')
    ]
    
    logger.info("=" * 80)
    logger.info("MoRF Interface Matching")
    logger.info("=" * 80)
    logger.info(f"Input file: {args.input}")
    logger.info(f"Species: {args.species}")
    logger.info(f"Output file: {output_path}")
    logger.info("")
    
    # Step 1: Compile MoRF data from FASTA files
    morf_df = compile_morf_data(fasta_files)
    
    if len(morf_df) == 0:
        logger.error("No MoRF data extracted. Exiting.")
        return
    
    logger.info("")
    
    # Step 2: Map Disprot/Ideal IDs to UniProt
    morf_df = map_all_ids_to_uniprot(morf_df)
    
    logger.info("")
    
    # Step 3: Load and group ortholog data
    ortholog_df = load_and_group_ortholog_data(args.input)
    
    logger.info("")
    
    # Step 4: Match MoRF data to ortholog data
    matched_df = match_morf_to_orthologs(morf_df, ortholog_df)
    
    logger.info("")
    
    # Step 5: Save results
    logger.info(f"Saving results to {output_path}...")
    matched_df.to_csv(output_path, index=False)
    logger.info(f"✓ Successfully saved {len(matched_df)} matched rows")
    
    # Print summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    logger.info(f"Total MoRF entries: {len(morf_df)}")
    logger.info(f"  - UniProt IDs: {len(morf_df[morf_df['database_type'] == 'uniprot'])}")
    logger.info(f"  - Disprot IDs: {len(morf_df[morf_df['database_type'] == 'disprot'])}")
    logger.info(f"  - Ideal IDs: {len(morf_df[morf_df['database_type'] == 'ideal'])}")
    logger.info(f"Successfully mapped IDs: {len(morf_df[morf_df['uniprot_id'].notna()])}")
    logger.info(f"Total matches: {len(matched_df)}")
    logger.info(f"  - Matches to protein_id: {len(matched_df[matched_df['morf_residues'].notna()])}")
    logger.info(f"  - Matches to human_ortholog_id: {len(matched_df[matched_df['human_morf_residues'].notna()])}")
    logger.info(f"Output file: {output_path}")
    logger.info("")


if __name__ == '__main__':
    main()


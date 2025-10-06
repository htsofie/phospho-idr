#!/usr/bin/env python3
"""
Protein name-based mapping using UniProt search.

This script searches UniProt by protein description to find the correct protein
and extract its full sequence when other mapping methods fail.

Usage:
    python name_search.py --input data.csv --species rat --output results.csv
"""

# Run this script on the command line with: python scripts/name_search.py -i data/processed/rat/cleaned_test_data_grouped_processed_aligned_remapped.csv -s rat

import argparse
import pandas as pd
import requests
import time
import re
from typing import Optional, Tuple, List
import sys
import os

# Add scripts directory to path
sys.path.append(os.path.join(os.path.dirname(__file__)))
from get_full_seq import get_uniprot_sequence

# Set up logging
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# UniProt search API
UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def search_uniprot_by_name(protein_name: str, species: str, max_results: int = 5) -> List[dict]:
    """
    Search UniProt by protein name/description for a specific species.
    
    Args:
        protein_name: Protein name/description to search for
        species: Species name (rat, mouse, human)
        max_results: Maximum number of results to return
        
    Returns:
        List of dictionaries with protein information
    """
    logger.info(f"Searching UniProt for: '{protein_name}' in {species}")
    
    # Map species names to UniProt organism codes
    species_map = {
        'rat': 'Rattus norvegicus',
        'mouse': 'Mus musculus', 
        'human': 'Homo sapiens'
    }
    
    if species.lower() not in species_map:
        logger.error(f"Unknown species: {species}")
        return []
    
    organism = species_map[species.lower()]
    
    # Clean up protein name for search
    # Remove common prefixes/suffixes that might interfere with search
    clean_name = protein_name.strip()
    
    # Build search query - use correct field names
    query = f'"{clean_name}" AND organism_name:"{organism}"'
    
    params = {
        'query': query,
        'format': 'json',
        'size': max_results
    }
    
    try:
        response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
        response.raise_for_status()
        
        data = response.json()
        results = []
        
        for entry in data.get('results', []):
            # Determine reviewed status using available fields
            entry_type = str(entry.get('entryType', '')).lower()
            reviewed_flag = bool(entry.get('reviewed', False))
            is_reviewed = reviewed_flag or (entry_type == 'reviewed')

            protein_info = {
                'accession': entry.get('primaryAccession', ''),
                'id': entry.get('uniProtkbId', ''),
                'name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                'organism': entry.get('organism', {}).get('scientificName', ''),
                'sequence': entry.get('sequence', {}).get('value', ''),
                'length': entry.get('sequence', {}).get('length', 0),
                'reviewed': is_reviewed,
            }
            results.append(protein_info)
        
        logger.info(f"Found {len(results)} UniProt entries for '{protein_name}'")
        return results
        
    except requests.exceptions.RequestException as e:
        logger.error(f"UniProt search failed: {e}")
        return []
    except Exception as e:
        logger.error(f"Error processing UniProt search results: {e}")
        return []


def extract_first_description(protein_description: str) -> str:
    """
    Extract the first protein description from a semicolon-separated list.
    
    Args:
        protein_description: Semicolon-separated protein descriptions
        
    Returns:
        First description, cleaned up
    """
    if pd.isna(protein_description) or not protein_description:
        return ""
    
    # Split by semicolon and take the first one
    first_desc = str(protein_description).split(';')[0].strip()
    
    # Clean up common issues
    first_desc = re.sub(r'\s+', ' ', first_desc)  # Remove extra whitespace
    first_desc = first_desc.strip()
    
    return first_desc


def get_all_protein_matches(protein_name: str, species: str) -> List[dict]:
    """
    Get all matching proteins from UniProt search. Keep UniProt's original
    relevance order within equal-length accessions, but order groups by
    accession length (shortest first). Only entries with sequences kept.
    """
    results = search_uniprot_by_name(protein_name, species, max_results=10)
    if not results:
        return []

    # Keep only entries that include a full sequence
    results_with_seq = [p for p in results if p.get('sequence')]

    # Stable sort by accession length; Python sort is stable so UniProt
    # original order is preserved within each length bucket
    ordered = sorted(results_with_seq, key=lambda p: len(p.get('accession', '')))

    logger.info(f"Found {len(ordered)} UniProt entries for '{protein_name}' (shortest accession first)")
    top_preview = ", ".join(p['accession'] for p in ordered[:3])
    if top_preview:
        logger.info(f"  Top accessions: {top_preview}")

    return ordered


def process_dataset_by_name(df: pd.DataFrame, species: str) -> pd.DataFrame:
    """
    Process dataset by searching protein names in UniProt.
    
    Args:
        df: Input DataFrame
        species: Species name
        
    Returns:
        Updated DataFrame with name-based mappings
    """
    logger.info(f"Processing dataset by protein name search for {species}")
    
    # Find rows that need name-based mapping
    # Look for rows with protein_description but no successful mapping
    needs_mapping = df[
        (df.get('protein_description', '').notna()) & 
        (df.get('protein_description', '') != '') &
        (df.get('uniprot_mapped_id', '').isna() | (df.get('uniprot_mapped_id', '') == ''))
    ]
    
    logger.info(f"Found {len(needs_mapping)} rows that need name-based mapping")
    
    if len(needs_mapping) == 0:
        logger.info("No rows need name-based mapping")
        return df
    
    successful = 0
    failed = 0
    
    # Process each unique protein group
    protein_groups = needs_mapping['Protein'].unique()
    
    for protein_id in protein_groups:
        logger.info(f"Processing protein group: {protein_id}")
        
        group_rows = df[df['Protein'] == protein_id]
        first_row = group_rows.iloc[0]
        
        # Extract first protein description
        protein_description = first_row.get('protein_description', '')
        if not protein_description:
            logger.warning(f"No protein description for group {protein_id}")
            failed += 1
            continue
        
        first_desc = extract_first_description(protein_description)
        if not first_desc:
            logger.warning(f"Could not extract description from: {protein_description}")
            failed += 1
            continue
        
        logger.info(f"Searching for: '{first_desc}'")
        
        # Search UniProt for all matches
        all_matches = get_all_protein_matches(first_desc, species)
        
        if all_matches:
            # Update all rows in this protein group with multiple candidate sequences
            group_idx = group_rows.index
            
            # Store all candidate UniProt IDs and sequences
            candidate_ids = []
            candidate_sequences = []
            candidate_lengths = []
            
            for match in all_matches:
                if match['sequence']:
                    candidate_ids.append(match['accession'])
                    candidate_sequences.append(match['sequence'])
                    candidate_lengths.append(match['length'])
            
            if candidate_ids:
                # Store as semicolon-separated lists for multiple candidates
                df.loc[group_idx, 'uniprot_mapped_id'] = ';'.join(candidate_ids)
                df.loc[group_idx, 'mapping_source'] = 'NAME_SEARCH_MULTIPLE'
                df.loc[group_idx, 'full_sequence'] = ';'.join(candidate_sequences)
                df.loc[group_idx, 'sequence_length'] = ';'.join(map(str, candidate_lengths))
                df.loc[group_idx, 'sequence_type'] = 'NAME_MAPPED_MULTIPLE'
                df.loc[group_idx, 'manual_review_flag'] = False  # Clear manual review flag
                df.loc[group_idx, 'alignment_error'] = ''  # Clear alignment errors
                
                successful += 1
                logger.info(f"Found {len(candidate_ids)} candidates for group {protein_id}: {candidate_ids[0]}{'...' if len(candidate_ids) > 1 else ''}")
            else:
                failed += 1
                logger.warning(f"No valid sequences found for: '{first_desc}'")
        else:
            failed += 1
            logger.warning(f"Failed to find UniProt match for: '{first_desc}'")
        
        # Add delay to avoid overwhelming UniProt API
        time.sleep(1)
    
    logger.info(f"Name-based mapping complete: {successful} successful, {failed} failed")
    return df


def main():
    parser = argparse.ArgumentParser(description='Protein name-based mapping using UniProt search')
    parser.add_argument('--input', '-i', required=True, help='Input CSV/Parquet file')
    parser.add_argument('--species', '-s', required=True, choices=['mouse', 'rat', 'human'], help='Species name')
    parser.add_argument('--output', '-o', help='Output CSV file (optional - if not provided, creates name_search_results.csv in same directory)')
    
    args = parser.parse_args()
    
    # Load data
    logger.info(f"Loading data from {args.input}")
    if args.input.endswith('.parquet'):
        df = pd.read_parquet(args.input)
    else:
        df = pd.read_csv(args.input)
    
    logger.info(f"Loaded {len(df)} rows")
    
    # Process dataset
    updated_df = process_dataset_by_name(df, args.species)
    
    # Save results
    if args.output:
        output_file = args.output
    else:
        # Create name_search_results.csv in the same directory as input file
        input_path = os.path.dirname(args.input)
        output_file = os.path.join(input_path, "name_search_results.csv")
    
    logger.info(f"Saving updated dataset to {output_file}")
    updated_df.to_csv(output_file, index=False)
    
    # Print summary
    name_mapped = updated_df[updated_df['mapping_source'] == 'NAME_SEARCH']['Protein'].nunique()
    total_groups = updated_df['Protein'].nunique()
    
    logger.info(f"Processing complete:")
    logger.info(f"  Total protein groups: {total_groups}")
    logger.info(f"  Groups mapped by name search: {name_mapped}")
    logger.info(f"  Next step: Run align_to_full_seq.py on the output file to perform alignment")


if __name__ == "__main__":
    main()

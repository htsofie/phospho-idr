#!/usr/bin/env python3
"""
Consolidate duplicate mouse and human interactions in interface data.

This script processes int_interfaces.csv files to merge mouse and human interactions
that represent the same interaction, keeping mouse data in specific columns and
human data in all other columns.

Usage:
    python scripts/consolidate_duplicate_interactions.py --input data/processed/mouse/int_interfaces.csv --species mouse --output data/processed/mouse/int_interfaces_consolidated.csv
    python scripts/consolidate_duplicate_interactions.py --input data/processed/rat/int_interfaces.csv --species rat
"""

import pandas as pd
import argparse
import logging
import os
import sys
import time
from typing import Optional

# Add scripts directory to path for imports
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

# Import functions from map_human_orthologs
from archive.map_human_orthologs import (
    get_uniprot_entry,
    extract_gene_id,
    get_human_ortholog_from_alliance,
    get_human_uniprot_id
)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# API rate limiting delay (seconds)
API_DELAY = 0.15


def get_human_ortholog_for_mouse_protein(mouse_protein_id: str, species: str, cache: dict) -> Optional[str]:
    """
    Get human ortholog UniProt ID for a mouse/rat protein ID.
    
    Uses the same pipeline as map_human_orthologs.py.
    
    Args:
        mouse_protein_id: Mouse/rat UniProt protein ID
        species: 'mouse' or 'rat'
        cache: Dictionary to cache results (key: protein_id, value: human_ortholog_id)
        
    Returns:
        Human ortholog UniProt ID, or None if not found
    """
    # Check cache first
    if mouse_protein_id in cache:
        return cache[mouse_protein_id]
    
    try:
        # Step 1: Get UniProt entry and extract gene ID
        uniprot_data = get_uniprot_entry(mouse_protein_id)
        time.sleep(API_DELAY)
        
        if not uniprot_data:
            cache[mouse_protein_id] = None
            return None
        
        gene_id = extract_gene_id(uniprot_data, species)
        if not gene_id:
            cache[mouse_protein_id] = None
            return None
        
        time.sleep(API_DELAY)
        
        # Step 2: Query Alliance Genome for human ortholog
        hgnc_id = get_human_ortholog_from_alliance(gene_id)
        time.sleep(API_DELAY)
        
        if not hgnc_id:
            cache[mouse_protein_id] = None
            return None
        
        # Step 3: Search UniProt for human protein
        human_uniprot_id = get_human_uniprot_id(hgnc_id)
        time.sleep(API_DELAY)
        
        cache[mouse_protein_id] = human_uniprot_id
        return human_uniprot_id
        
    except Exception as e:
        logger.error(f"Error getting human ortholog for {mouse_protein_id}: {e}")
        cache[mouse_protein_id] = None
        return None


def consolidate_duplicate_interactions(interfaces_df, species):
    """
    Consolidate mouse and human interactions that represent the same interaction.
    
    For rows with mouse_interactor_id:
    1. Find human ortholog of mouse_interactor_id
    2. Check if it matches any human_interactor_id for the same protein_id
    3. If match found, merge the rows
    
    Args:
        interfaces_df: DataFrame with mouse and human interactions
        species: 'mouse' or 'rat'
        
    Returns:
        Consolidated DataFrame
    """
    logger.info("Consolidating duplicate mouse/human interactions...")
    
    # Cache for human ortholog lookups
    ortholog_cache = {}
    
    # Group by protein_id
    consolidated_rows = []
    processed_indices = set()
    
    for protein_id, group in interfaces_df.groupby('protein_id'):
        group_rows = group.reset_index(drop=True)
        group_indices = group.index.tolist()
        
        # Find human orthologs for mouse interactors in this group
        mouse_rows_with_orthologs = []
        for idx, row in group_rows.iterrows():
            if pd.notna(row.get('mouse_interactor_id')):
                mouse_interactor = row['mouse_interactor_id']
                human_ortholog = get_human_ortholog_for_mouse_protein(
                    mouse_interactor, species, ortholog_cache
                )
                mouse_rows_with_orthologs.append((idx, row, human_ortholog))
        
        # Check for matches and merge
        for mouse_idx, mouse_row, mouse_interactor_human_ortholog in mouse_rows_with_orthologs:
            if mouse_interactor_human_ortholog is None:
                # No human ortholog found, keep mouse row as-is
                if group_indices[mouse_idx] not in processed_indices:
                    consolidated_rows.append(mouse_row)
                    processed_indices.add(group_indices[mouse_idx])
                continue
            
            # Look for matching human interaction in this group
            match_found = False
            for human_idx, human_row in group_rows.iterrows():
                if group_indices[human_idx] in processed_indices:
                    continue
                    
                if pd.notna(human_row.get('human_interactor_id')):
                    if human_row['human_interactor_id'] == mouse_interactor_human_ortholog:
                        # Same interaction! Merge the rows
                        # Start with human row, then copy specific columns from mouse row
                        merged_row = human_row.copy()
                        
                        # Copy these columns from mouse row
                        mouse_columns_to_copy = [
                            'protein_id',
                            'full_sequence',
                            'uniprot_protein_id_length',
                            'mouse_interactor_id',
                            'phosphoprotein_IRES',
                            'interactor_IRES'
                        ]
                        
                        for col in mouse_columns_to_copy:
                            if col in mouse_row and col in merged_row:
                                merged_row[col] = mouse_row[col]
                        
                        consolidated_rows.append(merged_row)
                        processed_indices.add(group_indices[human_idx])
                        processed_indices.add(group_indices[mouse_idx])
                        match_found = True
                        break
            
            # No match found, keep mouse row
            if not match_found and group_indices[mouse_idx] not in processed_indices:
                consolidated_rows.append(mouse_row)
                processed_indices.add(group_indices[mouse_idx])
        
        # Add any remaining human-only rows
        for idx, row in group_rows.iterrows():
            if group_indices[idx] not in processed_indices:
                consolidated_rows.append(row)
                processed_indices.add(group_indices[idx])
    
    result_df = pd.DataFrame(consolidated_rows)
    logger.info(f"Consolidation complete: {len(result_df)} rows (from {len(interfaces_df)} input rows)")
    
    return result_df


def main():
    """Main workflow function."""
    parser = argparse.ArgumentParser(
        description='Consolidate duplicate mouse and human interactions in interface data'
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Path to input int_interfaces.csv file'
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
        help='Path to output CSV file (default: input file with _consolidated suffix)'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        return
    
    # Determine output path
    if args.output:
        output_path = args.output
    else:
        # Add _consolidated before .csv extension
        base, ext = os.path.splitext(args.input)
        output_path = f"{base}_consolidated{ext}"
    
    logger.info("=" * 80)
    logger.info("Consolidate Duplicate Interactions")
    logger.info("=" * 80)
    logger.info(f"Input file: {args.input}")
    logger.info(f"Species: {args.species}")
    logger.info(f"Output file: {output_path}")
    logger.info("")
    
    # Load input CSV
    logger.info("Loading input file...")
    try:
        df = pd.read_csv(args.input)
        logger.info(f"Loaded {len(df)} rows")
    except Exception as e:
        logger.error(f"Error reading CSV file: {e}")
        return
    
    # Check required columns
    required_cols = ['protein_id']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.error(f"Error: Missing required columns: {missing_cols}")
        logger.error(f"Available columns: {list(df.columns)}")
        return
    
    # Only consolidate for mouse species
    if args.species != 'mouse':
        logger.warning("Consolidation is only applicable for mouse species. Rat data will be passed through unchanged.")
        result_df = df
    else:
        # Consolidate duplicate interactions
        result_df = consolidate_duplicate_interactions(df, args.species)
        logger.info("")
    
    # Remove any duplicate rows that might have been created
    logger.info("Removing duplicate rows...")
    initial_count = len(result_df)
    result_df = result_df.drop_duplicates()
    final_count = len(result_df)
    if initial_count != final_count:
        logger.info(f"Removed {initial_count - final_count} duplicate rows")
    logger.info("")
    
    # Save to output file
    logger.info(f"Saving results to {output_path}...")
    try:
        result_df.to_csv(output_path, index=False)
        logger.info(f"✓ Successfully saved {len(result_df)} rows")
    except Exception as e:
        logger.error(f"Error saving CSV file: {e}")
        return
    
    # Print summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    logger.info(f"Input rows: {len(df)}")
    logger.info(f"Output rows: {len(result_df)}")
    logger.info(f"Rows removed: {len(df) - len(result_df)}")
    logger.info(f"Output file: {output_path}")
    logger.info("")


if __name__ == '__main__':
    main()


#!/usr/bin/env python3
"""
Map interfaces from Interactome Insider to protein IDs and human orthologs.

This script processes human_ortholog_mapped.csv files to identify protein-protein
interactions from Interactome Insider database, mapping both source species proteins
(mouse only) and human orthologs to their interaction partners and interface residues.

Usage:
    python scripts/map_interfaces_int_insider.py --input data/processed/mouse/human_ortholog_mapped.csv --species mouse
    python scripts/map_interfaces_int_insider.py --input data/processed/rat/human_ortholog_mapped.csv --species rat
"""

import pandas as pd
import argparse
import logging
import os
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def is_empty_ires(ires_value):
    """
    Check if IRES value is empty.
    
    Args:
        ires_value: IRES value from interface file
        
    Returns:
        True if empty, False otherwise
    """
    if pd.isna(ires_value):
        return True
    if isinstance(ires_value, str):
        ires_value = ires_value.strip()
        if ires_value == '' or ires_value == '[]':
            return True
    return False


def filter_valid_interactions(interface_df):
    """
    Filter interface dataframe to keep only valid interactions.
    
    Filters out:
    - Interactions where Source == "ECLAIR"
    - Interactions where P1_IRES or P2_IRES is empty/null/"[]"
    
    Args:
        interface_df: DataFrame with interface data
        
    Returns:
        Filtered DataFrame
    """
    # Filter out ECLAIR sources
    filtered = interface_df[interface_df['Source'] != 'ECLAIR'].copy()
    
    # Filter out rows where P1_IRES or P2_IRES is empty
    mask = filtered.apply(
        lambda row: not is_empty_ires(row['P1_IRES']) and not is_empty_ires(row['P2_IRES']),
        axis=1
    )
    filtered = filtered[mask].copy()
    
    return filtered


def process_mouse_interfaces(interfaces_df, mouse_interface_file):
    """
    Process mouse interface matching.
    
    For each protein_id, find matches in M_musculus_interfacesHQ.csv and create
    interaction rows with mouse_interactor_id, phosphoprotein_IRES, and interactor_IRES.
    
    Args:
        interfaces_df: Base dataframe with protein information
        mouse_interface_file: Path to M_musculus_interfacesHQ.csv
        
    Returns:
        DataFrame with mouse interaction data added
    """
    logger.info("Loading mouse interface file...")
    mouse_interfaces = pd.read_csv(mouse_interface_file)
    logger.info(f"Loaded {len(mouse_interfaces)} mouse interface records")
    
    # Filter valid interactions
    mouse_interfaces = filter_valid_interactions(mouse_interfaces)
    logger.info(f"After filtering: {len(mouse_interfaces)} valid mouse interface records")
    
    # Initialize new columns
    interfaces_df['mouse_interactor_id'] = None
    interfaces_df['phosphoprotein_IRES'] = None
    interfaces_df['interactor_IRES'] = None
    # Initialize human columns to None (will be cleared for mouse-only rows)
    interfaces_df['human_interactor_id'] = None
    interfaces_df['human_phosphoprotein_IRES'] = None
    interfaces_df['human_interactor_IRES'] = None
    
    # List to store new rows (for multiple matches)
    new_rows = []
    # Track seen interactions to avoid duplicates
    seen_interactions = set()
    
    logger.info("Matching mouse protein IDs to interfaces...")
    for idx, row in interfaces_df.iterrows():
        protein_id = row['protein_id']
        
        # Find matches where protein_id is in P1 or P2
        p1_matches = mouse_interfaces[mouse_interfaces['P1'] == protein_id]
        p2_matches = mouse_interfaces[mouse_interfaces['P2'] == protein_id]
        
        matches_found = False
        
        # Process P1 matches
        for _, match in p1_matches.iterrows():
            # Create unique key for this interaction to avoid duplicates
            interaction_key = (protein_id, match['P2'], match['P1_IRES'], match['P2_IRES'])
            if interaction_key in seen_interactions:
                continue
            seen_interactions.add(interaction_key)
            
            matches_found = True
            new_row = row.copy()
            new_row['mouse_interactor_id'] = match['P2']
            new_row['phosphoprotein_IRES'] = match['P1_IRES']
            new_row['interactor_IRES'] = match['P2_IRES']
            # Clear human interaction columns for this mouse-only row
            new_row['human_interactor_id'] = None
            new_row['human_phosphoprotein_IRES'] = None
            new_row['human_interactor_IRES'] = None
            new_rows.append(new_row)
        
        # Process P2 matches
        for _, match in p2_matches.iterrows():
            # Create unique key for this interaction to avoid duplicates
            interaction_key = (protein_id, match['P1'], match['P2_IRES'], match['P1_IRES'])
            if interaction_key in seen_interactions:
                continue
            seen_interactions.add(interaction_key)
            
            matches_found = True
            new_row = row.copy()
            new_row['mouse_interactor_id'] = match['P1']
            new_row['phosphoprotein_IRES'] = match['P2_IRES']
            new_row['interactor_IRES'] = match['P1_IRES']
            # Clear human interaction columns for this mouse-only row
            new_row['human_interactor_id'] = None
            new_row['human_phosphoprotein_IRES'] = None
            new_row['human_interactor_IRES'] = None
            new_rows.append(new_row)
        
        # If no matches, keep original row with None values
        if not matches_found:
            new_rows.append(row)
    
    # Create new dataframe from all rows
    result_df = pd.DataFrame(new_rows)
    logger.info(f"Mouse interface matching complete: {len(result_df)} rows (from {len(interfaces_df)} unique proteins)")
    
    return result_df


def process_human_interfaces(interfaces_df, human_interface_file):
    """
    Process human interface matching.
    
    For each row (including those with mouse matches), find matches where
    human_ortholog_id is in H_sapiens_interfacesHQ.csv and create new rows
    with human_interactor_id, human_phosphoprotein_IRES, and human_interactor_IRES.
    
    Args:
        interfaces_df: DataFrame with protein information (may include mouse matches)
        human_interface_file: Path to H_sapiens_interfacesHQ.csv
        
    Returns:
        DataFrame with human interaction data added
    """
    logger.info("Loading human interface file...")
    human_interfaces = pd.read_csv(human_interface_file)
    logger.info(f"Loaded {len(human_interfaces)} human interface records")
    
    # Filter valid interactions
    human_interfaces = filter_valid_interactions(human_interfaces)
    logger.info(f"After filtering: {len(human_interfaces)} valid human interface records")
    
    # Initialize new columns if they don't exist
    if 'human_interactor_id' not in interfaces_df.columns:
        interfaces_df['human_interactor_id'] = None
        interfaces_df['human_phosphoprotein_IRES'] = None
        interfaces_df['human_interactor_IRES'] = None
    
    # List to store new rows (for multiple matches)
    new_rows = []
    # Track seen interactions to avoid duplicates
    seen_interactions = set()
    
    logger.info("Matching human ortholog IDs to interfaces...")
    for idx, row in interfaces_df.iterrows():
        human_ortholog_id = row['human_ortholog_id']
        
        # Skip if no human ortholog ID
        if pd.isna(human_ortholog_id):
            new_rows.append(row)
            continue
        
        # Find matches where human_ortholog_id is in P1 or P2
        p1_matches = human_interfaces[human_interfaces['P1'] == human_ortholog_id]
        p2_matches = human_interfaces[human_interfaces['P2'] == human_ortholog_id]
        
        matches_found = False
        
        # Process P1 matches
        for _, match in p1_matches.iterrows():
            # Create unique key for this interaction to avoid duplicates
            interaction_key = (human_ortholog_id, match['P2'], match['P1_IRES'], match['P2_IRES'])
            if interaction_key in seen_interactions:
                continue
            seen_interactions.add(interaction_key)
            
            matches_found = True
            new_row = row.copy()
            new_row['human_interactor_id'] = match['P2']
            new_row['human_phosphoprotein_IRES'] = match['P1_IRES']
            new_row['human_interactor_IRES'] = match['P2_IRES']
            # Clear mouse interaction columns for this human-only row (if they exist)
            if 'mouse_interactor_id' in new_row:
                new_row['mouse_interactor_id'] = None
                new_row['phosphoprotein_IRES'] = None
                new_row['interactor_IRES'] = None
            new_rows.append(new_row)
        
        # Process P2 matches
        for _, match in p2_matches.iterrows():
            # Create unique key for this interaction to avoid duplicates
            interaction_key = (human_ortholog_id, match['P1'], match['P2_IRES'], match['P1_IRES'])
            if interaction_key in seen_interactions:
                continue
            seen_interactions.add(interaction_key)
            
            matches_found = True
            new_row = row.copy()
            new_row['human_interactor_id'] = match['P1']
            new_row['human_phosphoprotein_IRES'] = match['P2_IRES']
            new_row['human_interactor_IRES'] = match['P1_IRES']
            # Clear mouse interaction columns for this human-only row (if they exist)
            if 'mouse_interactor_id' in new_row:
                new_row['mouse_interactor_id'] = None
                new_row['phosphoprotein_IRES'] = None
                new_row['interactor_IRES'] = None
            new_rows.append(new_row)
        
        # Always keep the original row if it has mouse interactions (even if human matches were found)
        # This preserves mouse interaction data when both mouse and human matches exist
        if not matches_found:
            # No human matches, keep original row as-is
            new_rows.append(row)
        else:
            # Human matches found, but we still want to keep the original row if it has mouse interactions
            # Check if original row has mouse interactions
            if 'mouse_interactor_id' in row and pd.notna(row.get('mouse_interactor_id')):
                # Keep original row with mouse interactions (human columns remain None)
                new_rows.append(row)
    
    # Create new dataframe from all rows
    result_df = pd.DataFrame(new_rows)
    logger.info(f"Human interface matching complete: {len(result_df)} rows (from {len(interfaces_df)} input rows)")
    
    return result_df


def main():
    """Main workflow function."""
    parser = argparse.ArgumentParser(
        description='Map interfaces from Interactome Insider to protein IDs and human orthologs'
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
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        return
    
    # Determine output path (same directory as input)
    input_dir = os.path.dirname(args.input)
    output_path = os.path.join(input_dir, 'int_interfaces.csv')
    
    # Determine interface file paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    mouse_interface_file = os.path.join(project_root, 'data', 'interactome_insider', 'M_musculus_interfacesHQ.csv')
    human_interface_file = os.path.join(project_root, 'data', 'interactome_insider', 'H_sapiens_interfacesHQ.csv')
    
    logger.info("=" * 80)
    logger.info("Interactome Insider Interface Mapping")
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
    required_cols = ['protein_id', 'full_sequence', 'human_ortholog_id', 'human_ortholog_sequence']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.error(f"Error: Missing required columns: {missing_cols}")
        logger.error(f"Available columns: {list(df.columns)}")
        return
    
    # Create interfaces_df with one row per unique protein_id
    logger.info("Creating base interfaces dataframe...")
    interfaces_df = df[required_cols].drop_duplicates(subset=['protein_id']).copy()
    logger.info(f"Created base dataframe with {len(interfaces_df)} unique proteins")
    logger.info("")
    
    # Process mouse interfaces (mouse only)
    if args.species == 'mouse':
        if not os.path.exists(mouse_interface_file):
            logger.error(f"Mouse interface file not found: {mouse_interface_file}")
            return
        
        interfaces_df = process_mouse_interfaces(interfaces_df, mouse_interface_file)
        logger.info("")
    
    # Process human interfaces (both mouse and rat)
    if not os.path.exists(human_interface_file):
        logger.error(f"Human interface file not found: {human_interface_file}")
        return
    
    interfaces_df = process_human_interfaces(interfaces_df, human_interface_file)
    logger.info("")
    
    # Remove duplicate rows
    logger.info("Removing duplicate rows...")
    initial_count = len(interfaces_df)
    interfaces_df = interfaces_df.drop_duplicates()
    final_count = len(interfaces_df)
    if initial_count != final_count:
        logger.info(f"Removed {initial_count - final_count} duplicate rows")
    logger.info("")
    
    # Reorder columns to match specified order
    base_cols = ['protein_id', 'full_sequence']
    mouse_cols = ['mouse_interactor_id', 'phosphoprotein_IRES', 'interactor_IRES'] if args.species == 'mouse' else []
    human_cols = ['human_ortholog_id', 'human_ortholog_sequence', 'human_interactor_id', 'human_phosphoprotein_IRES', 'human_interactor_IRES']
    
    # Build column order
    col_order = base_cols + mouse_cols + human_cols
    
    # Ensure all columns exist (add missing ones as None)
    for col in col_order:
        if col not in interfaces_df.columns:
            interfaces_df[col] = None
    
    # Reorder dataframe
    interfaces_df = interfaces_df[col_order]
    
    # Save to output file
    logger.info(f"Saving results to {output_path}...")
    try:
        interfaces_df.to_csv(output_path, index=False)
        logger.info(f"✓ Successfully saved {len(interfaces_df)} rows")
    except Exception as e:
        logger.error(f"Error saving CSV file: {e}")
        return
    
    # Print summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    logger.info(f"Total unique proteins processed: {len(df['protein_id'].unique())}")
    logger.info(f"Total interaction rows created: {len(interfaces_df)}")
    
    if args.species == 'mouse':
        mouse_matches = interfaces_df['mouse_interactor_id'].notna().sum()
        logger.info(f"Mouse interface matches: {mouse_matches}")
    
    human_matches = interfaces_df['human_interactor_id'].notna().sum()
    logger.info(f"Human interface matches: {human_matches}")
    logger.info(f"Output file: {output_path}")
    logger.info("")


if __name__ == '__main__':
    main()

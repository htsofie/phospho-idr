#!/usr/bin/env python3
"""
Match MoRF residues to proteins in human_ortholog_mapped.csv.

This script:
1. Uses precompiled CSV (default: data/MoRF2/morf_interfaces.csv) 
2. Matches MoRF interface residues to proteins in human_ortholog_mapped.csv (either protein_id or human_ortholog_id)
4. Generates matched_morfs.csv output file with

Usage:
    python scripts/match_morf_interfaces.py --input data/processed/mouse/human_ortholog_mapped.csv --species mouse
    python scripts/match_morf_interfaces.py --input data/processed/rat/human_ortholog_mapped.csv --species rat
"""

import argparse
import logging
import os
from typing import Optional, List

import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


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
    # Keep: protein_id, full_sequence (from full_sequence), length, human_ortholog_id, human_ortholog_sequence
    grouped = df.groupby(['protein_id', 'human_ortholog_id']).first().reset_index()
    
    # Extract required columns
    result_df = pd.DataFrame({
        'protein_id': grouped['protein_id'],
        'full_sequence': grouped.get('full_sequence', grouped.get('sequence', '')),
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
            row_dict["protein_id"] = uniprot_id

            result_rows.append({
                "protein_id": row_dict["protein_id"],
                "full_sequence": row_dict["full_sequence"],
                "length": row_dict.get("length", None),
                "morf_residues": morf_residues,
                "human_ortholog_id": row_dict["human_ortholog_id"],
                "human_ortholog_sequence": row_dict["human_ortholog_sequence"],
                "human_phosphoprotein_IRES": None,
            })
            logger.info(f"  Matched {uniprot_id} to protein_id")
        
        # Check if matches human_ortholog_id
        if uniprot_id in ortholog_by_human.index:
            ortholog_matches = ortholog_by_human.loc[uniprot_id]
            
            # Handle both Series (single row) and DataFrame (multiple rows) cases
            if isinstance(ortholog_matches, pd.Series):
                # Single match - convert to DataFrame for consistent processing
                ortholog_matches = ortholog_matches.to_frame().T
            # Now ortholog_matches is always a DataFrame
            
            # Create result rows for ALL protein_ids that have this human_ortholog_id
            for idx, ortholog_row in ortholog_matches.iterrows():
                # Convert to dict to access all values
                row_dict = ortholog_row.to_dict()
                # Add the human_ortholog_id (which is the index value)
                row_dict["human_ortholog_id"] = uniprot_id
                
                result_rows.append({
                    "protein_id": row_dict["protein_id"],
                    "full_sequence": row_dict["full_sequence"],
                    "length": row_dict.get("length", None),
                    "morf_residues": None,
                    "human_ortholog_id": row_dict["human_ortholog_id"],
                    "human_ortholog_sequence": row_dict["human_ortholog_sequence"],
                    "human_phosphoprotein_IRES": morf_residues,
                })
            
            logger.info(f"  Matched {uniprot_id} to human_ortholog_id ({len(ortholog_matches)} protein_ids)")
    
    result_df = pd.DataFrame(result_rows)
    logger.info(f"  Total matches: {len(result_df)}")
    
    return result_df


def main():
    """Main workflow function.

    This script assumes that MoRF interfaces have already been compiled into a CSV
    (shared between mouse and rat), e.g. via compile_morf_interfaces.py.
    It then matches those MoRFs to a species-specific human_ortholog_mapped.csv.
    """
    parser = argparse.ArgumentParser(
        description="Match MoRF residues to proteins in human_ortholog_mapped.csv using precompiled MoRF data."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to input human_ortholog_mapped.csv file",
    )
    parser.add_argument(
        "--species",
        type=str,
        required=True,
        choices=["mouse", "rat"],
        help="Species: mouse or rat",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to output matched_morfs.csv file (default: same directory as input)",
    )
    parser.add_argument(
        "--morf-csv",
        type=str,
        default=None,
        help="Path to compiled MoRF interfaces CSV (default: data/MoRF2/morf_interfaces.csv)",
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
    
    # Determine MoRF CSV path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    if args.morf_csv:
        morf_csv_path = args.morf_csv
    else:
        morf_csv_path = os.path.join(project_root, "data", "MoRF2", "morf_interfaces.csv")
    
    logger.info("=" * 80)
    logger.info("MoRF Interface Matching (using precompiled MoRF CSV)")
    logger.info("=" * 80)
    logger.info(f"Input file: {args.input}")
    logger.info(f"Species: {args.species}")
    logger.info(f"Output file: {output_path}")
    logger.info(f"MoRF CSV: {morf_csv_path}")
    logger.info("")
    
    # Step 1: Load precompiled MoRF data
    if not os.path.exists(morf_csv_path):
        logger.error(f"MoRF CSV not found: {morf_csv_path}")
        logger.error("Please run compile_morf_interfaces.py first.")
        return

    logger.info("Loading precompiled MoRF data...")
    morf_df = pd.read_csv(morf_csv_path)
    logger.info(f"  Loaded {len(morf_df)} MoRF entries")

    logger.info("")

    # Step 2: Load and group ortholog data
    ortholog_df = load_and_group_ortholog_data(args.input)
    
    logger.info("")
    
    # Step 3: Match MoRF data to ortholog data
    matched_df = match_morf_to_orthologs(morf_df, ortholog_df)
    
    logger.info("")
    
    # Step 4: Save results
    logger.info(f"Saving results to {output_path}...")
    matched_df.to_csv(output_path, index=False)
    logger.info(f"✓ Successfully saved {len(matched_df)} matched rows")
    
    # Print summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    logger.info(f"Total MoRF entries (from CSV): {len(morf_df)}")
    logger.info(f"Total matches: {len(matched_df)}")
    logger.info(f"  - Matches to protein_id: {len(matched_df[matched_df['morf_residues'].notna()])}")
    logger.info(f"  - Matches to human_ortholog_id: {len(matched_df[matched_df['human_phosphoprotein_IRES'].notna()])}")
    logger.info(f"Output file: {output_path}")
    logger.info("")


if __name__ == '__main__':
    main()


#!/usr/bin/env python3
"""
Extract UniProt data for protein IDs from specified columns in a CSV file.

This script queries the UniProt REST API to extract various aspects of proteins
(sequence, length, name, etc.) for protein IDs found in specified columns of an input CSV.

Usage:
    python scripts/extract_uniprot_data.py --input data.csv --columns protein_id,human_ortholog_id --aspects sequence,length,name --output output.csv
    python scripts/extract_uniprot_data.py --input data.csv --columns protein_id --aspects sequence,length --output output.csv
"""

import pandas as pd
import requests
import time
import argparse
import logging
import os
from typing import Dict, Optional, List, Set
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# API rate limiting delay (seconds)
API_DELAY = 0.15


def get_uniprot_entry(protein_id: str) -> Optional[dict]:
    """
    Fetch protein entry from UniProt REST API.
    
    Args:
        protein_id: UniProt protein ID (e.g., 'Q6NZQ8')
        
    Returns:
        Dictionary with UniProt entry data, or None if error
    """
    if pd.isna(protein_id) or not protein_id or str(protein_id).strip() == '':
        return None
    
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id.strip()}"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        logger.debug(f"Error fetching UniProt entry for {protein_id}: {e}")
        return None


def extract_sequence(uniprot_data: dict) -> Optional[str]:
    """Extract full protein sequence."""
    if not uniprot_data:
        return None
    try:
        sequence = uniprot_data.get('sequence', {}).get('value', '')
        return sequence if sequence else None
    except (KeyError, AttributeError):
        return None


def extract_length(uniprot_data: dict) -> Optional[int]:
    """Extract protein sequence length."""
    if not uniprot_data:
        return None
    try:
        length = uniprot_data.get('sequence', {}).get('length', None)
        return int(length) if length is not None else None
    except (KeyError, AttributeError, ValueError):
        return None


def extract_name(uniprot_data: dict) -> Optional[str]:
    """Extract recommended protein name."""
    if not uniprot_data:
        return None
    try:
        name = uniprot_data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '')
        return name if name else None
    except (KeyError, AttributeError):
        return None


def extract_gene_name(uniprot_data: dict) -> Optional[str]:
    """Extract gene name."""
    if not uniprot_data:
        return None
    try:
        genes = uniprot_data.get('genes', [])
        if genes:
            gene_name = genes[0].get('geneName', {}).get('value', '')
            return gene_name if gene_name else None
        return None
    except (KeyError, AttributeError, IndexError):
        return None


def extract_organism(uniprot_data: dict) -> Optional[str]:
    """Extract organism scientific name."""
    if not uniprot_data:
        return None
    try:
        organism = uniprot_data.get('organism', {}).get('scientificName', '')
        return organism if organism else None
    except (KeyError, AttributeError):
        return None


def extract_taxonomy_id(uniprot_data: dict) -> Optional[int]:
    """Extract NCBI taxonomy ID."""
    if not uniprot_data:
        return None
    try:
        taxon_id = uniprot_data.get('organism', {}).get('taxonId', None)
        return int(taxon_id) if taxon_id is not None else None
    except (KeyError, AttributeError, ValueError):
        return None


def extract_mass(uniprot_data: dict) -> Optional[float]:
    """Extract molecular mass (Da)."""
    if not uniprot_data:
        return None
    try:
        mass = uniprot_data.get('sequence', {}).get('molWeight', None)
        return float(mass) if mass is not None else None
    except (KeyError, AttributeError, ValueError):
        return None


def extract_reviewed_status(uniprot_data: dict) -> Optional[bool]:
    """Extract reviewed status (True for Swiss-Prot, False for TrEMBL)."""
    if not uniprot_data:
        return None
    try:
        entry_type = uniprot_data.get('entryType', '')
        return entry_type == 'reviewed'
    except (KeyError, AttributeError):
        return None


# Mapping of aspect names to extraction functions
ASPECT_EXTRACTORS = {
    'sequence': extract_sequence,
    'length': extract_length,
    'name': extract_name,
    'gene_name': extract_gene_name,
    'organism': extract_organism,
    'taxonomy_id': extract_taxonomy_id,
    'mass': extract_mass,
    'reviewed': extract_reviewed_status,
}


def process_protein_id(protein_id: str, aspects: List[str], cache: Dict[str, Dict[str, any]]) -> Dict[str, any]:
    """
    Process a single protein ID and extract requested aspects.
    
    Args:
        protein_id: UniProt protein ID
        aspects: List of aspect names to extract
        cache: Dictionary to cache UniProt entries
        
    Returns:
        Dictionary with extracted aspects
    """
    result = {}
    
    # Check cache first
    if protein_id in cache:
        uniprot_data = cache[protein_id]
    else:
        uniprot_data = get_uniprot_entry(protein_id)
        time.sleep(API_DELAY)
        if uniprot_data:
            cache[protein_id] = uniprot_data
    
    if not uniprot_data:
        # If no data, set all aspects to None
        for aspect in aspects:
            result[aspect] = None
        return result
    
    # Extract requested aspects
    for aspect in aspects:
        if aspect in ASPECT_EXTRACTORS:
            extractor = ASPECT_EXTRACTORS[aspect]
            result[aspect] = extractor(uniprot_data)
        else:
            logger.warning(f"Unknown aspect: {aspect}")
            result[aspect] = None
    
    return result


def main():
    """Main workflow function."""
    parser = argparse.ArgumentParser(
        description='Extract UniProt data for protein IDs from specified columns'
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Path to input CSV file'
    )
    parser.add_argument(
        '--columns',
        type=str,
        required=True,
        help='Comma-separated list of column names containing protein IDs'
    )
    parser.add_argument(
        '--aspects',
        type=str,
        required=True,
        help='Comma-separated list of aspects to extract (sequence, length, name, gene_name, organism, taxonomy_id, mass, reviewed)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Path to output CSV file (default: input file with _uniprot_data suffix)'
    )
    parser.add_argument(
        '--prefix',
        type=str,
        default='uniprot_',
        help='Prefix for output column names (default: uniprot_)'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        logger.error(f"Input file not found: {args.input}")
        return
    
    # Parse columns and aspects
    column_names = [col.strip() for col in args.columns.split(',')]
    aspect_names = [aspect.strip() for aspect in args.aspects.split(',')]
    
    # Validate aspects
    invalid_aspects = [a for a in aspect_names if a not in ASPECT_EXTRACTORS]
    if invalid_aspects:
        logger.error(f"Invalid aspects: {invalid_aspects}")
        logger.info(f"Valid aspects are: {', '.join(ASPECT_EXTRACTORS.keys())}")
        return
    
    # Determine output path
    if args.output:
        output_path = args.output
    else:
        input_base = os.path.splitext(args.input)[0]
        output_path = f"{input_base}_uniprot_data.csv"
    
    logger.info("=" * 80)
    logger.info("UniProt Data Extraction")
    logger.info("=" * 80)
    logger.info(f"Input file: {args.input}")
    logger.info(f"Columns to process: {', '.join(column_names)}")
    logger.info(f"Aspects to extract: {', '.join(aspect_names)}")
    logger.info(f"Output file: {output_path}")
    logger.info(f"Column prefix: {args.prefix}")
    logger.info("")
    
    # Load CSV file
    logger.info("Loading CSV file...")
    try:
        df = pd.read_csv(args.input)
        logger.info(f"Loaded {len(df)} rows")
    except Exception as e:
        logger.error(f"Error reading CSV file: {e}")
        return
    
    # Validate columns exist
    missing_cols = [col for col in column_names if col not in df.columns]
    if missing_cols:
        logger.error(f"Error: Columns not found in CSV: {missing_cols}")
        logger.error(f"Available columns: {list(df.columns)}")
        return
    
    # Collect all unique protein IDs from all specified columns
    logger.info("Collecting unique protein IDs...")
    all_protein_ids = set()
    for col in column_names:
        unique_ids = df[col].dropna().unique()
        all_protein_ids.update([str(id).strip() for id in unique_ids if str(id).strip()])
    
    logger.info(f"Found {len(all_protein_ids)} unique protein IDs to process")
    logger.info("")
    
    # Cache for UniProt entries
    cache = {}
    
    # Process each unique protein ID
    logger.info("Processing protein IDs...")
    protein_data = {}
    processed = 0
    failed = 0
    
    for idx, protein_id in enumerate(sorted(all_protein_ids), 1):
        if idx % 100 == 0:
            logger.info(f"  Processed {idx}/{len(all_protein_ids)} protein IDs...")
        
        result = process_protein_id(protein_id, aspect_names, cache)
        protein_data[protein_id] = result
        
        if any(v is None for v in result.values()):
            failed += 1
        processed += 1
    
    logger.info(f"Processed {processed} protein IDs ({failed} failed to fetch)")
    logger.info("")
    
    # Add new columns to dataframe
    logger.info("Adding extracted data to dataframe...")
    for col in column_names:
        for aspect in aspect_names:
            new_col_name = f"{args.prefix}{col}_{aspect}"
            
            # Create a mapping function
            def get_value(row, source_col=col, aspect_name=aspect):
                protein_id = row[source_col]
                if pd.isna(protein_id) or not protein_id:
                    return None
                protein_id = str(protein_id).strip()
                if protein_id in protein_data:
                    return protein_data[protein_id].get(aspect_name, None)
                return None
            
            df[new_col_name] = df.apply(get_value, axis=1)
    
    # Save to output file
    logger.info(f"Saving results to {output_path}...")
    try:
        df.to_csv(output_path, index=False)
        logger.info(f"✓ Successfully saved {len(df)} rows")
    except Exception as e:
        logger.error(f"Error saving CSV file: {e}")
        return
    
    # Print summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    logger.info(f"Total rows processed: {len(df)}")
    logger.info(f"Unique protein IDs processed: {len(all_protein_ids)}")
    logger.info(f"Successfully fetched: {processed - failed}")
    logger.info(f"Failed to fetch: {failed}")
    logger.info(f"Output file: {output_path}")
    logger.info("")
    logger.info("New columns added:")
    for col in column_names:
        for aspect in aspect_names:
            new_col_name = f"{args.prefix}{col}_{aspect}"
            non_null_count = df[new_col_name].notna().sum()
            logger.info(f"  {new_col_name}: {non_null_count} non-null values")
    logger.info("")


if __name__ == '__main__':
    main()


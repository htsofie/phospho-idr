#!/usr/bin/env python3
"""
Map human orthologs to mouse/rat proteins using Alliance Genome API and UniProt API.

This script processes full_disorder.csv files to add human ortholog information:
- Extracts MGI/RGD gene IDs from UniProt entries
- Queries Alliance Genome API for human orthologs
- Fetches human protein details from UniProt
- Updates all rows with the same protein_id

Usage:
    python scripts/map_human_orthologs.py --input data/processed/mouse/full_disorder.csv --species mouse
    python scripts/map_human_orthologs.py --input data/processed/rat/full_disorder.csv --species rat
"""

import pandas as pd
import requests
import time
import argparse
import logging
import os
from pathlib import Path
from typing import Dict, Optional, Tuple

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
    url = f"https://rest.uniprot.org/uniprotkb/{protein_id}"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        logger.error(f"Error fetching UniProt entry for {protein_id}: {e}")
        return None


def extract_gene_id(uniprot_data: dict, species: str) -> Optional[str]:
    """
    Extract MGI ID (mouse) or RGD ID (rat) from UniProt cross-references.
    
    Args:
        uniprot_data: UniProt entry JSON data
        species: 'mouse' or 'rat'
        
    Returns:
        Gene ID (MGI:XXXXX or RGD:XXXXX), or None if not found
    """
    if not uniprot_data:
        return None
    
    xrefs = uniprot_data.get('uniProtKBCrossReferences', [])
    
    if species == 'mouse':
        # Look for MGI database
        mgi_refs = [x for x in xrefs if x.get('database') == 'MGI']
        if mgi_refs:
            mgi_id = mgi_refs[0].get('id', '')
            return mgi_id if mgi_id else None
    
    elif species == 'rat':
        # Look for RGD database
        rgd_refs = [x for x in xrefs if x.get('database') == 'RGD']
        if rgd_refs:
            rgd_id = rgd_refs[0].get('id', '')
            # Format RGD ID as RGD:XXXXX if needed
            if rgd_id and not rgd_id.startswith('RGD:'):
                return f"RGD:{rgd_id}"
            return rgd_id if rgd_id else None
    
    return None


def get_human_ortholog_from_alliance(gene_id: str) -> Optional[str]:
    """
    Query Alliance Genome API to get human ortholog HGNC ID.
    
    Args:
        gene_id: MGI ID (mouse) or RGD ID (rat)
        
    Returns:
        HGNC ID (e.g., 'HGNC:26077'), or None if not found
    """
    url = f"https://www.alliancegenome.org/api/gene/{gene_id}/orthologs?limit=100"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        # Filter for human orthologs (NCBITaxon:9606)
        results = data.get('results', [])
        human_orthologs = [
            r for r in results
            if r.get('geneToGeneOrthologyGenerated', {}).get('objectGene', {}).get('taxon', {}).get('curie') == 'NCBITaxon:9606'
        ]
        
        if human_orthologs:
            hgnc_id = human_orthologs[0].get('geneToGeneOrthologyGenerated', {}).get('objectGene', {}).get('primaryExternalId', '')
            return hgnc_id if hgnc_id else None
        
        return None
        
    except requests.exceptions.RequestException as e:
        logger.error(f"Error querying Alliance Genome API for {gene_id}: {e}")
        return None


def get_human_uniprot_id(hgnc_id: str) -> Optional[str]:
    """
    Search UniProt for human protein using HGNC ID.
    
    Args:
        hgnc_id: HGNC gene ID (e.g., 'HGNC:26077')
        
    Returns:
        Human UniProt ID, or None if not found
    """
    # Search for reviewed human proteins with this HGNC ID
    query = f"database:hgnc+AND+{hgnc_id}+AND+organism_id:9606+AND+reviewed:true"
    url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=json&size=1"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        results = data.get('results', [])
        if results:
            return results[0].get('primaryAccession', '')
        
        return None
        
    except requests.exceptions.RequestException as e:
        logger.error(f"Error searching UniProt for {hgnc_id}: {e}")
        return None


def get_protein_info(uniprot_id: str) -> Optional[str]:
    """
    Get protein name/description from UniProt.
    
    Args:
        uniprot_id: UniProt protein ID
        
    Returns:
        Protein description, or None if error
    """
    uniprot_data = get_uniprot_entry(uniprot_id)
    if not uniprot_data:
        return None
    
    try:
        protein_desc = uniprot_data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '')
        return protein_desc if protein_desc else None
    except (KeyError, AttributeError):
        return None


def get_protein_sequence(uniprot_id: str) -> Optional[str]:
    """
    Get full sequence from UniProt.
    
    Args:
        uniprot_id: UniProt protein ID
        
    Returns:
        Protein sequence, or None if error
    """
    uniprot_data = get_uniprot_entry(uniprot_id)
    if not uniprot_data:
        return None
    
    try:
        sequence = uniprot_data.get('sequence', {}).get('value', '')
        return sequence if sequence else None
    except (KeyError, AttributeError):
        return None


def process_protein_group(group_df: pd.DataFrame, protein_id: str, species: str, cache: dict) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str], Optional[str]]:
    """
    Process one protein group to get human ortholog information.
    
    Args:
        group_df: DataFrame with all rows for this protein_id
        protein_id: UniProt protein ID
        species: 'mouse' or 'rat'
        cache: Dictionary to cache API results (key: protein_id, value: (ortholog_id, ortholog_desc, ortholog_seq, source_desc, error))
        
    Returns:
        Tuple of (human_ortholog_id, human_ortholog_description, human_ortholog_sequence, source_protein_description, api_error)
    """
    # Check cache first
    if protein_id in cache:
        return cache[protein_id]
    
    api_error = None
    human_ortholog_id = None
    human_ortholog_description = None
    human_ortholog_sequence = None
    source_protein_description = None
    
    try:
        # Step 1: Get UniProt entry and extract gene ID
        logger.info(f"Processing {protein_id} ({species})...")
        uniprot_data = get_uniprot_entry(protein_id)
        time.sleep(API_DELAY)
        
        if not uniprot_data:
            api_error = f"Failed to fetch UniProt entry for {protein_id}"
            result = (None, None, None, None, api_error)
            cache[protein_id] = result
            return result
        
        # Extract source protein description (always do this, regardless of ortholog mapping success)
        try:
            source_protein_description = uniprot_data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', '')
            if source_protein_description:
                logger.info(f"  Source protein: {source_protein_description}")
        except (KeyError, AttributeError):
            pass
        
        gene_id = extract_gene_id(uniprot_data, species)
        if not gene_id:
            api_error = f"No {species} gene ID (MGI/RGD) found in UniProt entry for {protein_id}"
            result = (None, None, None, source_protein_description, api_error)
            cache[protein_id] = result
            return result
        
        logger.info(f"  Found gene ID: {gene_id}")
        time.sleep(API_DELAY)
        
        # Step 2: Query Alliance Genome for human ortholog
        hgnc_id = get_human_ortholog_from_alliance(gene_id)
        time.sleep(API_DELAY)
        
        if not hgnc_id:
            api_error = f"No human ortholog found in Alliance Genome for {gene_id}"
            result = (None, None, None, source_protein_description, api_error)
            cache[protein_id] = result
            return result
        
        logger.info(f"  Found human HGNC ID: {hgnc_id}")
        time.sleep(API_DELAY)
        
        # Step 3: Search UniProt for human protein
        human_uniprot_id = get_human_uniprot_id(hgnc_id)
        time.sleep(API_DELAY)
        
        if not human_uniprot_id:
            api_error = f"No human UniProt ID found for {hgnc_id}"
            result = (None, None, None, source_protein_description, api_error)
            cache[protein_id] = result
            return result
        
        logger.info(f"  Found human UniProt ID: {human_uniprot_id}")
        human_ortholog_id = human_uniprot_id
        time.sleep(API_DELAY)
        
        # Step 4: Get human protein details
        human_ortholog_description = get_protein_info(human_uniprot_id)
        time.sleep(API_DELAY)
        
        human_ortholog_sequence = get_protein_sequence(human_uniprot_id)
        time.sleep(API_DELAY)
        
        if human_ortholog_description:
            logger.info(f"  Human protein: {human_ortholog_description}")
        if human_ortholog_sequence:
            logger.info(f"  Human sequence length: {len(human_ortholog_sequence)} amino acids")
        
        result = (human_ortholog_id, human_ortholog_description, human_ortholog_sequence, source_protein_description, None)
        cache[protein_id] = result
        return result
        
    except Exception as e:
        api_error = f"Unexpected error processing {protein_id}: {str(e)}"
        logger.error(api_error)
        result = (None, None, None, source_protein_description, api_error)
        cache[protein_id] = result
        return result


def main():
    """Main workflow function."""
    parser = argparse.ArgumentParser(
        description='Map human orthologs to mouse/rat proteins using Alliance Genome and UniProt APIs'
    )
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Path to input full_disorder.csv file'
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
        help='Path to output CSV file (default: {input_dir}/{species}/human_ortholog_mapped.csv)'
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
        input_dir = os.path.dirname(os.path.dirname(args.input))  # Go up from processed/{species}/
        output_path = os.path.join(input_dir, args.species, 'human_ortholog_mapped.csv')
    
    # Create output directory if needed
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    logger.info("=" * 80)
    logger.info("Human Ortholog Mapping")
    logger.info("=" * 80)
    logger.info(f"Input file: {args.input}")
    logger.info(f"Species: {args.species}")
    logger.info(f"Output file: {output_path}")
    logger.info("")
    
    # Load CSV file
    logger.info(f"Loading CSV file...")
    try:
        df = pd.read_csv(args.input)
        logger.info(f"Loaded {len(df)} rows")
    except Exception as e:
        logger.error(f"Error reading CSV file: {e}")
        return
    
    # Check required columns
    if 'protein_id' not in df.columns:
        logger.error("Error: 'protein_id' column not found in CSV file")
        logger.error(f"Available columns: {list(df.columns)}")
        return
    
    # Initialize new columns
    df['human_ortholog_id'] = None
    df['human_ortholog_description'] = None
    df['human_ortholog_sequence'] = None
    df['api_error'] = None
    
    # Also store UniProt description for source protein
    df['uniprot_protein_description'] = None
    
    # Group by protein_id
    unique_proteins = df['protein_id'].unique()
    logger.info(f"Found {len(unique_proteins)} unique protein IDs")
    logger.info("")
    
    # Cache to avoid duplicate API calls
    cache = {}
    
    # Process each unique protein_id
    successful = 0
    failed = 0
    
    for idx, protein_id in enumerate(unique_proteins, 1):
        logger.info(f"[{idx}/{len(unique_proteins)}] Processing protein: {protein_id}")
        
        # Get group of rows with this protein_id
        group_mask = df['protein_id'] == protein_id
        group_df = df[group_mask]
        
        # Process the protein group
        human_ortholog_id, human_ortholog_description, human_ortholog_sequence, source_protein_description, api_error = process_protein_group(
            group_df, protein_id, args.species, cache
        )
        
        # Update all rows in the group
        df.loc[group_mask, 'human_ortholog_id'] = human_ortholog_id
        df.loc[group_mask, 'human_ortholog_description'] = human_ortholog_description
        df.loc[group_mask, 'human_ortholog_sequence'] = human_ortholog_sequence
        df.loc[group_mask, 'api_error'] = api_error
        
        # Always update source protein description from UniProt (regardless of ortholog mapping success)
        if source_protein_description:
            df.loc[group_mask, 'uniprot_protein_description'] = source_protein_description
        
        if api_error:
            failed += 1
        else:
            successful += 1
        
        logger.info("")
    
    # Replace protein_description with UniProt description if available
    logger.info("Updating protein_description column with UniProt descriptions...")
    mask = df['uniprot_protein_description'].notna()
    df.loc[mask, 'protein_description'] = df.loc[mask, 'uniprot_protein_description']
    
    # Remove temporary column
    df = df.drop(columns=['uniprot_protein_description'])
    
    # Reorder columns: insert new columns after protein_description
    cols = list(df.columns)
    
    # Remove new columns from their current positions
    for col in ['human_ortholog_id', 'human_ortholog_description', 'human_ortholog_sequence', 'api_error']:
        if col in cols:
            cols.remove(col)
    
    # Find position of protein_description
    if 'protein_description' in cols:
        protein_desc_idx = cols.index('protein_description')
        # Insert new columns after protein_description
        cols.insert(protein_desc_idx + 1, 'human_ortholog_id')
        cols.insert(protein_desc_idx + 2, 'human_ortholog_description')
        cols.insert(protein_desc_idx + 3, 'human_ortholog_sequence')
    else:
        # If protein_description not found, add new columns at the beginning
        cols = ['human_ortholog_id', 'human_ortholog_description', 'human_ortholog_sequence'] + cols
    
    # Add api_error at the end
    cols.append('api_error')
    
    # Reorder dataframe
    df = df[cols]
    
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
    logger.info(f"Total proteins processed: {len(unique_proteins)}")
    logger.info(f"Successful mappings: {successful}")
    logger.info(f"Failed mappings: {failed}")
    logger.info(f"Output file: {output_path}")
    logger.info("")


if __name__ == '__main__':
    main()


#!/usr/bin/env python3
"""
Script to map various protein identifiers to UniProt IDs using BioPython.
Tries in order: uniprot_id → swissprot_id → trembl_id → ensembl_id
"""

import pandas as pd
import time
import argparse
import os
import requests
from typing import Optional, Dict, Any
import yaml
from Bio import ExPASy
from Bio import SwissProt


def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)


def map_refseq_to_uniprot(refseq_id: str) -> Optional[str]:
    """Map RefSeq protein ID to UniProt ID using UniProt REST API."""
    if pd.isna(refseq_id) or not refseq_id:
        return None
    
    try:
        clean_id = str(refseq_id).strip()
        
        # Handle multiple IDs separated by semicolons
        if ';' in clean_id:
            ids = [id.strip() for id in clean_id.split(';') if id.strip()]
            for single_id in ids:
                result = _map_single_refseq_id(single_id)
                if result:
                    return result
            return None
        else:
            return _map_single_refseq_id(clean_id)
    except Exception as e:
        print(f"  ✗ RefSeq mapping {refseq_id}: {e}")
        return None


def _map_single_refseq_id(refseq_id: str) -> Optional[str]:
    """Map a single RefSeq ID to UniProt ID using BioPython and REST API fallback."""
    
    # Approach 1: Try to find RefSeq in UniProt via cross-reference search (REST API)
    try:
        url = f"https://rest.uniprot.org/uniprotkb/search?query=xref:{refseq_id}&format=json&size=1"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            if data.get('results') and len(data['results']) > 0:
                result = data['results'][0]
                primary_accession = result.get('primaryAccession')
                if primary_accession:
                    return primary_accession
    except Exception as e:
        print(f"  ✗ UniProt cross-ref search failed for {refseq_id}: {e}")
    
    # Approach 2: Try direct RefSeq search in UniProt (REST API)
    try:
        url = f"https://rest.uniprot.org/uniprotkb/search?query={refseq_id}&format=json&size=1"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            if data.get('results') and len(data['results']) > 0:
                result = data['results'][0]
                primary_accession = result.get('primaryAccession')
                if primary_accession:
                    return primary_accession
    except Exception as e:
        print(f"  ✗ UniProt direct search failed for {refseq_id}: {e}")
    
    return None


def map_ensembl_to_uniprot(ensembl_id: str) -> Optional[str]:
    """Map Ensembl protein ID to UniProt ID using Ensembl REST API."""
    if pd.isna(ensembl_id) or not ensembl_id:
        return None
    
    try:
        clean_id = str(ensembl_id).strip()
        
        # Handle multiple IDs separated by semicolons
        if ';' in clean_id:
            ids = [id.strip() for id in clean_id.split(';') if id.strip()]
            for single_id in ids:
                result = _map_single_ensembl_id(single_id)
                if result:
                    return result
            return None
        else:
            return _map_single_ensembl_id(clean_id)
    except Exception as e:
        print(f"  ✗ Ensembl mapping {ensembl_id}: {e}")
        return None


def _map_single_ensembl_id(ensembl_id: str) -> Optional[str]:
    """Map a single Ensembl ID to UniProt ID."""
    try:
        # Get all cross-references for this Ensembl ID
        url = f"https://rest.ensembl.org/xrefs/id/{ensembl_id}"
        response = requests.get(url, headers={"Content-Type": "application/json"}, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            if data:
                # Look for UniProt/SWISS-PROT entries
                for entry in data:
                    dbname = entry.get('dbname', '')
                    if dbname in ['Uniprot/SWISSPROT', 'Uniprot/SPTREMBL', 'UniProt']:
                        uniprot_id = entry.get('primary_id')
                        if uniprot_id:
                            return uniprot_id
        return None
    except Exception as e:
        print(f"  ✗ Single Ensembl mapping {ensembl_id}: {e}")
        return None


def map_trembl_to_uniprot(trembl_id: str) -> Optional[str]:
    """Map TREMBL ID to UniProt ID using ExPASy."""
    if pd.isna(trembl_id) or not trembl_id:
        return None
    
    try:
        clean_id = str(trembl_id).strip()
        
        # Handle multiple IDs separated by semicolons
        if ';' in clean_id:
            ids = [id.strip() for id in clean_id.split(';') if id.strip()]
            for single_id in ids:
                try:
                    handle = ExPASy.get_sprot_raw(single_id)
                    record = SwissProt.read(handle)
                    handle.close()
                    return single_id  # Return the first valid ID
                except:
                    continue  # Try next ID
            return None
        else:
            # Single ID
            handle = ExPASy.get_sprot_raw(clean_id)
            record = SwissProt.read(handle)
            handle.close()
            return clean_id
    except Exception as e:
        print(f"  ✗ TREMBL mapping {trembl_id}: {e}")
        return None


def map_swissprot_to_uniprot(swissprot_id: str) -> Optional[str]:
    """Map SWISS-PROT ID to UniProt ID using ExPASy."""
    if pd.isna(swissprot_id) or not swissprot_id:
        return None
    
    try:
        clean_id = str(swissprot_id).strip()
        
        # Handle multiple IDs separated by semicolons
        if ';' in clean_id:
            ids = [id.strip() for id in clean_id.split(';') if id.strip()]
            for single_id in ids:
                try:
                    handle = ExPASy.get_sprot_raw(single_id)
                    record = SwissProt.read(handle)
                    handle.close()
                    return single_id  # Return the first valid ID
                except:
                    continue  # Try next ID
            return None
        else:
            # Single ID
            handle = ExPASy.get_sprot_raw(clean_id)
            record = SwissProt.read(handle)
            handle.close()
            return clean_id
    except Exception as e:
        print(f"  ✗ SWISS-PROT mapping {swissprot_id}: {e}")
        return None


def map_uniprot_id_to_uniprot(uniprot_id: str) -> Optional[str]:
    """Verify and return UniProt ID if it exists using ExPASy."""
    if pd.isna(uniprot_id) or not uniprot_id:
        return None
    
    try:
        clean_id = str(uniprot_id).strip()
        
        # Handle multiple IDs separated by semicolons
        if ';' in clean_id:
            ids = [id.strip() for id in clean_id.split(';') if id.strip()]
            for single_id in ids:
                try:
                    handle = ExPASy.get_sprot_raw(single_id)
                    record = SwissProt.read(handle)
                    handle.close()
                    return single_id  # Return the first valid ID
                except:
                    continue  # Try next ID
            return None
        else:
            # Single ID
            handle = ExPASy.get_sprot_raw(clean_id)
            record = SwissProt.read(handle)
            handle.close()
            return clean_id
    except Exception as e:
        print(f"  ✗ UniProt ID mapping {uniprot_id}: {e}")
        return None


def map_protein_to_uniprot(row: pd.Series, species: str) -> Dict[str, Any]:
    """
    Map a protein row to UniProt ID using the priority order:
    uniprot_id → swissprot_id → trembl_id → ensembl_id → refseq_id
    """
    result = {
        'uniprot_mapped_id': None,
        'mapping_source': None,
        'mapping_success': False
    }
    
    # Priority order for mapping
    if species == 'mouse':
        # For mouse data, we have swissprot_id, trembl_id, ensembl columns
        id_columns = [
            ('swissprot_id', 'SWISS-PROT'),
            ('trembl_id', 'TREMBL'),
            ('ensembl', 'ENSEMBL')
        ]
    else:
        # For rat data, we have uniprot, ensembl columns
        id_columns = [
            ('uniprot', 'UniProt'),
            ('ensembl', 'ENSEMBL')
        ]
    
    # Try the main ID columns first
    for col_name, source_name in id_columns:
        if col_name in row.index and not pd.isna(row[col_name]) and row[col_name]:
            id_value = str(row[col_name]).strip()
            if not id_value:
                continue
                
            print(f"  Trying {source_name}: {id_value}")
            
            # Map based on source type
            if source_name == 'UniProt':
                mapped_id = map_uniprot_id_to_uniprot(id_value)
            elif source_name == 'SWISS-PROT':
                mapped_id = map_swissprot_to_uniprot(id_value)
            elif source_name == 'TREMBL':
                mapped_id = map_trembl_to_uniprot(id_value)
            elif source_name == 'ENSEMBL':
                mapped_id = map_ensembl_to_uniprot(id_value)
            else:
                continue
            
            if mapped_id:
                result['uniprot_mapped_id'] = mapped_id
                result['mapping_source'] = source_name
                result['mapping_success'] = True
                print(f"  ✓ Mapped via {source_name}: {mapped_id}")
                return result
            else:
                print(f"  ✗ Failed to map via {source_name}")
    
    # If no main ID worked, try RefSeq as fallback
    # Extract RefSeq IDs from the original protein_IDs column if available
    if 'protein_IDs' in row.index and not pd.isna(row['protein_IDs']):
        protein_ids = str(row['protein_IDs'])
        if 'REFSEQ:' in protein_ids:
            # Extract RefSeq IDs
            refseq_ids = []
            parts = protein_ids.split('|')
            for part in parts:
                if part.startswith('REFSEQ:'):
                    refseq_id = part.split(':', 1)[1].strip()
                    if refseq_id:
                        refseq_ids.append(refseq_id)
            
            if refseq_ids:
                print(f"  Trying RefSeq: {refseq_ids[0]}")
                mapped_id = map_refseq_to_uniprot(refseq_ids[0])
                if mapped_id:
                    result['uniprot_mapped_id'] = mapped_id
                    result['mapping_source'] = 'RefSeq'
                    result['mapping_success'] = True
                    print(f"  ✓ Mapped via RefSeq: {mapped_id}")
                    return result
                else:
                    print(f"  ✗ Failed to map via RefSeq")
    
    if not result['mapping_success']:
        print(f"  ✗ No successful mapping for row")
    
    return result


def process_dataset(input_path: str, species: str, output_path: str) -> None:
    """Process the dataset and add UniProt mapping columns."""
    print(f"Loading dataset: {input_path}")
    df = pd.read_csv(input_path)
    print(f"Loaded {len(df)} rows")
    
    # Add new columns for mapping results
    df['uniprot_mapped_id'] = None
    df['mapping_source'] = None
    df['mapping_success'] = False
    
    print(f"Starting UniProt mapping for {species} data...")
    
    # Process each row
    for idx, row in df.iterrows():
        print(f"\nRow {idx + 1}/{len(df)}")
        
        # Map protein to UniProt ID
        mapping_result = map_protein_to_uniprot(row, species)
        
        # Update the dataframe
        df.at[idx, 'uniprot_mapped_id'] = mapping_result['uniprot_mapped_id']
        df.at[idx, 'mapping_source'] = mapping_result['mapping_source']
        df.at[idx, 'mapping_success'] = mapping_result['mapping_success']
        
        # Add a small delay to be respectful to the APIs
        time.sleep(0.1)
    
    # Reorder columns to put mapping results at the beginning
    mapping_cols = ['uniprot_mapped_id', 'mapping_source', 'mapping_success']
    other_cols = [col for col in df.columns if col not in mapping_cols]
    df_reordered = df[mapping_cols + other_cols]
    
    # Save the results
    print(f"\nSaving results to: {output_path}")
    df_reordered.to_csv(output_path, index=False)
    
    # Print summary statistics
    successful_mappings = df['mapping_success'].sum()
    print(f"\nMapping Summary:")
    print(f"  Total rows: {len(df)}")
    print(f"  Successful mappings: {successful_mappings}")
    print(f"  Failed mappings: {len(df) - successful_mappings}")
    print(f"  Success rate: {successful_mappings/len(df)*100:.1f}%")
    
    # Print mapping source breakdown
    if successful_mappings > 0:
        print(f"\nMapping source breakdown:")
        source_counts = df[df['mapping_success']]['mapping_source'].value_counts()
        for source, count in source_counts.items():
            print(f"  {source}: {count} ({count/successful_mappings*100:.1f}%)")


def main():
    parser = argparse.ArgumentParser(description='Map protein identifiers to UniProt IDs using BioPython, Ensembl REST API, and UniProt REST API')
    parser.add_argument('--input', '-i', required=True, help='Path to input CSV file (cleaned data)')
    parser.add_argument('--species', '-s', required=True, choices=['mouse', 'rat'], help='Species name')
    parser.add_argument('--output', '-o', help='Path to output CSV file (default: data/processed/{species}/input_filename_mapped.csv)')
    
    args = parser.parse_args()
    
    # Determine output path
    if args.output:
        output_path = args.output
    else:
        # Create output directory
        output_dir = os.path.join('data', 'processed', args.species)
        os.makedirs(output_dir, exist_ok=True)
        
        # Get input filename and create mapped filename
        input_filename = os.path.basename(args.input)
        base_name = os.path.splitext(input_filename)[0]
        output_filename = f"{base_name}_mapped.csv"
        output_path = os.path.join(output_dir, output_filename)
    
    # Process the dataset
    process_dataset(args.input, args.species, output_path)


if __name__ == "__main__":
    main()

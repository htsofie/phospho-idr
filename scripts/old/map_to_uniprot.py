#!/usr/bin/env python3
"""
Script to map various protein identifiers to UniProt IDs using BioPython.
Tries in order: UniProt → SWISS-PROT → TREMBL. ENSEMBL is not used.
"""

import pandas as pd
import time
import argparse
import os
import requests
from typing import Optional, Dict, Any, Tuple
import yaml
from Bio import ExPASy
from Bio import SwissProt


def load_config(config_path: str) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)


## Removed RefSeq mapping to comply with "only UniProt/Swiss-Prot/TREMBL"


## Removed ENSEMBL mapping per request


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
    UniProt → SWISS-PROT → TREMBL
    Verifies species match and tries other IDs if species doesn't match.
    """
    result = {
        'uniprot_mapped_id': None,
        'mapping_source': None,
        'mapping_success': False,
        'species_mismatch': False
    }
    
    # Priority order for mapping
    if species == 'mouse':
        # For mouse data, prefer direct UniProt if present; else Swiss-Prot, then TREMBL
        id_columns = [
            ('uniprot', 'UniProt'),
            ('swissprot_id', 'SWISS-PROT'),
            ('trembl_id', 'TREMBL')
        ]
    else:
        # For rat data, use UniProt only
        id_columns = [
            ('uniprot', 'UniProt')
        ]
    
    # Try the main ID columns first
    for col_name, source_name in id_columns:
        if col_name in row.index and not pd.isna(row[col_name]) and row[col_name]:
            id_value = str(row[col_name]).strip()
            if not id_value:
                continue
                
            print(f"  Trying {source_name}: {id_value}")
            
            # Handle multiple IDs separated by semicolons
            if ';' in id_value:
                ids = [id.strip() for id in id_value.split(';') if id.strip()]
                for single_id in ids:
                    print(f"    Trying individual ID: {single_id}")
                    
                    # Map based on source type
                    if source_name == 'UniProt':
                        mapped_id = map_uniprot_id_to_uniprot(single_id)
                    elif source_name == 'SWISS-PROT':
                        mapped_id = map_swissprot_to_uniprot(single_id)
                    elif source_name == 'TREMBL':
                        mapped_id = map_trembl_to_uniprot(single_id)
                    else:
                        continue
                    
                    if mapped_id:
                        # Verify species match
                        if verify_species_match(mapped_id, species):
                            result['uniprot_mapped_id'] = mapped_id
                            result['mapping_source'] = source_name
                            result['mapping_success'] = True
                            result['species_mismatch'] = False
                            print(f"  ✓ Mapped via {source_name}: {mapped_id} (species verified)")
                            return result
                        else:
                            print(f"  ⚠ Species mismatch for {mapped_id}, trying next ID...")
                            continue
                    else:
                        print(f"    ✗ Failed to map {single_id}")
            else:
                # Single ID
                # Map based on source type
                if source_name == 'UniProt':
                    mapped_id = map_uniprot_id_to_uniprot(id_value)
                elif source_name == 'SWISS-PROT':
                    mapped_id = map_swissprot_to_uniprot(id_value)
                elif source_name == 'TREMBL':
                    mapped_id = map_trembl_to_uniprot(id_value)
                else:
                    continue
                
                if mapped_id:
                    # Verify species match
                    if verify_species_match(mapped_id, species):
                        result['uniprot_mapped_id'] = mapped_id
                        result['mapping_source'] = source_name
                        result['mapping_success'] = True
                        result['species_mismatch'] = False
                        print(f"  ✓ Mapped via {source_name}: {mapped_id} (species verified)")
                        return result
                    else:
                        print(f"  ⚠ Species mismatch for {mapped_id}, trying next source...")
                        result['species_mismatch'] = True
                        continue
                else:
                    print(f"  ✗ Failed to map via {source_name}")
    
    # If we get here, either no mapping worked or all had species mismatches
    if result['species_mismatch']:
        print(f"  ⚠ All mappings had species mismatches")
    else:
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
    df['species_mismatch'] = False
    
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
        df.at[idx, 'species_mismatch'] = mapping_result['species_mismatch']
        
        # Add a small delay to be respectful to the APIs
        time.sleep(0.1)
    
    # Reorder columns to put mapping results at the beginning
    mapping_cols = ['uniprot_mapped_id', 'mapping_source', 'mapping_success', 'species_mismatch']
    other_cols = [col for col in df.columns if col not in mapping_cols]
    df_reordered = df[mapping_cols + other_cols]
    
    # Save the results
    print(f"\nSaving results to: {output_path}")
    df_reordered.to_csv(output_path, index=False)
    
    # Print summary statistics
    successful_mappings = df['mapping_success'].sum()
    species_mismatches = df['species_mismatch'].sum()
    print(f"\nMapping Summary:")
    print(f"  Total rows: {len(df)}")
    print(f"  Successful mappings: {successful_mappings}")
    print(f"  Species mismatches: {species_mismatches}")
    print(f"  Failed mappings: {len(df) - successful_mappings - species_mismatches}")
    print(f"  Success rate: {successful_mappings/len(df)*100:.1f}%")
    print(f"  Species mismatch rate: {species_mismatches/len(df)*100:.1f}%")
    
    # Print mapping source breakdown
    if successful_mappings > 0:
        print(f"\nMapping source breakdown:")
        source_counts = df[df['mapping_success']]['mapping_source'].value_counts()
        for source, count in source_counts.items():
            print(f"  {source}: {count} ({count/successful_mappings*100:.1f}%)")


def main():
    parser = argparse.ArgumentParser(description='Map protein identifiers to UniProt IDs using BioPython (UniProt, SWISS-PROT, TREMBL only)')
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

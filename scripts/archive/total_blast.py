#!/usr/bin/env python3
"""
Total BLAST phosphosite mapping using full UniProtKB species databases.

This script filters for rows with manual_review == TRUE, then BLASTs the first
sequence of each protein group against the complete species UniProtKB database.

Usage:
    python total_blast.py --input data/processed/rat/blast_results.csv --species rat
"""

import argparse
import os
import sys
import subprocess
import tempfile
import pandas as pd
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML

# Add scripts directory to path
sys.path.append(os.path.join(os.path.dirname(__file__)))

# Set up logging
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Excel cell limit (32,767 characters)
EXCEL_CELL_LIMIT = 32767
# Maximum number of matches to extract per protein group
MAX_MATCHES = 10
# Maximum number of IDs to search per group
MAX_IDS_PER_GROUP = 50
# Maximum number of sequences to BLAST per group
MAX_SEQUENCES_PER_GROUP = 20

# Species name mapping (common name -> scientific name)
SPECIES_MAPPING = {
    'mouse': 'Mus musculus',
    'rat': 'Rattus norvegicus',
    'human': 'Homo sapiens'
}


def load_full_database(species: str) -> Dict[str, Dict[str, Any]]:
    """Load the full UniProtKB species database FASTA file and parse into dictionary."""
    # Map species to taxonomy ID filenames
    species_files = {
        'mouse': 'uniprotkb_taxonomy_id_10090_2025_10_06.fasta',
        'rat': 'uniprotkb_taxonomy_id_10116_2025_10_06.fasta'
    }
    
    if species not in species_files:
        raise ValueError(f"Unsupported species: {species}")
    
    db_path = f"data/blast_dbs/{species_files[species]}"
    
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Full UniProtKB database not found: {db_path}")
    
    logger.info(f"Loading full UniProtKB database: {db_path}")
    
    database = {}
    try:
        # Read the file line by line to handle multi-line headers
        with open(db_path, 'r') as f:
            lines = f.readlines()
        
        current_header = ""
        current_sequence = []
        
        for line in lines:
            line = line.strip()
            if line.startswith('>'):
                # Process previous entry if exists
                if current_header and current_sequence:
                    sequence = ''.join(current_sequence)
                    # Parse header
                    if current_header.startswith('>sp|'):
                        db_type = 'swissprot'
                        parts = current_header.split('|')
                        if len(parts) >= 2:
                            uniprot_id = parts[1]
                        else:
                            current_header = line
                            current_sequence = []
                            continue
                    elif current_header.startswith('>tr|'):
                        db_type = 'trembl'
                        parts = current_header.split('|')
                        if len(parts) >= 2:
                            uniprot_id = parts[1]
                        else:
                            current_header = line
                            current_sequence = []
                            continue
                    else:
                        current_header = line
                        current_sequence = []
                        continue
                    
                    # Extract species information
                    species_match = re.search(r'OS=([^=]+)', current_header)
                    if species_match:
                        protein_species = species_match.group(1).strip()
                        
                        # Store in database
                        database[uniprot_id] = {
                            'sequence': sequence,
                            'length': len(sequence),
                            'db_type': db_type,
                            'species': protein_species,
                            'header': current_header
                        }
                
                # Start new entry
                current_header = line
                current_sequence = []
            else:
                # Add sequence line
                current_sequence.append(line)
        
        # Process last entry
        if current_header and current_sequence:
            sequence = ''.join(current_sequence)
            # Parse header
            if current_header.startswith('>sp|'):
                db_type = 'swissprot'
                parts = current_header.split('|')
                if len(parts) >= 2:
                    uniprot_id = parts[1]
                else:
                    pass
            elif current_header.startswith('>tr|'):
                db_type = 'trembl'
                parts = current_header.split('|')
                if len(parts) >= 2:
                    uniprot_id = parts[1]
                else:
                    pass
            else:
                pass
            
            # Extract species information
            species_match = re.search(r'OS=([^=]+)', current_header)
            if species_match:
                protein_species = species_match.group(1).strip()
                
                # Store in database
                database[uniprot_id] = {
                    'sequence': sequence,
                    'length': len(sequence),
                    'db_type': db_type,
                    'species': protein_species,
                    'header': current_header
                }
    
    except Exception as e:
        logger.error(f"Failed to load database: {e}")
        raise
    
    logger.info(f"Loaded {len(database)} proteins from paper database")
    return database


def create_blast_database(fasta_path: str, db_path: str) -> bool:
    """Create BLAST database from FASTA file."""
    logger.info(f"Creating BLAST database from {fasta_path}")
    
    try:
        cmd = ['makeblastdb', '-in', fasta_path, '-dbtype', 'prot', '-out', db_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"makeblastdb failed: {result.stderr}")
            return False
        
        logger.info("BLAST database created successfully")
        return True
        
    except Exception as e:
        logger.error(f"Failed to create BLAST database: {e}")
        return False


def blast_sequence(sequence: str, db_path: str) -> List[Dict[str, Any]]:
    """Run BLAST search and return hits."""
    logger.info(f"BLASTing sequence: {sequence[:13]}...")
    
    try:
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
            f.write(f">query\n{sequence}\n")
            temp_input = f.name
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
            temp_output = f.name
        
        # Run BLAST with parameters optimized for short sequences
        cmd = [
            'blastp', '-query', temp_input, '-db', db_path, '-out', temp_output,
            '-outfmt', '5', '-evalue', '1000', '-word_size', '2',
            '-gapopen', '9', '-gapextend', '1', '-matrix', 'BLOSUM62',
            '-comp_based_stats', '0', '-max_target_seqs', '50'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"BLAST failed: {result.stderr}")
            return []
        
        # Parse results
        hits = []
        try:
            with open(temp_output, 'r') as f:
                blast_records = NCBIXML.parse(f)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            # Extract UniProt ID from title - need to parse the FASTA header format
                            title = alignment.title
                            logger.debug(f"BLAST title: {title}")
                            
                            # Parse FASTA header: gnl|BL_ORD_ID|2388 tr|M0R9X8|M0R9X8_RAT
                            # We need to extract the UniProt ID from the tr| or sp| part
                            if 'tr|' in title:
                                # TrEMBL format: tr|UNIPROT_ID|PROTEIN_NAME
                                parts = title.split('tr|')[1].split('|')
                                uniprot_id = parts[0]
                            elif 'sp|' in title:
                                # SwissProt format: sp|UNIPROT_ID|PROTEIN_NAME
                                parts = title.split('sp|')[1].split('|')
                                uniprot_id = parts[0]
                            elif '|' in title:
                                # Generic format with pipes
                                parts = title.split('|')
                                if len(parts) >= 2:
                                    uniprot_id = parts[1]
                                else:
                                    uniprot_id = parts[0]
                            else:
                                uniprot_id = title.split()[0]
                            
                            logger.debug(f"Extracted UniProt ID: {uniprot_id}")
                            
                            hit = {
                                'protein_id': uniprot_id,
                                'identity_percent': (hsp.identities / hsp.align_length) * 100,
                                'evalue': hsp.expect,
                                'bit_score': hsp.bits,
                                'hsp_query': str(hsp.query),
                                'hsp_subject': str(hsp.sbjct),
                                'hsp_sbjct_start': int(hsp.sbjct_start),
                                'hsp_sbjct_end': int(hsp.sbjct_end),
                                'hsp_align_length': int(hsp.align_length),
                                'hsp_identities': int(hsp.identities)
                            }
                            # Only accept exact full-length matches: 100% identity over the entire query length
                            if hit["identity_percent"] == 100.0 and int(hsp.align_length) == len(sequence):
                                hits.append(hit)
        except Exception as e:
            logger.error(f"Failed to parse BLAST results: {e}")
        
        # Clean up
        os.unlink(temp_input)
        os.unlink(temp_output)
        
        logger.info(f"Found {len(hits)} BLAST hits")
        return hits
        
    except Exception as e:
        logger.error(f"BLAST failed: {e}")
        return []


def check_cell_limit(ids: List[str], sequences: List[str], lengths: List[int], types: List[str]) -> Tuple[List[str], List[str], List[int], List[str]]:
    """Limit results to fit within Excel cell limits (32,767 characters)."""
    if not ids:
        return [], [], [], []
    
    # Try with MAX_MATCHES first
    limited_ids = ids[:MAX_MATCHES]
    limited_seqs = sequences[:MAX_MATCHES]
    limited_lengths = lengths[:MAX_MATCHES]
    limited_types = types[:MAX_MATCHES]
    
    # Check if it fits within cell limit
    id_str = ';'.join(limited_ids)
    seq_str = ';'.join(limited_seqs)
    length_str = ';'.join(map(str, limited_lengths))
    type_str = ';'.join(limited_types)
    
    # Find the longest string
    max_len = max(len(id_str), len(seq_str), len(length_str), len(type_str))
    
    if max_len <= EXCEL_CELL_LIMIT:
        logger.info(f"Using {len(limited_ids)} matches (fits within cell limit)")
        return limited_ids, limited_seqs, limited_lengths, limited_types
    
    # If MAX_MATCHES doesn't fit, reduce further
    logger.warning(f"MAX_MATCHES ({MAX_MATCHES}) exceeds cell limit, reducing...")
    for i in range(MAX_MATCHES - 1, 0, -1):
        limited_ids = ids[:i]
        limited_seqs = sequences[:i]
        limited_lengths = lengths[:i]
        limited_types = types[:i]
        
        id_str = ';'.join(limited_ids)
        seq_str = ';'.join(limited_seqs)
        length_str = ';'.join(map(str, limited_lengths))
        type_str = ';'.join(limited_types)
        
        max_len = max(len(id_str), len(seq_str), len(length_str), len(type_str))
        
        if max_len <= EXCEL_CELL_LIMIT:
            logger.info(f"Using {i} matches to fit within cell limit")
            return limited_ids, limited_seqs, limited_lengths, limited_types
    
    # If even 1 doesn't fit, return empty (mark for manual review)
    logger.error("Cannot fit even 1 match within cell limit, marking for manual review")
    return [], [], [], []



def main():
    parser = argparse.ArgumentParser(description='Total BLAST phosphosite mapping using full UniProtKB databases')
    parser.add_argument('--input', '-i', required=True, help='Input CSV file path')
    parser.add_argument('--species', '-s', required=True, choices=['rat', 'mouse'], help='Species name')
    parser.add_argument('--output', '-o', help='Output CSV file path (optional)')
    
    args = parser.parse_args()
    
    # Load data
    logger.info(f"Loading data from {args.input}")
    df = pd.read_csv(args.input)
    logger.info(f"Loaded {len(df)} rows")
    
    # Check if manual_review column exists
    if 'manual_review' not in df.columns:
        logger.error("Error: 'manual_review' column not found. This script should run after blast.py or grouped_pipeline.py.")
        return
    
    # Filter to only manual_review == TRUE
    manual_review_mask = df['manual_review'] == True
    rows_to_process = df[manual_review_mask]
    logger.info(f"Found {len(rows_to_process)} rows with manual_review == TRUE")
    
    if len(rows_to_process) == 0:
        logger.info("No rows to process. Exiting.")
        return
    
    # Load full UniProtKB database
    try:
        database = load_full_database(args.species)
    except Exception as e:
        logger.error(f"Failed to load full UniProtKB database: {e}")
        return
    
    # Map species to taxonomy ID filenames for BLAST database
    species_files = {
        'mouse': 'uniprotkb_taxonomy_id_10090_2025_10_06.fasta',
        'rat': 'uniprotkb_taxonomy_id_10116_2025_10_06.fasta'
    }
    
    # Create BLAST database
    db_path = f"data/blast_dbs/{args.species}_full_blast"
    fasta_path = f"data/blast_dbs/{species_files[args.species]}"
    
    if not os.path.exists(f"{db_path}.phr"):
        if not create_blast_database(fasta_path, db_path):
            logger.error("Failed to create BLAST database")
            return
    else:
        logger.info(f"Using existing BLAST database: {db_path}")
    
    # Process protein groups (only those with manual_review == TRUE)
    protein_groups = rows_to_process.groupby('Protein')
    logger.info(f"Found {len(protein_groups)} protein groups to process")
    
    for protein_id, group_df in protein_groups:
        logger.info(f"Processing protein group: {protein_id} with {len(group_df)} rows")
        
        # Get group indices
        group_indices = group_df.index.tolist()
        group_errors: List[str] = []
        
        # BLAST all sequences in the group (one per row)
        logger.info(f"Starting BLAST for all rows in group {protein_id}")
        
        # Store BLAST results for each row
        row_blast_results = []  # List of sets of protein IDs
        all_rows_successful = True
        
        for row_idx, (idx, row) in enumerate(group_df.iterrows()):
            logger.info(f"  Processing row {row_idx + 1}/{len(group_df)}")
            
            # Get phosphosite sequences for this row
            phosphosite_sequence = None
            if pd.notna(row.get('cleaned_site_motif')) and row.get('cleaned_site_motif') != '':
                motif_str = str(row['cleaned_site_motif']).strip()
                segments = [seg.strip() for seg in motif_str.split(';')] if ';' in motif_str else [motif_str]
                
                # Try first sequence, then second if first fails
                for seg_idx, seg in enumerate(segments[:2]):  # Try up to 2 sequences
                    if len(seg) < 7:
                        continue
                    
                    logger.info(f"    BLASTing sequence {seg_idx + 1}: {seg}")
                    hits = blast_sequence(seg, db_path)
                    
                    if hits:
                        logger.info(f"    Found {len(hits)} BLAST hits")
                        # Collect protein IDs from this row's BLAST
                        row_ids = set()
                        for hit in hits:
                            protein_id_hit = hit['protein_id']
                            if protein_id_hit in database:
                                protein_data = database[protein_id_hit]
                                # Check species match using scientific name
                                scientific_name = SPECIES_MAPPING.get(args.species.lower(), args.species)
                                if scientific_name.lower() in protein_data['species'].lower():
                                    row_ids.add(protein_id_hit)
                        
                        if row_ids:
                            logger.info(f"    Row {row_idx + 1}: Found {len(row_ids)} matching IDs")
                            row_blast_results.append(row_ids)
                            phosphosite_sequence = seg
                            break  # Success, don't try next sequence
                
                if phosphosite_sequence is None:
                    logger.warning(f"    Row {row_idx + 1}: No valid BLAST results")
                    all_rows_successful = False
                    break
            else:
                logger.warning(f"    Row {row_idx + 1}: No cleaned_site_motif")
                all_rows_successful = False
                break
        
        # Find intersection of all row BLAST results (IDs that appear in ALL rows)
        if not all_rows_successful or not row_blast_results:
            logger.warning(f"Not all rows in group {protein_id} had successful BLAST results")
            group_errors.append('incomplete_group_blast')
            for idx in group_indices:
                df.at[idx, 'error_reason'] = ';'.join(group_errors)
            continue
        
        # Get intersection of all sets
        common_ids = set.intersection(*row_blast_results)
        logger.info(f"Found {len(common_ids)} IDs common to all {len(row_blast_results)} rows")
        
        if not common_ids:
            logger.warning(f"No common IDs found across all rows in group {protein_id}")
            
            # Fallback: Find ID with most matches across all rows
            all_ids = set()
            for row_ids in row_blast_results:
                all_ids.update(row_ids)
            
            if all_ids:
                # Count occurrences of each ID across all rows
                id_counts = {}
                for row_ids in row_blast_results:
                    for protein_id_hit in row_ids:
                        id_counts[protein_id_hit] = id_counts.get(protein_id_hit, 0) + 1
                
                # Find ID with maximum count
                best_id = max(id_counts, key=id_counts.get)
                max_count = id_counts[best_id]
                
                logger.info(f"Fallback: Using ID {best_id} found in {max_count}/{len(row_blast_results)} rows")
                common_ids = {best_id}
                group_errors.append('no_common_blast_ids_used_fallback')
            else:
                group_errors.append('no_common_blast_ids')
                for idx in group_indices:
                    df.at[idx, 'error_reason'] = ';'.join(group_errors)
                continue
        
        # Convert to lists and get sequences
        blast_matched_ids = list(common_ids)
        blast_matched_sequences = []
        blast_matched_lengths = []
        blast_matched_types = []
        
        for protein_id_hit in blast_matched_ids:
            protein_data = database[protein_id_hit]
            blast_matched_sequences.append(protein_data['sequence'])
            blast_matched_lengths.append(protein_data['length'])
            # Map db_type to sp/tr format
            id_type = 'sp' if protein_data['db_type'] == 'swissprot' else 'tr'
            blast_matched_types.append(id_type)
        
        if blast_matched_ids:
            logger.info(f"Found {len(blast_matched_ids)} BLAST matches")
            # Apply cell limit and top 10 restriction
            limited_ids, limited_seqs, limited_lengths, limited_types = check_cell_limit(blast_matched_ids, blast_matched_sequences, blast_matched_lengths, blast_matched_types)
            
            if limited_ids:
                # Update all rows in the group
                for idx in group_indices:
                    df.at[idx, 'ID_matches'] = ';'.join(limited_ids)
                    df.at[idx, 'ID_types'] = ';'.join(limited_types)
                    df.at[idx, 'length'] = ';'.join(map(str, limited_lengths))
                    df.at[idx, 'full_sequence'] = ';'.join(limited_seqs)
                    df.at[idx, 'manual_review'] = False
                    df.at[idx, 'match_method'] = 'total_blast'
                    df.at[idx, 'error_reason'] = ''  # Clear any previous errors
            else:
                # If even limited results don't fit, keep manual review flag
                logger.warning(f"BLAST results exceed cell limit")
                group_errors.append('exceeds_excel_cell_limit_total_blast')
                for idx in group_indices:
                    df.at[idx, 'error_reason'] = ';'.join(group_errors)
        else:
            logger.warning(f"No BLAST matches found for group {protein_id}")
            group_errors.append('no_exact_total_blast_match')
            for idx in group_indices:
                df.at[idx, 'error_reason'] = ';'.join(group_errors)
    
    # Use the updated dataframe
    updated_df = df
    
    # Save results
    if args.output:
        output_file = args.output
    else:
        # Create output in same directory as input with total_blast_results.csv format
        input_path = Path(args.input)
        output_file = input_path.parent / "total_blast_results.csv"
    
    logger.info(f"Saving results to {output_file}")
    updated_df.to_csv(output_file, index=False)
    
    # Print summary
    manual_review_count = updated_df['manual_review'].sum()
    successful_count = len(updated_df) - manual_review_count
    
    logger.info(f"Processing complete:")
    logger.info(f"  Total rows: {len(updated_df)}")
    logger.info(f"  Successfully mapped: {successful_count}")
    logger.info(f"  Manual review needed: {manual_review_count}")


if __name__ == "__main__":
    main()


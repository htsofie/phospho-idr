#!/usr/bin/env python3
"""
BLAST phosphosite mapping using downloaded paper databases.

This script performs BLAST searches against species-specific UniProtKB databases
downloaded from literature papers to map phosphosites to protein sequences.

Usage:
    python paper_blast_mapper.py --input data/processed/rat/cleaned_test_data.csv --species rat
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


def load_paper_database(species: str) -> Dict[str, Dict[str, Any]]:
    """Load the paper database FASTA file and parse into dictionary."""
    db_path = f"data/blast_dbs/UniprotKB_{species}_paper.fasta"
    
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Paper database not found: {db_path}")
    
    logger.info(f"Loading paper database: {db_path}")
    
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


def check_cell_limit(ids: List[str], sequences: List[str], lengths: List[int]) -> Tuple[List[str], List[str], List[int]]:
    """Limit results to fit within Excel cell limits (32,767 characters)."""
    if not ids:
        return [], [], []
    
    # Try with MAX_MATCHES first
    limited_ids = ids[:MAX_MATCHES]
    limited_seqs = sequences[:MAX_MATCHES]
    limited_lengths = lengths[:MAX_MATCHES]
    
    # Check if it fits within cell limit
    id_str = ';'.join(limited_ids)
    seq_str = ';'.join(limited_seqs)
    length_str = ';'.join(map(str, limited_lengths))
    
    # Find the longest string
    max_len = max(len(id_str), len(seq_str), len(length_str))
    
    if max_len <= EXCEL_CELL_LIMIT:
        logger.info(f"Using {len(limited_ids)} matches (fits within cell limit)")
        return limited_ids, limited_seqs, limited_lengths
    
    # If MAX_MATCHES doesn't fit, reduce further
    logger.warning(f"MAX_MATCHES ({MAX_MATCHES}) exceeds cell limit, reducing...")
    for i in range(MAX_MATCHES - 1, 0, -1):
        limited_ids = ids[:i]
        limited_seqs = sequences[:i]
        limited_lengths = lengths[:i]
        
        id_str = ';'.join(limited_ids)
        seq_str = ';'.join(limited_seqs)
        length_str = ';'.join(map(str, limited_lengths))
        
        max_len = max(len(id_str), len(seq_str), len(length_str))
        
        if max_len <= EXCEL_CELL_LIMIT:
            logger.info(f"Using {i} matches to fit within cell limit")
            return limited_ids, limited_seqs, limited_lengths
    
    # If even 1 doesn't fit, return empty (mark for manual review)
    logger.error("Cannot fit even 1 match within cell limit, marking for manual review")
    return [], [], []



def main():
    parser = argparse.ArgumentParser(description='BLAST phosphosite mapping using paper databases')
    parser.add_argument('--input', '-i', required=True, help='Input CSV file path')
    parser.add_argument('--species', '-s', required=True, choices=['human', 'rat', 'mouse'], help='Species name')
    parser.add_argument('--output', '-o', help='Output CSV file path (optional)')
    
    args = parser.parse_args()
    
    # Load data
    logger.info(f"Loading data from {args.input}")
    df = pd.read_csv(args.input)
    logger.info(f"Loaded {len(df)} rows")
    
    # Load paper database
    try:
        database = load_paper_database(args.species)
    except Exception as e:
        logger.error(f"Failed to load paper database: {e}")
        return
    
    # Create BLAST database
    db_path = f"data/blast_dbs/{args.species}_paper_blast"
    fasta_path = f"data/blast_dbs/UniprotKB_{args.species}_paper.fasta"
    
    if not os.path.exists(f"{db_path}.phr"):
        if not create_blast_database(fasta_path, db_path):
            logger.error("Failed to create BLAST database")
            return
    else:
        logger.info(f"Using existing BLAST database: {db_path}")
    
    # Initialize new columns in the main dataframe
    df['ID_matches'] = ''
    df['length'] = ''
    df['full_sequence'] = ''
    df['manual_review'] = False
    df['match_method'] = ''
    df['error_reason'] = ''
    # (Removed temporary BLAST HSP debug columns)
    
    # Process protein groups
    protein_groups = df.groupby('Protein')
    logger.info(f"Found {len(protein_groups)} protein groups to process")
    
    for protein_id, group_df in protein_groups:
        logger.info(f"Processing protein group: {protein_id}")
        
        # Get group indices
        group_indices = group_df.index.tolist()
        group_errors: List[str] = []
        
        # BLAST on phosphosite sequences (ID search removed)
        logger.info(f"Starting BLAST for group {protein_id}")
        
        # Get all phosphosite sequences from the group
        phosphosite_sequences = []
        for idx, row in group_df.iterrows():
            if pd.notna(row.get('cleaned_site_motif')) and row.get('cleaned_site_motif') != '':
                motif_str = str(row['cleaned_site_motif']).strip()
                segments = [seg.strip() for seg in motif_str.split(';')] if ';' in motif_str else [motif_str]
                phosphosite_sequences.extend(segments)
        
        # Limit number of sequences to BLAST
        if len(phosphosite_sequences) > MAX_SEQUENCES_PER_GROUP:
            logger.warning(f"Too many sequences ({len(phosphosite_sequences)}), limiting to {MAX_SEQUENCES_PER_GROUP}")
            phosphosite_sequences = phosphosite_sequences[:MAX_SEQUENCES_PER_GROUP]
        
        if not phosphosite_sequences:
            logger.warning(f"No phosphosite sequences found for group {protein_id}")
            group_errors.append('no_phosphosite_sequences')
            for idx in group_indices:
                df.at[idx, 'manual_review'] = True
                df.at[idx, 'error_reason'] = ';'.join(group_errors)
            continue
        
        # Try BLAST on each phosphosite sequence
        blast_matched_ids = []
        blast_matched_sequences = []
        blast_matched_lengths = []
        
        for seq in phosphosite_sequences:
            if len(seq) < 7:  # Skip too short sequences
                continue
                
            logger.info(f"BLASTing sequence: {seq}")
            hits = blast_sequence(seq, db_path)
            
            if hits:
                logger.info(f"Found {len(hits)} BLAST hits for sequence {seq}")
                # Filter hits by species and get sequences from database
                # Get scientific name for species matching
                scientific_name = SPECIES_MAPPING.get(args.species.lower(), args.species)
                for hit in hits:
                    protein_id_hit = hit['protein_id']
                    logger.info(f"  Checking BLAST hit: {protein_id_hit}")
                    if protein_id_hit in database:
                        protein_data = database[protein_id_hit]
                        logger.info(f"    Found in database, species: {protein_data['species']}")
                        # Check species match using scientific name
                        if scientific_name.lower() in protein_data['species'].lower():
                            logger.info(f"    Species match! Adding {protein_id_hit}")
                            if protein_id_hit not in blast_matched_ids:
                                blast_matched_ids.append(protein_id_hit)
                                blast_matched_sequences.append(protein_data['sequence'])
                                blast_matched_lengths.append(protein_data['length'])
                        else:
                            logger.info(f"    Species mismatch: expected {scientific_name}, got {protein_data['species']}")
                    else:
                        logger.info(f"    {protein_id_hit} not found in database")
                
                if blast_matched_ids:
                    logger.info(f"Found {len(blast_matched_ids)} BLAST matches, stopping search")
                    break  # Found matches, no need to try other sequences
        
        if blast_matched_ids:
            logger.info(f"Found {len(blast_matched_ids)} BLAST matches")
            # Apply cell limit and top 10 restriction
            limited_ids, limited_seqs, limited_lengths = check_cell_limit(blast_matched_ids, blast_matched_sequences, blast_matched_lengths)
            
            if limited_ids:
                # Update all rows in the group
                for idx in group_indices:
                    df.at[idx, 'ID_matches'] = ';'.join(limited_ids)
                    df.at[idx, 'length'] = ';'.join(map(str, limited_lengths))
                    df.at[idx, 'full_sequence'] = ';'.join(limited_seqs)
                    df.at[idx, 'manual_review'] = False
                    df.at[idx, 'match_method'] = 'blast'
            else:
                # If even limited results don't fit, flag for review
                logger.warning(f"BLAST results exceed cell limit, flagging for manual review")
                group_errors.append('exceeds_excel_cell_limit_blast')
                for idx in group_indices:
                    df.at[idx, 'manual_review'] = True
                    df.at[idx, 'error_reason'] = ';'.join(group_errors)
        else:
            logger.warning(f"No BLAST matches found for group {protein_id}")
            group_errors.append('no_exact_blast_match')
            for idx in group_indices:
                df.at[idx, 'manual_review'] = True
                df.at[idx, 'error_reason'] = ';'.join(group_errors)
    
    # Use the updated dataframe
    updated_df = df
    
    # Save results
    if args.output:
        output_file = args.output
    else:
        # Create output in same directory as input with paperdb_ID.csv format
        input_path = Path(args.input)
        output_file = input_path.parent / "blast_results.csv"
    
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


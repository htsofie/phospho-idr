#!/usr/bin/env python3
"""
Simplified BLAST-based phosphosite mapping for manually flagged protein groups.

This script uses BLAST only for protein mapping (finding the correct UniProt ID and full sequence),
then relies on the existing align_to_full_seq.py script for actual alignment.

Usage:
    python blast_phosphosite_mapper_simple.py --input data.csv --species mouse --output results.csv
"""

import argparse
import os
import sys
import subprocess
import tempfile
import pandas as pd
import requests
import gzip
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Align import PairwiseAligner

# Add scripts directory to path
sys.path.append(os.path.join(os.path.dirname(__file__)))
from get_full_seq import get_uniprot_sequence
import requests

# Database URLs
SWISSPROT_URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/release/complete/complete.protein.faa.gz"

# Set up logging
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def _is_species_match(header: str, species: str) -> bool:
    """Check if a FASTA header matches the target species."""
    header_lower = header.lower()
    species_lower = species.lower()
    
    # Define species-specific patterns
    species_patterns = {
        'rat': ['rattus norvegicus', 'rattus rattus', 'rattus', 'rn_'],
        'mouse': ['mus musculus', 'mus sp.', 'mm_', 'mouse'],
        'human': ['homo sapiens', 'hs_', 'human']
    }
    
    if species_lower in species_patterns:
        patterns = species_patterns[species_lower]
        return any(pattern in header_lower for pattern in patterns)
    
    # Fallback to simple substring match
    return species_lower in header_lower


def download_database(url: str, output_path: str, species: str) -> bool:
    """Download and filter UniProt database by species."""
    logger.info(f"Downloading database from {url}")
    
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        species_count = 0
        total_count = 0
        
        with gzip.open(response.raw, 'rt') as f_in, open(output_path, 'w') as f_out:
            current_seq = []
            current_header = None
            
            for line in f_in:
                if line.startswith('>'):
                    # Process previous sequence
                    if current_header and current_seq:
                        total_count += 1
                        # More specific species matching
                        if _is_species_match(current_header, species):
                            f_out.write(current_header + '\n')
                            f_out.write(''.join(current_seq) + '\n')
                            species_count += 1
                    
                    # Start new sequence
                    current_header = line.strip()
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            
            # Process last sequence
            if current_header and current_seq:
                total_count += 1
                if _is_species_match(current_header, species):
                    f_out.write(current_header + '\n')
                    f_out.write(''.join(current_seq) + '\n')
                    species_count += 1
        
        logger.info(f"Downloaded {species_count} sequences for {species} (out of {total_count} total)")
        return species_count > 0
        
    except Exception as e:
        logger.error(f"Failed to download database: {e}")
        return False


def get_sequence_from_uniprot(uniprot_id: str) -> Optional[Tuple[str, str]]:
    """Get full protein sequence from UniProt, handling both SwissProt and TrEMBL IDs."""
    try:
        # First try the existing SwissProt function
        result = get_uniprot_sequence(uniprot_id)
        if result:
            return result
        
        # If that fails, try fetching directly from UniProt REST API
        logger.info(f"SwissProt fetch failed, trying UniProt REST API for {uniprot_id}")
        
        # UniProt REST API endpoint
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            # Parse FASTA response
            lines = response.text.strip().split('\n')
            if len(lines) >= 2:
                header = lines[0]
                sequence = ''.join(lines[1:])
                
                # Determine if it's SwissProt or TrEMBL based on header
                if 'sp|' in header:
                    seq_type = "swissprot"
                elif 'tr|' in header:
                    seq_type = "trembl"
                else:
                    seq_type = "uniprot"
                
                logger.info(f"Successfully fetched sequence via REST API (length: {len(sequence)}, type: {seq_type})")
                return sequence, seq_type
        else:
            logger.warning(f"UniProt REST API failed for {uniprot_id}: HTTP {response.status_code}")
            
    except Exception as e:
        logger.error(f"Failed to fetch sequence for {uniprot_id}: {e}")
    
    return None


def check_database_content(fasta_path: str, species: str) -> None:
    """Check what's in the database file for debugging."""
    try:
        logger.info(f"Checking database content: {fasta_path}")
        with open(fasta_path, 'r') as f:
            lines = f.readlines()
        
        # Count sequences and check headers
        sequence_count = 0
        species_matches = 0
        sample_headers = []
        
        for line in lines:
            if line.startswith('>'):
                sequence_count += 1
                header = line.strip()
                if species.lower() in header.lower():
                    species_matches += 1
                if len(sample_headers) < 5:  # Keep first 5 headers as samples
                    sample_headers.append(header)
        
        logger.info(f"Database contains {sequence_count} sequences")
        logger.info(f"Sequences with '{species}' in header: {species_matches}")
        logger.info("Sample headers:")
        for header in sample_headers:
            logger.info(f"  {header}")
            
    except Exception as e:
        logger.error(f"Failed to check database content: {e}")


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


def blast_sequence(sequence: str, db_path: str, evalue: float = 10.0) -> List[Dict[str, Any]]:
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
            '-comp_based_stats', '0', '-max_target_seqs', '50'  # Increased from 10 to 50
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
                    logger.debug(f"BLAST record query length: {blast_record.query_length}")
                    for alignment in blast_record.alignments:
                        logger.debug(f"Alignment title: {alignment.title}")
                        for hsp in alignment.hsps:
                            # Extract UniProt ID from title (handle both SwissProt and TrEMBL)
                            title = alignment.title
                            if 'sp|' in title:
                                # SwissProt format: sp|UNIPROT_ID|PROTEIN_NAME
                                parts = title.split('sp|')[1].split('|')
                                uniprot_id = parts[0]
                            elif 'tr|' in title:
                                # TrEMBL format: tr|UNIPROT_ID|PROTEIN_NAME
                                parts = title.split('tr|')[1].split('|')
                                uniprot_id = parts[0]
                            elif '|' in title:
                                # Generic format with pipes
                                parts = title.split('|')
                                if len(parts) >= 2:
                                    uniprot_id = parts[1]
                                else:
                                    uniprot_id = parts[0]
                            else:
                                # Fallback to first word
                                uniprot_id = title.split()[0]
                        
                            hit = {
                                'protein_id': uniprot_id,
                                'query_start': hsp.query_start,
                                'query_end': hsp.query_end,
                                'subject_start': hsp.sbjct_start,
                                'subject_end': hsp.sbjct_end,
                                'query_seq': hsp.query,
                                'subject_seq': hsp.sbjct,
                                'match_seq': hsp.match,
                                'identity_percent': (hsp.identities / hsp.align_length) * 100,
                                'coverage_percent': (hsp.align_length / blast_record.query_length) * 100,
                                'evalue': hsp.expect,
                                'bit_score': hsp.bits
                            }
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


def process_flagged_groups(df: pd.DataFrame, species: str, db_type: str = "swissprot") -> pd.DataFrame:
    """Process groups flagged for manual review using BLAST for mapping only."""
    logger.info(f"Processing flagged groups for species: {species}")
    
    # Check if manual_review_flag column exists
    if 'manual_review_flag' not in df.columns:
        logger.error("Error: 'manual_review_flag' column not found in dataset. This column is required.")
        raise ValueError("manual_review_flag column is required but not found in dataset")
    
    # Filter groups with manual_review_flag = True
    flagged_groups = df[df['manual_review_flag'] == True]['Protein'].unique()
    logger.info(f"Found {len(flagged_groups)} flagged groups")
    
    if len(flagged_groups) == 0:
        logger.warning("No groups flagged for manual review")
        return df
    
    # Process each flagged group
    successful = 0
    failed = 0
    
    for protein_id in flagged_groups:
        logger.info(f"Processing group: {protein_id}")
        group_idx = df[df['Protein'] == protein_id].index
        group_df = df.loc[group_idx]
        
        # Collect all phosphosite sequences
        phosphosite_sequences = []
        for idx, row in group_df.iterrows():
            if not pd.isna(row.get('cleaned_site_motif')):
                motif_str = str(row['cleaned_site_motif']).strip()
                segments = [seg.strip() for seg in motif_str.split(';')] if ';' in motif_str else [motif_str]
                phosphosite_sequences.extend(segments)
        
        if not phosphosite_sequences:
            logger.warning(f"No phosphosite sequences found for group {protein_id}")
            failed += 1
            continue
        
        # Try SwissProt first, then RefSeq if needed
        mapping_success = False
        for current_db_type in ["swissprot", "refseq"]:
            logger.info(f"Trying {current_db_type.upper()} database for group {protein_id}")
            
            # Download and prepare database
            db_name = f"{current_db_type}_{species}"
            db_path = f"data/blast_dbs/{db_name}"
            fasta_path = f"data/blast_dbs/{db_name}.fasta"
            
            os.makedirs("data/blast_dbs", exist_ok=True)
            
            # Check if BLAST database already exists (look for .phr file)
            if not os.path.exists(f"{db_path}.phr"):
                # Download database if FASTA file doesn't exist
                if not os.path.exists(fasta_path):
                    url = SWISSPROT_URL if current_db_type == "swissprot" else REFSEQ_URL
                    if not download_database(url, fasta_path, species):
                        logger.error(f"Failed to download {current_db_type} database")
                        continue
                
                # Create BLAST database
                if not create_blast_database(fasta_path, db_path):
                    logger.error(f"Failed to create BLAST database")
                    continue
            else:
                logger.info(f"Using existing BLAST database: {db_path}")
            
            # Check database content for debugging
            check_database_content(fasta_path, species)
            
            # BLAST each sequence and collect protein scores
            protein_scores = {}
            sequences_with_hits = 0
            
            logger.info(f"BLASTing {len(phosphosite_sequences)} sequences against {current_db_type.upper()} database")
            
            for seq in phosphosite_sequences:
                if len(seq) < 7:  # Skip too short sequences
                    logger.warning(f"Skipping sequence '{seq}' - too short (< 7 characters)")
                    continue
                    
                hits = blast_sequence(seq, db_path)
                logger.info(f"Sequence '{seq}': found {len(hits)} hits")
                
                if hits:
                    sequences_with_hits += 1
                    for hit in hits:
                        protein_id_hit = hit['protein_id']
                        logger.info(f"  Hit: {protein_id_hit} (identity: {hit['identity_percent']:.1f}%, evalue: {hit['evalue']:.2e})")
                        if protein_id_hit not in protein_scores:
                            protein_scores[protein_id_hit] = {'hits': [], 'total_identity': 0, 'hit_count': 0, 'sequences_matched': set()}
                        protein_scores[protein_id_hit]['hits'].append(hit)
                        protein_scores[protein_id_hit]['total_identity'] += hit['identity_percent']
                        protein_scores[protein_id_hit]['hit_count'] += 1
                        protein_scores[protein_id_hit]['sequences_matched'].add(seq)
            
            if not protein_scores:
                logger.warning(f"No BLAST hits found for group {protein_id} in {current_db_type}")
                continue
            
            # Find proteins that have high identity for ALL phosphosite sequences
            valid_proteins = []
            logger.info(f"Evaluating {len(protein_scores)} proteins for mapping")
            
            for protein_id_hit, scores in protein_scores.items():
                # Check if this protein has high identity (>= 95%) for all sequences
                all_high_identity = all(hit['identity_percent'] >= 95.0 for hit in scores['hits'])
                # Check if this protein matched all phosphosite sequences
                matched_all_sequences = len(scores['sequences_matched']) == len(phosphosite_sequences)
                
                avg_identity = scores['total_identity'] / scores['hit_count']
                logger.info(f"  Protein {protein_id_hit}: avg_identity={avg_identity:.1f}%, all_high_identity={all_high_identity}, matched_all={matched_all_sequences}")
                
                if all_high_identity and matched_all_sequences:
                    valid_proteins.append((protein_id_hit, scores))
                    logger.info(f"    -> VALID CANDIDATE")
            
            if not valid_proteins:
                logger.warning(f"No protein found with >=95% identity for all phosphosites in group {protein_id} in {current_db_type}")
                continue
            
            # If multiple proteins meet criteria, pick the one with highest average identity
            best_protein, best_scores = max(valid_proteins, 
                                          key=lambda x: x[1]['total_identity'] / x[1]['hit_count'])
            best_identity = best_scores['total_identity'] / best_scores['hit_count']
            
            # Get full sequence for best protein
            logger.info(f"Fetching sequence for UniProt ID: {best_protein}")
            seq_result = get_sequence_from_uniprot(best_protein)
            if not seq_result:
                logger.warning(f"Failed to fetch sequence for {best_protein}")
                continue
                
            protein_sequence, seq_type = seq_result
            logger.info(f"Retrieved sequence (length: {len(protein_sequence)}, type: {seq_type})")
            
            # Update the dataframe with BLAST mapping results only
            df.loc[group_idx, 'uniprot_mapped_id'] = best_protein
            df.loc[group_idx, 'mapping_source'] = f'BLAST_{current_db_type.upper()}'
            df.loc[group_idx, 'full_sequence'] = protein_sequence
            df.loc[group_idx, 'sequence_length'] = len(protein_sequence)
            df.loc[group_idx, 'sequence_type'] = 'BLAST_MAPPED'
            # Note: Do NOT clear manual_review_flag here - let alignment script handle it
            df.loc[group_idx, 'alignment_error'] = ''  # Clear alignment errors
            df.loc[group_idx, 'blast_identity'] = best_identity  # Add BLAST identity column
            df.loc[group_idx, 'blast_sequences_matched'] = len(best_scores['sequences_matched'])  # Add count of matched sequences
            
            successful += 1
            mapping_success = True
            logger.info(f"Successfully mapped group {protein_id} to {best_protein} (identity: {best_identity:.1f}%) using {current_db_type.upper()}")
            break  # Exit the database loop if successful
        
        if not mapping_success:
            failed += 1
            logger.warning(f"Failed to map group {protein_id} with 90%+ identity in both SwissProt and RefSeq")
    
    logger.info(f"BLAST processing complete: {successful} successful, {failed} failed")
    logger.info("Note: Use align_to_full_seq.py script to perform actual alignment on the mapped sequences")
    return df


def main():
    parser = argparse.ArgumentParser(description='BLAST phosphosite mapping for manually flagged groups')
    parser.add_argument('--input', '-i', required=True, help='Input CSV/Parquet file from remapping script')
    parser.add_argument('--species', '-s', required=True, choices=['mouse', 'rat', 'human'], help='Species name')
    parser.add_argument('--output', '-o', help='Output CSV file (optional - if not provided, creates blast_results.csv in same directory)')
    parser.add_argument('--db-type', choices=['swissprot', 'trembl'], default='swissprot', help='Database type to use')
    
    args = parser.parse_args()
    
    # Load data
    logger.info(f"Loading data from {args.input}")
    if args.input.endswith('.parquet'):
        df = pd.read_parquet(args.input)
    else:
        df = pd.read_csv(args.input)
    
    logger.info(f"Loaded {len(df)} rows")
    
    # Process flagged groups and update dataframe
    updated_df = process_flagged_groups(df, args.species, args.db_type)
    
    # Save updated dataframe
    if args.output:
        output_file = args.output
    else:
        # Create blast_results.csv in the same directory as input file
        input_path = Path(args.input)
        output_file = input_path.parent / "blast_results.csv"
    
    logger.info(f"Saving updated dataset to {output_file}")
    updated_df.to_csv(output_file, index=False)
    
    # Print summary
    flagged_groups = updated_df[updated_df['manual_review_flag'] == True]['Protein'].nunique()
    successful_groups = updated_df[updated_df['mapping_source'].str.contains('BLAST', na=False)]['Protein'].nunique()
    
    logger.info(f"Processing complete:")
    logger.info(f"  Groups still flagged for manual review: {flagged_groups}")
    logger.info(f"  Groups successfully mapped with BLAST: {successful_groups}")
    logger.info(f"  Next step: Run align_to_full_seq.py on the output file to perform alignment")


if __name__ == "__main__":
    main()

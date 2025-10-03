#!/usr/bin/env python3
"""
BLAST utilities for phosphosite sequence mapping.

This module provides functions to use NCBI BLAST+ for mapping phosphosite sequences
when standard alignment methods fail.
"""

import os
import tempfile
import subprocess
import pandas as pd
from typing import Dict, List, Optional, Tuple, Any
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def create_blast_database(protein_sequences: Dict[str, str], db_path: str) -> bool:
    """Create a BLAST database from protein sequences.
    
    Args:
        protein_sequences: Dictionary mapping protein IDs to sequences
        db_path: Path where to create the BLAST database
        
    Returns:
        True if successful, False otherwise
    """
    try:
        # Create FASTA file
        fasta_path = f"{db_path}.fasta"
        with open(fasta_path, 'w') as f:
            for protein_id, sequence in protein_sequences.items():
                f.write(f">{protein_id}\n{sequence}\n")
        
        # Create BLAST database
        cmd = [
            'makeblastdb',
            '-in', fasta_path,
            '-dbtype', 'prot',
            '-out', db_path,
            '-title', 'Phosphosite Protein Database'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"BLAST database created successfully at {db_path}")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"Error creating BLAST database: {e}")
        print(f"stderr: {e.stderr}")
        return False
    except Exception as e:
        print(f"Unexpected error creating BLAST database: {e}")
        return False


def blast_peptide_against_database(peptide_seq: str, db_path: str, 
                                 evalue: float = 1e-5, max_target_seqs: int = 10) -> List[Dict[str, Any]]:
    """BLAST a peptide sequence against a protein database.
    
    Args:
        peptide_seq: Peptide sequence to search
        db_path: Path to BLAST database
        evalue: E-value threshold
        max_target_seqs: Maximum number of target sequences to return
        
    Returns:
        List of BLAST hits with alignment information
    """
    try:
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_file:
            temp_file.write(f">query\n{peptide_seq}\n")
            temp_input = temp_file.name
        
        # Create temporary output file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as temp_file:
            temp_output = temp_file.name
        
        # Run BLAST
        cmd = [
            'blastp',
            '-query', temp_input,
            '-db', db_path,
            '-out', temp_output,
            '-outfmt', '5',  # XML output
            '-evalue', str(evalue),
            '-max_target_seqs', str(max_target_seqs),
            '-word_size', '2',  # Short peptides
            '-gapopen', '9',
            '-gapextend', '1'
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        # Parse XML output
        hits = _parse_blast_xml(temp_output)
        
        # Clean up temporary files
        os.unlink(temp_input)
        os.unlink(temp_output)
        
        return hits
        
    except subprocess.CalledProcessError as e:
        print(f"Error running BLAST: {e}")
        print(f"stderr: {e.stderr}")
        return []
    except Exception as e:
        print(f"Unexpected error running BLAST: {e}")
        return []
    finally:
        # Ensure cleanup
        for temp_file in [temp_input, temp_output]:
            if 'temp_file' in locals() and os.path.exists(temp_file):
                try:
                    os.unlink(temp_file)
                except:
                    pass


def _parse_blast_xml(xml_path: str) -> List[Dict[str, Any]]:
    """Parse BLAST XML output to extract hit information."""
    try:
        from Bio.Blast import NCBIXML
        
        hits = []
        with open(xml_path, 'r') as handle:
            blast_records = NCBIXML.parse(handle)
            
            for blast_record in blast_records:
                for alignment in blast_record.alignments:
                    for hsp in alignment.hsps:
                        hit = {
                            'subject_id': alignment.title.split()[0],
                            'subject_length': alignment.length,
                            'query_length': blast_record.query_length,
                            'evalue': hsp.expect,
                            'bit_score': hsp.bits,
                            'identity': hsp.identities,
                            'positive': hsp.positives,
                            'gaps': hsp.gaps,
                            'query_start': hsp.query_start,
                            'query_end': hsp.query_end,
                            'subject_start': hsp.sbjct_start,
                            'subject_end': hsp.sbjct_end,
                            'query_seq': hsp.query,
                            'subject_seq': hsp.sbjct,
                            'match_seq': hsp.match,
                            'identity_percent': (hsp.identities / hsp.align_length) * 100,
                            'coverage_percent': (hsp.align_length / blast_record.query_length) * 100
                        }
                        hits.append(hit)
        
        return hits
        
    except Exception as e:
        print(f"Error parsing BLAST XML: {e}")
        return []


def find_phosphosite_in_blast_hit(peptide_seq: str, phospho_pos_in_peptide: int, 
                                 hit: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """Find phosphosite position in BLAST hit alignment.
    
    Args:
        peptide_seq: Original peptide sequence
        phospho_pos_in_peptide: Position of phosphosite in peptide (1-based)
        hit: BLAST hit information
        
    Returns:
        Dictionary with phosphosite mapping information or None
    """
    try:
        query_seq = hit['query_seq'].replace('-', '')
        subject_seq = hit['subject_seq'].replace('-', '')
        match_seq = hit['match_seq']
        
        # Find phosphosite in alignment
        query_pos = 0
        subject_pos = 0
        phosphosite_found = False
        phosphosite_subject_pos = None
        
        for i, (q_char, s_char, m_char) in enumerate(zip(hit['query_seq'], hit['subject_seq'], match_seq)):
            if q_char != '-':
                query_pos += 1
            if s_char != '-':
                subject_pos += 1
                
            if query_pos == phospho_pos_in_peptide and q_char != '-':
                phosphosite_found = True
                phosphosite_subject_pos = subject_pos
                break
        
        if not phosphosite_found:
            return None
        
        # Calculate actual position in subject sequence
        actual_subject_pos = hit['subject_start'] + phosphosite_subject_pos - 1
        
        return {
            'blast_hit': hit,
            'phosphosite_subject_position': actual_subject_pos,
            'phosphosite_conserved': m_char == '|' if phosphosite_found else False,
            'alignment_quality': hit['identity_percent'],
            'coverage': hit['coverage_percent']
        }
        
    except Exception as e:
        print(f"Error finding phosphosite in BLAST hit: {e}")
        return None


def blast_map_phosphosite(peptide_seq: str, phospho_pos_in_peptide: int, 
                         protein_sequences: Dict[str, str], 
                         db_path: Optional[str] = None,
                         evalue: float = 1e-5) -> Optional[Dict[str, Any]]:
    """Map a phosphosite using BLAST against protein sequences.
    
    Args:
        peptide_seq: Peptide sequence containing phosphosite
        phospho_pos_in_peptide: Position of phosphosite in peptide (1-based)
        protein_sequences: Dictionary mapping protein IDs to sequences
        db_path: Path to BLAST database (if None, creates temporary one)
        evalue: E-value threshold for BLAST
        
    Returns:
        Best BLAST mapping result or None
    """
    try:
        # Create temporary database if not provided
        temp_db = None
        if db_path is None:
            temp_db = tempfile.mktemp(prefix='blast_db_')
            if not create_blast_database(protein_sequences, temp_db):
                return None
            db_path = temp_db
        
        # Run BLAST
        hits = blast_peptide_against_database(peptide_seq, db_path, evalue)
        
        if not hits:
            return None
        
        # Find best hit with phosphosite mapping
        best_result = None
        best_score = -1
        
        for hit in hits:
            phosphosite_info = find_phosphosite_in_blast_hit(
                peptide_seq, phospho_pos_in_peptide, hit
            )
            
            if phosphosite_info is None:
                continue
            
            # Score based on identity, coverage, and phosphosite conservation
            score = (
                hit['identity_percent'] * 0.4 +
                hit['coverage_percent'] * 0.3 +
                (100 if phosphosite_info['phosphosite_conserved'] else 0) * 0.3
            )
            
            if score > best_score:
                best_score = score
                best_result = {
                    'protein_id': hit['subject_id'],
                    'aligned_position': phosphosite_info['phosphosite_subject_position'],
                    'identity': hit['identity_percent'],
                    'coverage': hit['coverage_percent'],
                    'evalue': hit['evalue'],
                    'bit_score': hit['bit_score'],
                    'phosphosite_conserved': phosphosite_info['phosphosite_conserved'],
                    'blast_hit': hit,
                    'mapping_method': 'blast'
                }
        
        # Clean up temporary database
        if temp_db:
            for ext in ['.fasta', '.phr', '.pin', '.psq']:
                temp_file = temp_db + ext
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
        
        return best_result
        
    except Exception as e:
        print(f"Error in BLAST phosphosite mapping: {e}")
        return None


def integrate_blast_into_remapping(group_df: pd.DataFrame, species: str, 
                                 min_identity: float = 70.0) -> Optional[Dict[str, Any]]:
    """Integrate BLAST mapping into the remapping pipeline.
    
    This function should be called when standard alignment methods fail.
    
    Args:
        group_df: DataFrame containing the protein group
        species: Species name
        min_identity: Minimum identity threshold
        
    Returns:
        BLAST mapping result or None
    """
    try:
        # Extract protein sequences from the group
        protein_sequences = {}
        for idx, row in group_df.iterrows():
            if not pd.isna(row.get('full_sequence')) and row.get('full_sequence'):
                protein_id = f"protein_{idx}"
                protein_sequences[protein_id] = str(row['full_sequence'])
        
        if not protein_sequences:
            return None
        
        # Try BLAST mapping for each row in the group
        blast_results = {}
        all_successful = True
        
        for idx, row in group_df.iterrows():
            if (pd.isna(row.get('cleaned_site_motif')) or 
                pd.isna(row.get('motif_position')) or 
                not row.get('cleaned_site_motif') or 
                not row.get('motif_position')):
                all_successful = False
                break
            
            # Get peptide and position
            motif_str = str(row['cleaned_site_motif']).strip()
            pos_str = str(row['motif_position']).strip()
            peptide_segments = [seg.strip() for seg in motif_str.split(';')] if ';' in motif_str else [motif_str]
            position_segments = [p.strip() for p in pos_str.split(';')] if ';' in pos_str else [pos_str]
            
            # Try each peptide segment
            best_blast_result = None
            for pep, pos in zip(peptide_segments, position_segments):
                if not pep or not pos:
                    continue
                try:
                    phospho_pos = int(pos)
                    blast_result = blast_map_phosphosite(pep, phospho_pos, protein_sequences)
                    if blast_result and blast_result['identity'] >= min_identity:
                        if (best_blast_result is None or 
                            blast_result['identity'] > best_blast_result['identity']):
                            best_blast_result = blast_result
                except ValueError:
                    continue
            
            if best_blast_result is None:
                all_successful = False
                break
            
            blast_results[idx] = best_blast_result
        
        if not all_successful or not blast_results:
            return None
        
        # Return successful BLAST mapping
        return {
            'mapping_method': 'blast',
            'blast_results': blast_results,
            'success': True
        }
        
    except Exception as e:
        print(f"Error in BLAST integration: {e}")
        return None

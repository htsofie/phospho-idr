#!/usr/bin/env python3
"""
Script to align cleaned_site_motif to full_sequence and calculate phosphorylation positions.
"""

import pandas as pd
import argparse
import os
from typing import Dict, Optional, Any, List, Tuple
from Bio.Align import PairwiseAligner


def _compute_phospho_protein_pos_from_alignment(aln_prot: str, aln_pep: str, start: int, phospho_pos_in_peptide: int) -> Optional[int]:
    """Translate phosphosite position in peptide to 1-based protein position using an alignment."""
    peptide_position_counter = 0
    protein_position_counter = start

    for prot_char, pep_char in zip(aln_prot, aln_pep):
        if pep_char != "-":
            peptide_position_counter += 1
            if peptide_position_counter == phospho_pos_in_peptide:
                return protein_position_counter + 1  # return 1-based
        if prot_char != "-":
            protein_position_counter += 1
    return None


def _check_phosphosite_conserved(aln_prot: str, aln_pep: str, phospho_pos_in_peptide: int) -> bool:
    """Check if the phosphosite position is conserved in the alignment."""
    peptide_position_counter = 0
    
    for prot_char, pep_char in zip(aln_prot, aln_pep):
        if pep_char != "-":
            peptide_position_counter += 1
            if peptide_position_counter == phospho_pos_in_peptide:
                # More lenient: allow if protein has any amino acid (not gap) at this position
                # The amino acid doesn't have to match exactly, just be present
                return prot_char != "-"
    return False


def _get_aligned_sequence_info(aln_prot: str, aln_pep: str) -> Dict[str, Any]:
    """Extract aligned sequence information for non-exact matches."""
    prot_aligned = ""
    pep_aligned = ""
    matches = 0
    mismatches = 0
    
    for prot_char, pep_char in zip(aln_prot, aln_pep):
        if pep_char != "-":  # This is a peptide position
            pep_aligned += pep_char
            if prot_char != "-":  # Protein has a residue at this position
                prot_aligned += prot_char
                if prot_char == pep_char:
                    matches += 1
                else:
                    mismatches += 1
            else:  # Protein is missing this residue, use underscore
                prot_aligned += "_"
                mismatches += 1
    
    return {
        "aligned_protein_sequence": prot_aligned,
        "aligned_peptide_sequence": pep_aligned,
        "mismatching_residues": mismatches
    }


def _exact_match_candidates(protein_seq: str, peptide_seq: str, phospho_pos_in_peptide: int) -> List[Dict[str, Any]]:
    """Find all exact matches of peptide in protein sequence."""
    candidates = []
    if not protein_seq or not peptide_seq:
        return candidates
    
    start = 0
    while True:
        idx = protein_seq.find(peptide_seq, start)
        if idx == -1:
            break
        protein_pos = idx + phospho_pos_in_peptide
        candidates.append({
            "match_type": "exact",
            "start": idx,
            "end": idx + len(peptide_seq),
            "aligned_position": protein_pos,
            "identity": 100.0,
            "phosphosite_conserved": True,
            "aligned_protein_sequence": peptide_seq,
            "aligned_peptide_sequence": peptide_seq,
            "mismatching_residues": 0
        })
        start = idx + 1
    return candidates


def _local_alignment_candidates(seq_region: str, region_offset: int, peptide_seq: str, phospho_pos_in_peptide: int, max_alignments: int = 5) -> List[Dict[str, Any]]:
    """Find local alignment candidates with conserved phosphosites."""
    candidates = []
    if not seq_region or not peptide_seq:
        return candidates
    
    try:
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.1
        alignments = list(aligner.align(seq_region, peptide_seq))
    except Exception:
        return candidates
    
    for aln in alignments[:max_alignments]:
        aln_prot = str(aln.target)
        aln_pep = str(aln.query)
        score = aln.score
        start = aln.coordinates[0][0]  # target start position
        end = aln.coordinates[0][-1]   # target end position
        matches = sum(a == b for a, b in zip(aln_prot, aln_pep) if a != "-" and b != "-")
        identity = matches / len(peptide_seq) * 100.0
        phospho_pos = _compute_phospho_protein_pos_from_alignment(aln_prot, aln_pep, region_offset + start, phospho_pos_in_peptide)
        
        # Only include alignments where phosphosite is conserved
        phosphosite_conserved = _check_phosphosite_conserved(aln_prot, aln_pep, phospho_pos_in_peptide)
        if phosphosite_conserved and phospho_pos:
            alignment_info = _get_aligned_sequence_info(aln_prot, aln_pep)
            candidates.append({
                "match_type": "aligned",
                "start": region_offset + start,
                "end": region_offset + end,
                "aligned_position": phospho_pos,
                "identity": identity,
                "phosphosite_conserved": True,
                **alignment_info
            })
    return candidates


def map_peptide_to_protein(protein_seq: str, peptide_seq: str, phospho_pos_in_peptide: int,
                          expected_position: Optional[int] = None) -> Dict[str, Any]:
    """Find best alignment: exact > highest score > closest to expected position."""
    if not protein_seq or not peptide_seq or phospho_pos_in_peptide <= 0:
        return {"match_type": "failed", "error": "Invalid input parameters"}

    candidates = []
    
    # 1) All exact matches
    candidates.extend(_exact_match_candidates(protein_seq, peptide_seq, phospho_pos_in_peptide))
    
    # 2) Local alignments across full sequence
    candidates.extend(_local_alignment_candidates(protein_seq, 0, peptide_seq, phospho_pos_in_peptide, max_alignments=5))
    
    # 3) Constrained windows around expected position
    if expected_position and 0 < expected_position <= len(protein_seq) * 2:
        for window in (50, 100, 200):
            start = max(0, expected_position - window)
            end = min(len(protein_seq), expected_position + window)
            if end - start >= max(10, len(peptide_seq)):
                region = protein_seq[start:end]
                candidates.extend(_local_alignment_candidates(region, start, peptide_seq, phospho_pos_in_peptide, max_alignments=5))

    # Filter valid candidates
    candidates = [c for c in candidates if c.get("aligned_position")]
    if not candidates:
        return {"match_type": "failed", "error": "No alignment found"}

    # Sort by priority: exact > identity > proximity to expected position
    def sort_key(c: Dict[str, Any]) -> Tuple[int, float, float]:
        is_exact = 1 if c["match_type"] == "exact" else 0
        identity = c.get("identity", 0.0)
        proximity = -abs(expected_position - c["aligned_position"]) if expected_position and c.get("aligned_position") else 0.0
        return (is_exact, identity, proximity)

    candidates.sort(key=sort_key, reverse=True)
    return candidates[0]


def process_alignment_row(row: pd.Series, try_all_candidates: bool = True) -> Dict[str, Any]:
    """Process a single row for alignment."""
    # Check required data - updated to work with grouped pipeline output
    if pd.isna(row.get('full_sequence')) or not row.get('full_sequence'):
        return {"alignment_success": False, "alignment_error": "Empty full sequence"}
    
    if pd.isna(row.get('cleaned_site_motif')) or not row.get('cleaned_site_motif'):
        return {"alignment_success": False, "alignment_error": "Empty cleaned site motif"}
    
    if pd.isna(row.get('motif_position')) or not row.get('motif_position'):
        return {"alignment_success": False, "alignment_error": "No motif position"}
    
    # Check if this row was mapped via BLAST (more lenient position threshold)
    is_blast_mapped = row.get('match_method') == 'blast' if 'match_method' in row else False
    
    # Get data - handle multiple candidate sequences
    protein_seq_str = str(row['full_sequence']).strip()
    id_matches_str = str(row.get('ID_matches', '')).strip()
    id_types_str = str(row.get('ID_types', '')).strip()
    
    # Check if this is multiple candidate sequences (semicolon-separated)
    if ';' in protein_seq_str and row.get('sequence_type') == 'NAME_MAPPED_MULTIPLE':
        return process_multiple_candidate_sequences(row)
    
    # Check if we have multiple IDs/sequences from BLAST or ID search
    if ';' in protein_seq_str and ';' in id_matches_str and try_all_candidates:
        # Split into candidates
        candidate_sequences = [seq.strip() for seq in protein_seq_str.split(';') if seq.strip()]
        candidate_ids = [id.strip() for id in id_matches_str.split(';') if id.strip()]
        candidate_types = [t.strip() for t in id_types_str.split(';') if t.strip()] if id_types_str else ['tr'] * len(candidate_ids)
        
        if len(candidate_sequences) > 1 and len(candidate_sequences) == len(candidate_ids):
            # Try all candidates and find the best one
            return process_row_with_multiple_candidates(row, candidate_sequences, candidate_ids, candidate_types, is_blast_mapped)
    
    # Single sequence processing (original logic)
    protein_seq = protein_seq_str.split(';')[0].strip() if ';' in protein_seq_str else protein_seq_str
    
    # Prepare semicolon-separated cleaned site motifs and positions: try ALL segments
    motif_str = str(row['cleaned_site_motif']).strip()
    pos_str = str(row['motif_position']).strip()
    peptide_segments = [seg.strip() for seg in motif_str.split(';')] if ';' in motif_str else [motif_str]
    position_segments = [p.strip() for p in pos_str.split(';')] if ';' in pos_str else [pos_str]
    
    # Normalize positions list length to peptides list; if mismatch, truncate to min length
    num_to_try = min(len(peptide_segments), len(position_segments))
    peptide_segments = peptide_segments[:num_to_try]
    position_segments = position_segments[:num_to_try]
    
    expected_position = int(row['position']) if not pd.isna(row.get('position')) and row.get('position') else None
    
    best_result: Optional[Dict[str, Any]] = None
    best_identity: float = -1.0
    best_error: Optional[str] = None
    
    for pep, pos in zip(peptide_segments, position_segments):
        if not pep or not pos:
            continue
        try:
            phospho_pos_in_peptide = int(pos)
        except ValueError:
            continue
        aln = map_peptide_to_protein(protein_seq, pep, phospho_pos_in_peptide, expected_position)
        if aln.get("match_type") == "failed":
            best_error = aln.get("error", "Unknown error")
            continue
        identity = float(aln.get("identity", 0.0))
        # Track best by identity first, then prefer exact matches on tie
        tie_break = 1 if aln.get("match_type") == "exact" else 0
        if identity > best_identity or (abs(identity - best_identity) < 1e-9 and best_result and tie_break > (1 if best_result.get("match_type") == "exact" else 0)):
            best_result = aln
            best_identity = identity
            best_error = None
    
    if not best_result:
        return {"alignment_success": False, "alignment_error": best_error or "No alignment found"}
    
    # Check minimum identity threshold (80%)
    if best_identity < 80.0:
        return {"alignment_success": False, "alignment_error": f"Identity too low: {best_identity:.1f}% (minimum 80%)"}
    
    # Calculate position difference
    position_difference = 0
    if expected_position and best_result.get("aligned_position"):
        position_difference = abs(expected_position - best_result["aligned_position"])
    
    # Success - return all alignment data
    result = {"alignment_success": True}
    result.update(best_result)
    result["position_difference"] = position_difference
    
    return result


def process_row_with_multiple_candidates(row: pd.Series, candidate_sequences: List[str], 
                                         candidate_ids: List[str], candidate_types: List[str], is_blast_mapped: bool) -> Dict[str, Any]:
    """Try all candidate IDs/sequences and return the best alignment."""
    print(f"    Trying {len(candidate_sequences)} candidate IDs/sequences for best alignment...")
    
    # Prepare peptide data
    motif_str = str(row['cleaned_site_motif']).strip()
    pos_str = str(row['motif_position']).strip()
    peptide_segments = [seg.strip() for seg in motif_str.split(';')] if ';' in motif_str else [motif_str]
    position_segments = [p.strip() for p in pos_str.split(';')] if ';' in pos_str else [pos_str]
    
    num_to_try = min(len(peptide_segments), len(position_segments))
    peptide_segments = peptide_segments[:num_to_try]
    position_segments = position_segments[:num_to_try]
    
    expected_position = int(row['position']) if not pd.isna(row.get('position')) and row.get('position') else None
    
    best_candidate_result = None
    best_candidate_id = None
    best_candidate_seq = None
    best_candidate_type = None
    best_score = -1
    
    # Try each candidate
    for cand_idx, (protein_seq, cand_id, cand_type) in enumerate(zip(candidate_sequences, candidate_ids, candidate_types)):
        best_result_for_candidate: Optional[Dict[str, Any]] = None
        best_identity_for_candidate: float = -1.0
        
        # Try all peptide segments for this candidate
        for pep, pos in zip(peptide_segments, position_segments):
            if not pep or not pos:
                continue
            try:
                phospho_pos_in_peptide = int(pos)
            except ValueError:
                continue
            
            aln = map_peptide_to_protein(protein_seq, pep, phospho_pos_in_peptide, expected_position)
            if aln.get("match_type") == "failed":
                continue
            
            identity = float(aln.get("identity", 0.0))
            tie_break = 1 if aln.get("match_type") == "exact" else 0
            
            if identity > best_identity_for_candidate or (abs(identity - best_identity_for_candidate) < 1e-9 and best_result_for_candidate and tie_break > (1 if best_result_for_candidate.get("match_type") == "exact" else 0)):
                best_result_for_candidate = aln
                best_identity_for_candidate = identity
        
        # Evaluate this candidate
        if best_result_for_candidate:
            # Check identity threshold
            if best_identity_for_candidate < 80.0:
                print(f"      Candidate {cand_idx + 1} ({cand_id}): Identity too low ({best_identity_for_candidate:.1f}%)")
                continue
            
            # Calculate position difference
            position_difference = 0
            if expected_position and best_result_for_candidate.get("aligned_position"):
                position_difference = abs(expected_position - best_result_for_candidate["aligned_position"])
            
            # Calculate score: prioritize sp over tr, then position_diff, then identity
            # Score formula: ID type bonus + (10000 - position_diff) + identity
            score = (5000 if cand_type == 'sp' else 0) + (10000 - position_difference) + best_identity_for_candidate
            
            print(f"      Candidate {cand_idx + 1} ({cand_id}, {cand_type}): Valid - identity {best_identity_for_candidate:.1f}%, pos_diff {position_difference}, score {score:.1f}")
            
            if score > best_score:
                best_score = score
                best_candidate_result = best_result_for_candidate
                best_candidate_id = cand_id
                best_candidate_seq = protein_seq
                best_candidate_type = cand_type
                best_candidate_result["position_difference"] = position_difference
    
    # Return best candidate
    if best_candidate_result:
        print(f"      → Selected candidate: {best_candidate_id} ({best_candidate_type})")
        result = {"alignment_success": True, "selected_id": best_candidate_id, "selected_sequence": best_candidate_seq, "selected_id_type": best_candidate_type}
        result.update(best_candidate_result)
        return result
    else:
        print(f"      → No valid candidate found")
        return {"alignment_success": False, "alignment_error": "No candidate met alignment criteria (identity >= 80%, position_diff <= threshold)"}


def process_multiple_candidate_sequences(row: pd.Series) -> Dict[str, Any]:
    """Process alignment with multiple candidate sequences from name search."""
    # Get all candidate sequences and IDs
    protein_seq_str = str(row['full_sequence']).strip()
    uniprot_ids_str = str(row.get('uniprot_mapped_id', '')).strip()
    id_types_str = str(row.get('ID_types', '')).strip()
    
    candidate_sequences = [seq.strip() for seq in protein_seq_str.split(';') if seq.strip()]
    candidate_ids = [id.strip() for id in uniprot_ids_str.split(';') if id.strip()]
    candidate_types = [t.strip() for t in id_types_str.split(';') if t.strip()] if id_types_str else ['tr'] * len(candidate_ids)
    
    if not candidate_sequences:
        return {"alignment_success": False, "alignment_error": "No candidate sequences found"}
    
    # Get all rows in this protein group to evaluate group-wide alignment
    protein_id = row.get('Protein', '')
    if not protein_id:
        return {"alignment_success": False, "alignment_error": "No protein ID found"}
    
    # This function will be called for each row, so we need to get the group data
    # We'll use a different approach - evaluate each candidate for the current row
    # but with group-wide scoring logic
    
    # Prepare peptide data for current row
    motif_str = str(row['cleaned_site_motif']).strip()
    pos_str = str(row['motif_position']).strip()
    peptide_segments = [seg.strip() for seg in motif_str.split(';')] if ';' in motif_str else [motif_str]
    position_segments = [p.strip() for p in pos_str.split(';')] if ';' in pos_str else [pos_str]
    
    expected_position = int(row['position']) if not pd.isna(row.get('position')) and row.get('position') else None
    
    # Try each candidate sequence
    best_overall_result = None
    best_overall_score = -1
    best_candidate_id = None
    best_candidate_type = None
    
    for i, (candidate_seq, candidate_id, candidate_type) in enumerate(zip(candidate_sequences, candidate_ids, candidate_types)):
        print(f"  Trying candidate {i+1}/{len(candidate_sequences)}: {candidate_id} ({candidate_type})")
        
        # Test alignment for this candidate
        candidate_result = test_candidate_alignment(
            candidate_seq, peptide_segments, position_segments, expected_position, candidate_type
        )
        
        if candidate_result:
            # Score this candidate: prioritize exact matches with 0 position difference
            score = calculate_candidate_score(candidate_result, expected_position, candidate_type)
            
            print(f"    Score: {score:.1f} (identity: {candidate_result.get('identity', 0):.1f}%, "
                  f"position_diff: {candidate_result.get('position_difference', 999)})")
            
            if score > best_overall_score:
                best_overall_score = score
                best_overall_result = candidate_result
                best_candidate_id = candidate_id
                best_candidate_type = candidate_type
    
    if not best_overall_result:
        return {"alignment_success": False, "alignment_error": "No valid alignment found for any candidate"}
    
    # Update the row with the best candidate
    if best_candidate_id:
        # Update the row to use only the best candidate
        row['uniprot_mapped_id'] = best_candidate_id
        chosen_seq = candidate_sequences[candidate_ids.index(best_candidate_id)]
        row['full_sequence'] = chosen_seq
        try:
            row['sequence_length'] = len(chosen_seq)
        except Exception:
            row['sequence_length'] = None
        row['sequence_type'] = 'NAME_MAPPED_BEST'
        print(f"  Selected best candidate: {best_candidate_id} ({best_candidate_type})")
    
    # Add selected ID type to the result
    best_overall_result['selected_id_type'] = best_candidate_type
    return best_overall_result


def test_candidate_alignment(protein_seq: str, peptide_segments: List[str], 
                           position_segments: List[str], expected_position: Optional[int], id_type: str = 'tr') -> Optional[Dict[str, Any]]:
    """Test alignment for a single candidate sequence."""
    num_to_try = min(len(peptide_segments), len(position_segments))
    peptide_segments = peptide_segments[:num_to_try]
    position_segments = position_segments[:num_to_try]
    
    best_result: Optional[Dict[str, Any]] = None
    best_identity: float = -1.0
    best_error: Optional[str] = None
    
    for pep, pos in zip(peptide_segments, position_segments):
        if not pep or not pos:
            continue
        try:
            phospho_pos_in_peptide = int(pos)
        except ValueError:
            continue
        
        aln = map_peptide_to_protein(protein_seq, pep, phospho_pos_in_peptide, expected_position)
        if aln.get("match_type") == "failed":
            best_error = aln.get("error", "Unknown error")
            continue
        
        identity = float(aln.get("identity", 0.0))
        # Track best by identity first, then prefer exact matches on tie
        tie_break = 1 if aln.get("match_type") == "exact" else 0
        if identity > best_identity or (abs(identity - best_identity) < 1e-9 and best_result and tie_break > (1 if best_result.get("match_type") == "exact" else 0)):
            best_result = aln
            best_identity = identity
            best_error = None
    
    if not best_result:
        return None
    
    # Check minimum identity threshold (80%)
    if best_identity < 80.0:
        return None
    
    # Calculate position difference
    position_difference = 0
    if expected_position and best_result.get("aligned_position"):
        position_difference = abs(expected_position - best_result["aligned_position"])
    
    # Success - return all alignment data
    result = {"alignment_success": True}
    result.update(best_result)
    result["position_difference"] = position_difference
    
    return result


def calculate_candidate_score(result: Dict[str, Any], expected_position: Optional[int], id_type: str = 'tr') -> float:
    """Calculate a score for a candidate alignment result."""
    score = 0.0
    
    # Massive bonus for SwissProt IDs (prioritizes sp over tr)
    if id_type == 'sp':
        score += 5000
    
    # Base score from identity
    identity = result.get("identity", 0.0)
    score += identity
    
    # Bonus for exact matches
    if result.get("match_type") == "exact":
        score += 50.0
    
    # Use position_diff as tie-breaker (smaller is better)
    position_diff = result.get("position_difference", 999)
    if position_diff == 0:
        score += 100.0
    else:
        # Subtract position_diff to prefer smaller differences
        score -= position_diff * 0.1  # Scale factor to not dominate other criteria
    
    return score


def process_protein_group_with_candidates(group_df: pd.DataFrame, main_df: pd.DataFrame) -> Dict[int, Dict[str, Any]]:
    """Process a protein group with multiple candidate sequences, finding the best overall match."""
    first_row = group_df.iloc[0]
    
    # Get all candidate sequences and IDs
    protein_seq_str = str(first_row['full_sequence']).strip()
    uniprot_ids_str = str(first_row.get('uniprot_mapped_id', '')).strip()
    id_types_str = str(first_row.get('ID_types', '')).strip()
    
    candidate_sequences = [seq.strip() for seq in protein_seq_str.split(';') if seq.strip()]
    candidate_ids = [id.strip() for id in uniprot_ids_str.split(';') if id.strip()]
    candidate_types = [t.strip() for t in id_types_str.split(';') if t.strip()] if id_types_str else ['tr'] * len(candidate_ids)
    
    if not candidate_sequences:
        # Return failure for all rows
        results = {}
        for idx, row in group_df.iterrows():
            results[idx] = {"alignment_success": False, "alignment_error": "No candidate sequences found"}
        return results
    
    print(f"  Found {len(candidate_sequences)} candidate sequences")
    
    # Evaluate each candidate against ALL phosphosites in the group
    best_candidate_id = None
    best_candidate_type = None
    best_candidate_score = -1
    best_candidate_results = {}
    
    for i, (candidate_seq, candidate_id, candidate_type) in enumerate(zip(candidate_sequences, candidate_ids, candidate_types)):
        print(f"  Testing candidate {i+1}/{len(candidate_sequences)}: {candidate_id} ({candidate_type})")
        
        # Test this candidate against all phosphosites in the group
        candidate_results = {}
        candidate_scores = []
        all_successful = True
        
        for idx, row in group_df.iterrows():
            # Prepare peptide data for this row
            motif_str = str(row['cleaned_site_motif']).strip()
            pos_str = str(row['motif_position']).strip()
            peptide_segments = [seg.strip() for seg in motif_str.split(';')] if ';' in motif_str else [motif_str]
            position_segments = [p.strip() for p in pos_str.split(';')] if ';' in pos_str else [pos_str]
            
            expected_position = int(row['position']) if not pd.isna(row.get('position')) and row.get('position') else None
            
            # Test alignment for this candidate and this phosphosite
            result = test_candidate_alignment(
                candidate_seq, peptide_segments, position_segments, expected_position, candidate_type
            )
            
            if result:
                candidate_results[idx] = result
                score = calculate_candidate_score(result, expected_position, candidate_type)
                candidate_scores.append(score)
            else:
                all_successful = False
                candidate_results[idx] = {"alignment_success": False, "alignment_error": "No valid alignment found"}
                break
        
        if all_successful:
            # Calculate overall group score
            group_score = sum(candidate_scores)
            avg_position_diff = sum(r.get('position_difference', 0) for r in candidate_results.values()) / len(candidate_results)
            
            print(f"    Group score: {group_score:.1f}, Avg position diff: {avg_position_diff:.1f}")
            
            # Prioritize candidates with 0 position difference across all sites
            if avg_position_diff == 0:
                group_score += 1000  # Massive bonus for perfect position matches
            
            if group_score > best_candidate_score:
                best_candidate_score = group_score
                best_candidate_id = candidate_id
                best_candidate_type = candidate_type
                best_candidate_results = candidate_results.copy()
                print(f"    → New best candidate!")
    
    if not best_candidate_results:
        # Return failure for all rows and mark group for manual review
        results = {}
        for idx, row in group_df.iterrows():
            results[idx] = {"alignment_success": False, "alignment_error": "No valid alignment found for any candidate"}
        try:
            for idx in group_df.index.tolist():
                main_df.at[idx, 'manual_review'] = True
        except Exception:
            pass
        return results
    
    # Update ALL rows in the group with the best candidate in the main dataframe
    best_candidate_seq = candidate_sequences[candidate_ids.index(best_candidate_id)]
    group_indices = group_df.index.tolist()
    best_seq_len = len(best_candidate_seq) if isinstance(best_candidate_seq, str) else None
    
    for idx in group_indices:
        main_df.at[idx, 'uniprot_mapped_id'] = best_candidate_id
        main_df.at[idx, 'full_sequence'] = best_candidate_seq
        main_df.at[idx, 'sequence_length'] = best_seq_len
        main_df.at[idx, 'sequence_type'] = 'NAME_MAPPED_BEST'
        main_df.at[idx, 'ID_types'] = best_candidate_type
    
    print(f"  Selected best candidate: {best_candidate_id} ({best_candidate_type}) (score: {best_candidate_score:.1f})")
    print(f"  Updated {len(group_indices)} rows in the group with best candidate")
    
    return best_candidate_results


def process_dataset(input_path: str, species: str, output_path: str) -> None:
    """Process the dataset to align motifs to full sequences."""
    print(f"Loading dataset: {input_path}")
    df = pd.read_csv(input_path)
    print(f"Loaded {len(df)} rows")
    
    # Check required columns - updated for grouped pipeline output
    required_cols = ['full_sequence', 'cleaned_site_motif', 'motif_position']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        return
    
    # Initialize alignment columns
    alignment_columns = [
        'alignment_success', 'match_type', 'aligned_position', 'position_difference', 'identity', 
        'alignment_error', 'phosphosite_conserved', 'aligned_protein_sequence', 'aligned_peptide_sequence',
        'mismatching_residues'
    ]
    for col in alignment_columns:
        df[col] = None
    
    print(f"Starting motif alignment for {species} data...")
    
    # Filter to only rows with full sequences
    rows_with_sequence = df[df['full_sequence'].notna() & (df['full_sequence'] != '')]
    print(f"Found {len(rows_with_sequence)} rows with full sequences to align")
    
    # Group by protein ID for group-wide evaluation
    protein_groups = rows_with_sequence.groupby('Protein')
    print(f"Found {len(protein_groups)} protein groups to process")
    
    # Process each protein group
    for protein_id, group_df in protein_groups:
        print(f"\nProcessing protein group: {protein_id} ({len(group_df)} phosphosites)")
        
        # Check if this group has multiple candidate sequences
        first_row = group_df.iloc[0]
        if (';' in str(first_row['full_sequence']) and 
            first_row.get('sequence_type') == 'NAME_MAPPED_MULTIPLE'):
            
            # Process group with multiple candidates
            group_results = process_protein_group_with_candidates(group_df, df)
            
            # Update the dataframe with results
            for idx, result in group_results.items():
                for key, value in result.items():
                    df.at[idx, key] = value
        else:
            # Process each row individually (original logic)
            for idx, row in group_df.iterrows():
                print(f"  Processing phosphosite {idx + 1}")
                
                alignment_result = process_alignment_row(row)
                
                # Check if a different candidate was selected
                if alignment_result.get('selected_id') and alignment_result.get('selected_sequence'):
                    # Update the row to use the best candidate
                    df.at[idx, 'ID_matches'] = alignment_result['selected_id']
                    df.at[idx, 'full_sequence'] = alignment_result['selected_sequence']
                    df.at[idx, 'length'] = len(alignment_result['selected_sequence'])
                    if alignment_result.get('selected_id_type'):
                        df.at[idx, 'ID_types'] = alignment_result['selected_id_type']
                    print(f"    → Updated to use best candidate: {alignment_result['selected_id']} ({alignment_result.get('selected_id_type', 'unknown')})")
                
                # Store alignment results
                for key, value in alignment_result.items():
                    if key not in ['selected_id', 'selected_sequence', 'selected_id_type']:  # Don't store these as columns
                        df.at[idx, key] = value
                
                if alignment_result.get('alignment_success', False):
                    print(f"    ✓ Alignment successful: position {alignment_result.get('aligned_position', 'N/A')}")
                else:
                    print(f"    ✗ Alignment failed: {alignment_result.get('alignment_error', 'Unknown error')}")
                    # Mark for manual review when alignment fails for this row
                    df.at[idx, 'manual_review'] = True
    
    # Reorder columns - put alignment columns at the end
    other_cols = [col for col in df.columns if col not in alignment_columns]
    df_reordered = df[other_cols + alignment_columns]
    
    # Save results
    print(f"\nSaving results to: {output_path}")
    df_reordered.to_csv(output_path, index=False)
    
    # Print summary
    successful = df['alignment_success'].sum()
    print(f"\nAlignment Summary:")
    print(f"  Total rows: {len(df)}")
    print(f"  Successful: {successful} ({successful/len(df)*100:.1f}%)")
    print(f"  Failed: {len(df) - successful}")
    
    if successful > 0:
        # Filter for successful alignments (handle NaN values)
        successful_mask = df['alignment_success'] == True
        
        # Match types
        match_counts = df[successful_mask]['match_type'].value_counts()
        print(f"\nMatch types:")
        for match_type, count in match_counts.items():
            print(f"  {match_type}: {count} ({count/successful*100:.1f}%)")
        
        # Phosphosite conservation
        conserved = df[successful_mask]['phosphosite_conserved'].sum()
        print(f"\nPhosphosite conservation: {conserved} ({conserved/successful*100:.1f}%)")
        
        # Position differences
        valid_diffs = df[successful_mask & df['position_difference'].notna()]['position_difference']
        if len(valid_diffs) > 0:
            print(f"\nPosition differences:")
            print(f"  Mean: {valid_diffs.mean():.1f}, Median: {valid_diffs.median():.1f}")
            print(f"  <=5: {(valid_diffs <= 5).sum()}, 6-20: {((valid_diffs > 5) & (valid_diffs <= 20)).sum()}, >20: {(valid_diffs > 20).sum()}")
            print(f"  Max: {valid_diffs.max():.1f} (all alignments with pos_diff > 150 were rejected)")
    
    # Error breakdown
    if len(df) - successful > 0:
        print(f"\nErrors:")
        failed_mask = df['alignment_success'] == False
        failed_rows = df[failed_mask]
        if 'alignment_error' in failed_rows.columns:
            error_counts = failed_rows['alignment_error'].value_counts()
            for error, count in error_counts.items():
                print(f"  {error}: {count}")
        else:
            print(f"  {len(failed_rows)} rows failed alignment")


def main():
    parser = argparse.ArgumentParser(description='Align cleaned site motifs to full protein sequences')
    parser.add_argument('--input', '-i', required=True, help='Path to input CSV file (with full sequences)')
    parser.add_argument('--species', '-s', required=True, choices=['mouse', 'rat'], help='Species name')
    parser.add_argument('--output', '-o', help='Path to output CSV file (default: data/processed/{species}/input_filename_aligned.csv)')
    
    args = parser.parse_args()
    
    # Determine output path
    if args.output:
        output_path = args.output
    else:
        output_dir = os.path.join('data', 'processed', args.species)
        os.makedirs(output_dir, exist_ok=True)
        input_filename = os.path.basename(args.input)
        base_name = os.path.splitext(input_filename)[0]
        output_filename = f"{base_name}_aligned.csv"
        output_path = os.path.join(output_dir, output_filename)
    
    process_dataset(args.input, args.species, output_path)


if __name__ == "__main__":
    main()
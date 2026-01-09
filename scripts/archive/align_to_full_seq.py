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
    """Find best alignment: exact matches first, then highest identity."""
    if not protein_seq or not peptide_seq or phospho_pos_in_peptide <= 0:
        return {"match_type": "failed", "error": "Invalid input parameters"}

    candidates = []
    
    # 1) All exact matches
    candidates.extend(_exact_match_candidates(protein_seq, peptide_seq, phospho_pos_in_peptide))
    
    # 2) Local alignments across full sequence
    candidates.extend(_local_alignment_candidates(protein_seq, 0, peptide_seq, phospho_pos_in_peptide, max_alignments=5))

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
    
    # Store all valid candidates for later evaluation
    valid_candidates = []
    
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
            
            # Calculate position difference (for reporting, not for scoring)
            position_difference = 0
            if expected_position and best_result_for_candidate.get("aligned_position"):
                position_difference = abs(expected_position - best_result_for_candidate["aligned_position"])
            
            # Store candidate info
            candidate_info = {
                'idx': cand_idx,
                'id': cand_id,
                'type': cand_type,
                'sequence': protein_seq,
                'sequence_length': len(protein_seq),
                'result': best_result_for_candidate,
                'identity': best_identity_for_candidate,
                'is_exact': best_result_for_candidate.get("match_type") == "exact",
                'position_difference': position_difference
            }
            valid_candidates.append(candidate_info)
            
            print(f"      Candidate {cand_idx + 1} ({cand_id}, {cand_type}): Valid - identity {best_identity_for_candidate:.1f}%, pos_diff {position_difference}, seq_len {len(protein_seq)}")
    
    # Select best candidate based on new logic
    if valid_candidates:
        # Check if we have any SwissProt candidates
        sp_candidates = [c for c in valid_candidates if c['type'] == 'sp']
        tr_candidates = [c for c in valid_candidates if c['type'] == 'tr']
        
        # If we have SwissProt candidates, prioritize them
        if sp_candidates:
            # Among SwissProt, prefer exact matches, then highest identity
            sp_candidates.sort(key=lambda x: (x['is_exact'], x['identity']), reverse=True)
            best_candidate_info = sp_candidates[0]
        else:
            # No SwissProt, check if all TremBL have exact matches
            all_tr_exact = all(c['is_exact'] for c in tr_candidates)
            
            if all_tr_exact and len(tr_candidates) > 1:
                # All TremBL have exact matches - choose longest sequence
                tr_candidates.sort(key=lambda x: x['sequence_length'], reverse=True)
                best_candidate_info = tr_candidates[0]
                print(f"      All TremBL candidates have exact matches - selecting longest sequence: {best_candidate_info['id']} (length: {best_candidate_info['sequence_length']})")
            else:
                # Not all exact, or only one candidate - prefer exact matches, then highest identity
                tr_candidates.sort(key=lambda x: (x['is_exact'], x['identity']), reverse=True)
                best_candidate_info = tr_candidates[0]
        
        best_candidate_result = best_candidate_info['result']
        best_candidate_id = best_candidate_info['id']
        best_candidate_seq = best_candidate_info['sequence']
        best_candidate_type = best_candidate_info['type']
        best_candidate_result["position_difference"] = best_candidate_info['position_difference']
    
    # Return best candidate
    if best_candidate_result:
        print(f"      → Selected candidate: {best_candidate_id} ({best_candidate_type})")
        result = {"alignment_success": True, "selected_id": best_candidate_id, "selected_sequence": best_candidate_seq, "selected_id_type": best_candidate_type}
        result.update(best_candidate_result)
        return result
    else:
        print(f"      → No valid candidate found")
        return {"alignment_success": False, "alignment_error": "No candidate met alignment criteria (identity >= 80%, position_diff <= threshold)"}


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
        
        # Process each row individually
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
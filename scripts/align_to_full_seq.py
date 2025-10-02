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


def process_alignment_row(row: pd.Series) -> Dict[str, Any]:
    """Process a single row for alignment."""
    # Check required data - updated to work with grouped pipeline output
    if pd.isna(row.get('full_sequence')) or not row.get('full_sequence'):
        return {"alignment_success": False, "alignment_error": "Empty full sequence"}
    
    if pd.isna(row.get('cleaned_site_motif')) or not row.get('cleaned_site_motif'):
        return {"alignment_success": False, "alignment_error": "Empty cleaned site motif"}
    
    if pd.isna(row.get('motif_position')) or not row.get('motif_position'):
        return {"alignment_success": False, "alignment_error": "No motif position"}
    
    # Get data
    protein_seq = str(row['full_sequence']).strip()
    
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
    
    # Check minimum identity threshold (70%)
    if best_identity < 70.0:
        return {"alignment_success": False, "alignment_error": f"Identity too low: {best_identity:.1f}% (minimum 70%)"}
    
    # Success - return all alignment data
    result = {"alignment_success": True}
    result.update(best_result)
    
    # Calculate position difference
    if expected_position and best_result.get("aligned_position"):
        result["position_difference"] = abs(expected_position - best_result["aligned_position"])
    
    return result


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
    
    # Process each row with full sequence
    for idx, row in rows_with_sequence.iterrows():
        print(f"Processing row {idx + 1} (Protein: {row.get('Protein', 'N/A')})")
        
        alignment_result = process_alignment_row(row)
        for key, value in alignment_result.items():
            df.at[idx, key] = value
        
        if alignment_result.get('alignment_success', False):
            print(f"  ✓ Alignment successful: position {alignment_result.get('aligned_position', 'N/A')}")
        else:
            print(f"  ✗ Alignment failed: {alignment_result.get('alignment_error', 'Unknown error')}")
    
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
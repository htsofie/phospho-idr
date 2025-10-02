#!/usr/bin/env python3
"""
Remap protein IDs for groups with failed alignments and re-align all sites per group.

Workflow:
- Input: CSV produced by grouped pipeline + alignment step (has full sequences and alignment columns)
- For each Protein group with any alignment failure:
  * Try alternative identifiers available in the group's rows (UniProt → SWISS-PROT → TREMBL; no ENSEMBL)
  * For each candidate mapping: verify species, fetch full sequence, and attempt to align ALL rows in the group
  * If all sites align (≥70% identity), accept this mapping and update the group rows
  * If none work, flag group for manual review

Output: writes a new CSV with updated mapping and alignment columns
"""

import argparse
import os
import sys
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

# Add scripts directory to path to import existing modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import existing functions from other scripts
from grouped_pipeline import get_uniprot_organism, verify_species_match, get_uniprot_sequence
from align_to_full_seq import process_alignment_row


# ---------------- Remapping driver ----------------

def verify_uniprot_exists(uniprot_id: str) -> bool:
    """Check if a UniProt ID exists by trying to fetch it."""
    try:
        from Bio import ExPASy
        from Bio import SwissProt
        handle = ExPASy.get_sprot_raw(uniprot_id)
        SwissProt.read(handle)
        handle.close()
        return True
    except Exception:
        return False


def get_all_valid_uniprot_ids(id_value: str) -> List[str]:
    """Get all valid UniProt IDs from a semicolon-separated list."""
    if not id_value or pd.isna(id_value):
        return []
    ids = [x.strip() for x in str(id_value).split(';')] if ';' in str(id_value) else [str(id_value).strip()]
    valid_ids = []
    for single in ids:
        if not single:
            continue
        if verify_uniprot_exists(single):
            valid_ids.append(single)
    return valid_ids


def find_best_alignment_with_thresholds(row: pd.Series, min_identity: float, max_position_difference: Optional[int]) -> Dict[str, Any]:
    """Find the best alignment that meets both identity and position difference thresholds.
    
    If the best alignment fails position difference threshold, try alternative alignment positions.
    """
    from align_to_full_seq import map_peptide_to_protein
    
    # Check required data
    if pd.isna(row.get('full_sequence')) or not row.get('full_sequence'):
        return {"alignment_success": False, "alignment_error": "Empty full sequence"}
    
    if pd.isna(row.get('cleaned_site_motif')) or not row.get('cleaned_site_motif'):
        return {"alignment_success": False, "alignment_error": "Empty cleaned site motif"}
    
    if pd.isna(row.get('motif_position')) or not row.get('motif_position'):
        return {"alignment_success": False, "alignment_error": "No motif position"}
    
    # Get data
    protein_seq = str(row['full_sequence']).strip()
    motif_str = str(row['cleaned_site_motif']).strip()
    pos_str = str(row['motif_position']).strip()
    peptide_segments = [seg.strip() for seg in motif_str.split(';')] if ';' in motif_str else [motif_str]
    position_segments = [p.strip() for p in pos_str.split(';')] if ';' in pos_str else [pos_str]
    
    # Normalize positions list length to peptides list
    num_to_try = min(len(peptide_segments), len(position_segments))
    peptide_segments = peptide_segments[:num_to_try]
    position_segments = position_segments[:num_to_try]
    
    expected_position = int(row['position']) if not pd.isna(row.get('position')) and row.get('position') else None
    
    best_result = None
    best_identity = -1.0
    best_error = None
    
    # Try each peptide segment
    for pep, pos in zip(peptide_segments, position_segments):
        if not pep or not pos:
            continue
        try:
            phospho_pos_in_peptide = int(pos)
        except ValueError:
            continue
            
        # Get all alignment candidates for this peptide
        candidates = _get_all_alignment_candidates(protein_seq, pep, phospho_pos_in_peptide, expected_position)
        
        # Try each candidate to find one that meets thresholds
        for candidate in candidates:
            identity = candidate.get("identity", 0.0)
            
            # Check identity threshold
            if identity < min_identity:
                continue
                
            # Check position difference threshold
            if max_position_difference is not None and expected_position:
                pos_diff = abs(expected_position - candidate.get("aligned_position", 0))
                if pos_diff > max_position_difference:
                    continue
            
            # This candidate meets all thresholds
            if identity > best_identity:
                best_result = candidate
                best_identity = identity
                best_error = None
    
    if not best_result:
        return {"alignment_success": False, "alignment_error": best_error or "No alignment found meeting thresholds"}
    
    # Success - return all alignment data
    result = {"alignment_success": True}
    result.update(best_result)
    
    # Calculate position difference
    if expected_position and best_result.get("aligned_position"):
        result["position_difference"] = abs(expected_position - best_result["aligned_position"])
    
    return result


def _get_all_alignment_candidates(protein_seq: str, peptide_seq: str, phospho_pos_in_peptide: int, expected_position: Optional[int] = None) -> List[Dict[str, Any]]:
    """Get all possible alignment candidates for a peptide, sorted by quality."""
    from align_to_full_seq import _exact_match_candidates, _local_alignment_candidates
    
    candidates = []
    
    # 1) All exact matches
    candidates.extend(_exact_match_candidates(protein_seq, peptide_seq, phospho_pos_in_peptide))
    
    # 2) Local alignments across full sequence (more candidates)
    candidates.extend(_local_alignment_candidates(protein_seq, 0, peptide_seq, phospho_pos_in_peptide, max_alignments=10))
    
    # 3) Constrained windows around expected position
    if expected_position and 0 < expected_position <= len(protein_seq) * 2:
        for window in (50, 100, 200, 500):  # Added larger window
            start = max(0, expected_position - window)
            end = min(len(protein_seq), expected_position + window)
            if end - start >= max(10, len(peptide_seq)):
                region = protein_seq[start:end]
                candidates.extend(_local_alignment_candidates(region, start, peptide_seq, phospho_pos_in_peptide, max_alignments=10))
    
    # 4) Additional windows across the protein for better coverage
    if len(protein_seq) > 1000:  # Only for longer proteins
        step_size = max(200, len(peptide_seq) * 3)
        for start in range(0, len(protein_seq) - len(peptide_seq), step_size):
            end = min(len(protein_seq), start + len(peptide_seq) * 4)
            region = protein_seq[start:end]
            candidates.extend(_local_alignment_candidates(region, start, peptide_seq, phospho_pos_in_peptide, max_alignments=5))
    
    # Filter valid candidates
    candidates = [c for c in candidates if c.get("aligned_position")]
    
    # Sort by priority: exact > identity > proximity to expected position
    def sort_key(c: Dict[str, Any]) -> Tuple[int, float, float]:
        is_exact = 1 if c["match_type"] == "exact" else 0
        identity = c.get("identity", 0.0)
        proximity = -abs(expected_position - c["aligned_position"]) if expected_position and c.get("aligned_position") else 0.0
        return (is_exact, identity, proximity)
    
    candidates.sort(key=sort_key, reverse=True)
    return candidates


def perform_alignment_with_threshold(row: pd.Series, min_identity: float, max_position_difference: Optional[int]) -> Dict[str, Any]:
    """Wrapper for backward compatibility - now uses the enhanced alignment function."""
    return find_best_alignment_with_thresholds(row, min_identity, max_position_difference)


def try_group_mapping_and_alignment(group_df: pd.DataFrame, species: str, min_identity: float, max_position_difference: Optional[int]) -> Optional[Dict[str, Any]]:
    """Try alternate IDs for a group. Return dict with keys if success, else None.
    
    Logic:
    1. For each column type (uniprot → swissprot → trembl), try ALL semicolon-separated IDs from ALL rows
    2. Only move to next column type if ALL IDs from current column fail
    3. For each ID, verify species and try to align ALL rows in the group
    4. Accept first ID where ALL rows align successfully
    """
    # Define column priority. Try UniProt first, then SWISS-PROT, then TREMBL for all species
    column_priority = [('uniprot', 'UniProt'), ('swissprot_id', 'SWISS-PROT'), ('trembl_id', 'TREMBL')]
    
    # Try each column type in priority order
    for col, label in column_priority:
        if col not in group_df.columns:
            continue
            
        print(f"  Trying {label} IDs...")
        
        # Build ordered candidate lists per row (preserve within-row order)
        row_to_ids: Dict[int, List[str]] = {}
        max_ids_len = 0
        for idx, row in group_df.iterrows():
            cell_value = row.get(col)
            ids_for_row = get_all_valid_uniprot_ids(str(cell_value)) if not pd.isna(cell_value) and cell_value else []
            row_to_ids[idx] = ids_for_row
            if len(ids_for_row) > max_ids_len:
                max_ids_len = len(ids_for_row)

        if max_ids_len == 0:
            print(f"    No valid {label} IDs found in group")
            continue

        # Iterate position-wise: first values across rows, then second, etc.
        from collections import OrderedDict
        for position_index in range(max_ids_len):
            # Gather the position-th ID from each row (deduplicated, preserve encounter order)
            ordered_candidates: "OrderedDict[str, None]" = OrderedDict()
            for idx in group_df.index:
                ids_for_row = row_to_ids.get(idx, [])
                if position_index < len(ids_for_row):
                    ordered_candidates.setdefault(ids_for_row[position_index], None)

            if not ordered_candidates:
                continue

            print(f"    Trying {label} candidate set at position {position_index + 1}: {list(ordered_candidates.keys())}")

            # Try each ID from this position across the group
            for uniprot_id in ordered_candidates.keys():
                print(f"    Trying {label} ID: {uniprot_id}")
                
                # Verify species match
                if not verify_species_match(uniprot_id, species):
                    print(f"      ✗ Species mismatch")
                    continue
                
                # Fetch sequence
                seq = get_uniprot_sequence(uniprot_id)
                if not seq:
                    print(f"      ✗ Failed to fetch sequence")
                    continue
                
                full_sequence, seq_type = seq
                print(f"      ✓ Got sequence (length: {len(full_sequence)}, type: {seq_type})")
                
                # Try to align ALL rows in the group
                per_row_results: Dict[int, Dict[str, Any]] = {}
                all_ok = True
                
                for idx, row in group_df.iterrows():
                    temp_row = row.copy()
                    temp_row['full_sequence'] = full_sequence
                    aln = perform_alignment_with_threshold(temp_row, min_identity, max_position_difference)
                    per_row_results[idx] = aln
                    
                    if not aln.get('alignment_success', False):
                        print(f"      ✗ Alignment failed for row {idx}: {aln.get('alignment_error', 'Unknown error')}")
                        all_ok = False
                        break
                    else:
                        identity = aln.get('identity', 0.0)
                        print(f"      ✓ Row {idx} aligned (identity: {identity:.1f}%)")
                
                if all_ok:
                    print(f"      ✓ SUCCESS: All rows aligned with {label} ID {uniprot_id}")
                    # Clear alignment errors for successful alignments
                    for row_idx, aln in per_row_results.items():
                        if aln.get('alignment_success', False):
                            aln['alignment_error'] = ''
                    return {
                        'uniprot_mapped_id': uniprot_id,
                        'mapping_source': label,
                        'full_sequence': full_sequence,
                        'sequence_length': len(full_sequence),
                        'sequence_type': seq_type,
                        'per_row_results': per_row_results,
                    }
                else:
                    print(f"      ✗ Failed: Not all rows aligned")
    
    print(f"  ✗ All ID types failed for group")
    return None


def process_dataset(input_path: str, species: str, output_path: str, min_identity: float = 70.0, max_position_difference: Optional[int] = 150) -> None:
    print(f"Loading dataset: {input_path}")
    df = pd.read_csv(input_path)
    print(f"Loaded {len(df)} rows")

    # Quick checks
    required_cols = ['Protein', 'cleaned_site_motif', 'motif_position']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"Missing required columns: {missing}")
        return

    # Identify groups that need remapping: any alignment_success == False (or NaN) in group
    needs_remap_groups: List[str] = []
    if 'alignment_success' in df.columns:
        group_bad = df.groupby('Protein')['alignment_success'].apply(lambda s: (~s.fillna(False)).any())
        needs_remap_groups = group_bad[group_bad].index.tolist()
    else:
        # If no alignment columns, remap nothing
        print("No alignment_success column found; nothing to remap.")
        needs_remap_groups = []

    print(f"Groups needing remap: {len(needs_remap_groups)}")

    remapped_groups = 0
    failed_groups = 0

    for protein_id in needs_remap_groups:
        group_idx = df[df['Protein'] == protein_id].index
        group_df = df.loc[group_idx]
        print(f"Remapping group Protein={protein_id} (rows={len(group_df)})")
        result = try_group_mapping_and_alignment(group_df, species, min_identity, max_position_difference)
        if result:
            remapped_groups += 1
            # Apply mapping and per-row alignment results
            df.loc[group_idx, 'uniprot_mapped_id'] = result['uniprot_mapped_id']
            df.loc[group_idx, 'mapping_source'] = result['mapping_source']
            df.loc[group_idx, 'full_sequence'] = result['full_sequence']
            df.loc[group_idx, 'sequence_length'] = result['sequence_length']
            df.loc[group_idx, 'sequence_type'] = result['sequence_type']
            for row_idx, aln in result['per_row_results'].items():
                for key, val in aln.items():
                    df.at[row_idx, key] = val
            # Clear manual review flag if present
            if 'manual_review_flag' in df.columns:
                df.loc[group_idx, 'manual_review_flag'] = False
            if 'processing_notes' in df.columns:
                df.loc[group_idx, 'processing_notes'] = df.loc[group_idx, 'processing_notes'].fillna('').apply(
                    lambda x: (str(x) + '; remapped_ok') if x else 'remapped_ok')
            print("  ✓ Remap successful for group")
        else:
            failed_groups += 1
            if 'manual_review_flag' in df.columns:
                df.loc[group_idx, 'manual_review_flag'] = True
            if 'processing_notes' in df.columns:
                df.loc[group_idx, 'processing_notes'] = df.loc[group_idx, 'processing_notes'].fillna('').apply(
                    lambda x: (str(x) + '; remap_failed') if x else 'remap_failed')
            print("  ✗ Remap failed; flagged for manual review")

    # Save
    print(f"\nSaving results to: {output_path}")
    df.to_csv(output_path, index=False)

    # Summary
    print("\nRemap summary:")
    print(f"  Groups needing remap: {len(needs_remap_groups)}")
    print(f"  Remapped successfully: {remapped_groups}")
    print(f"  Failed to remap: {failed_groups}")


def main() -> None:
    parser = argparse.ArgumentParser(description='Remap ID and re-align groups that failed alignment')
    parser.add_argument('--input', '-i', required=True, help='Path to input CSV (aligned output)')
    parser.add_argument('--species', '-s', required=True, choices=['mouse', 'rat'], help='Species name')
    parser.add_argument('--min-identity', type=float, default=70.0, help='Minimum alignment identity threshold')
    parser.add_argument('--max-pos-diff', type=int, default=150, help='Maximum allowed absolute position difference before rejecting a candidate')
    parser.add_argument('--output', '-o', help='Path to output CSV (default: *_remapped.csv)')
    args = parser.parse_args()

    if args.output:
        output_path = args.output
    else:
        base = os.path.splitext(os.path.basename(args.input))[0]
        output_dir = os.path.dirname(args.input) or '.'
        output_path = os.path.join(output_dir, f"{base}_remapped.csv")

    process_dataset(args.input, args.species, output_path, min_identity=args.min_identity, max_position_difference=args.max_pos_diff)


if __name__ == '__main__':
    main()



#!/usr/bin/env python3
"""
Script to calculate identity at specific positions between two aligned sequences.
Parses CLUSTAL format alignment files and position lists.

Usage:
    python scripts/calculate_alignment_identity.py \
        --alignment data/processed/mouse/aligned.txt \
        --positions "[33-35,60,62-63,94-96,144-146,149-151,171,173-179,186,231-232,245-251,258-259]" \
        --reference P42209 \
        --output results.csv
"""

import re
import argparse
import csv
from typing import List, Tuple, Dict, Any
from pathlib import Path


def parse_position_list(position_str: str) -> List[int]:
    """
    Parse position list in format: [33-35,60,62-63,94-96,...]
    
    Args:
        position_str: String like "[33-35,60,62-63,94-96]"
        
    Returns:
        List of all positions (e.g., [33, 34, 35, 60, 62, 63, 94, 95, 96])
    """
    # Remove brackets and whitespace
    position_str = position_str.strip().strip('[]')
    
    if not position_str:
        return []
    
    positions = []
    # Split by comma
    for part in position_str.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            # Range like "33-35"
            try:
                start, end = map(int, part.split('-'))
                positions.extend(range(start, end + 1))
            except ValueError:
                print(f"Warning: Could not parse range '{part}', skipping")
        else:
            # Single position
            try:
                positions.append(int(part))
            except ValueError:
                print(f"Warning: Could not parse position '{part}', skipping")
    
    return sorted(positions)


def parse_clustal_alignment(file_path: str) -> Dict[str, str]:
    """
    Parse CLUSTAL format alignment file.
    
    Args:
        file_path: Path to CLUSTAL alignment file
        
    Returns:
        Dictionary mapping sequence IDs to aligned sequences (with gaps)
    """
    sequences = {}
    current_id = None
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.rstrip()  # Keep right whitespace for alignment
            
            # Skip header and empty lines
            if not line or line.startswith('CLUSTAL'):
                continue
            
            # Skip alignment quality lines (asterisks, colons, dots)
            if line.strip().startswith(('*', ':', '.')) or not line.strip():
                continue
            
            # Check if this is a sequence line
            # Format: ID      SEQUENCE     POSITION
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                seq_part = parts[1]
                
                # Check if second part looks like a sequence (contains amino acids or gaps)
                if any(c in seq_part for c in 'ACDEFGHIKLMNPQRSTVWY-'):
                    if seq_id not in sequences:
                        sequences[seq_id] = []
                    sequences[seq_id].append(seq_part)
    
    # Concatenate all sequence parts for each ID
    for seq_id in sequences:
        sequences[seq_id] = ''.join(sequences[seq_id])
    
    return sequences


def map_position_to_alignment(seq_with_gaps: str, position: int) -> int:
    """
    Map a position in the original sequence to its position in the aligned sequence.
    
    Args:
        seq_with_gaps: Aligned sequence with gaps ('-')
        position: Position in original sequence (1-indexed)
        
    Returns:
        Position in aligned sequence (0-indexed), or -1 if not found
    """
    # Count non-gap characters
    non_gap_count = 0
    for i, char in enumerate(seq_with_gaps):
        if char != '-':
            non_gap_count += 1
            if non_gap_count == position:
                return i
    
    return -1


def map_alignment_to_position(seq_with_gaps: str, aligned_pos: int) -> int:
    """
    Map an aligned position (in alignment string) to position in original sequence (non-gap counting).
    
    Args:
        seq_with_gaps: Aligned sequence with gaps ('-')
        aligned_pos: Position in aligned sequence (0-indexed)
        
    Returns:
        Position in original sequence (1-indexed, non-gap counting), or -1 if gap or invalid
    """
    if aligned_pos < 0 or aligned_pos >= len(seq_with_gaps):
        return -1
    
    if seq_with_gaps[aligned_pos] == '-':
        return -1
    
    # Count non-gap characters up to and including this position
    non_gap_count = 0
    for i in range(aligned_pos + 1):
        if seq_with_gaps[i] != '-':
            non_gap_count += 1
    
    return non_gap_count


def calculate_identity_at_positions(
    seq1_aligned: str,
    seq2_aligned: str,
    positions: List[int],
    reference_seq_id: str,
    aligned_seq_id: str,
    sequences: Dict[str, str]
) -> Dict[str, Any]:
    """
    Calculate identity at specific positions between two aligned sequences.
    
    Args:
        seq1_aligned: First aligned sequence
        seq2_aligned: Second aligned sequence
        positions: List of positions in reference sequence (1-indexed)
        reference_seq_id: ID of the reference sequence (positions refer to this)
        aligned_seq_id: ID of the aligned sequence (non-reference, to show position)
        sequences: Dictionary of all sequences (to get reference and aligned)
        
    Returns:
        Dictionary with match results and identity percentage
    """
    if reference_seq_id not in sequences:
        raise ValueError(f"Reference sequence ID '{reference_seq_id}' not found in sequences")
    if aligned_seq_id not in sequences:
        raise ValueError(f"Aligned sequence ID '{aligned_seq_id}' not found in sequences")
    
    reference_seq = sequences[reference_seq_id]
    aligned_seq = sequences[aligned_seq_id]
    
    matches = 0
    mismatches = 0
    gaps = 0
    not_found = 0
    position_results = []
    
    for pos in positions:
        # Map position to alignment coordinates
        aligned_pos = map_position_to_alignment(reference_seq, pos)
        
        if aligned_pos == -1:
            not_found += 1
            position_results.append({
                'position': pos,
                'aligned_position': None,
                'status': 'not_found',
                'ref_aa': None,
                'aligned_aa': None
            })
            continue
        
        # Get amino acids at this aligned position
        ref_aa = reference_seq[aligned_pos] if aligned_pos < len(reference_seq) else None
        aligned_aa = aligned_seq[aligned_pos] if aligned_pos < len(aligned_seq) else None
        
        # Map aligned position back to position in the aligned sequence (1-indexed, non-gap counting)
        aligned_pos_in_seq = map_alignment_to_position(aligned_seq, aligned_pos)
        
        # Check for gaps
        if ref_aa == '-' or aligned_aa == '-':
            gaps += 1
            status = 'gap'
        elif ref_aa == aligned_aa:
            matches += 1
            status = 'match'
        else:
            mismatches += 1
            status = 'mismatch'
        
        position_results.append({
            'position': pos,  # Position in reference sequence (1-indexed)
            'aligned_position': aligned_pos_in_seq,  # Position in aligned sequence (1-indexed, non-gap)
            'status': status,
            'ref_aa': ref_aa,
            'aligned_aa': aligned_aa
        })
    
    # Calculate identity (excluding gaps and not found)
    total_compared = matches + mismatches
    identity = (matches / total_compared * 100.0) if total_compared > 0 else 0.0
    
    return {
        'matches': matches,
        'mismatches': mismatches,
        'gaps': gaps,
        'not_found': not_found,
        'total_positions': len(positions),
        'identity_percentage': identity,
        'position_results': position_results
    }


def print_results(results: Dict[str, Any], reference_id: str, aligned_id: str):
    """Print results to console."""
    print(f"\n{'='*70}")
    print("IDENTITY CALCULATION RESULTS")
    print(f"{'='*70}")
    print(f"Reference sequence: {reference_id}")
    print(f"Aligned sequence: {aligned_id}")
    print(f"\nSummary:")
    print(f"  Matches: {results['matches']}")
    print(f"  Mismatches: {results['mismatches']}")
    print(f"  Gaps: {results['gaps']}")
    print(f"  Not found: {results['not_found']}")
    print(f"  Total positions checked: {results['total_positions']}")
    print(f"  Identity: {results['identity_percentage']:.2f}%")
    
    print(f"\nDetailed position results:")
    print(f"{'RefPos':<8} {'AlignedPos':<12} {'Status':<12} {'Ref_AA':<8} {'Aligned_AA':<10}")
    print("-" * 60)
    for res in results['position_results']:
        aligned_pos = res.get('aligned_position', 'N/A')
        aligned_pos_str = str(aligned_pos) if aligned_pos != 'N/A' and aligned_pos is not None else 'N/A'
        print(f"{res['position']:<8} {aligned_pos_str:<12} {res['status']:<12} "
              f"{res.get('ref_aa', 'N/A'):<8} {res.get('aligned_aa', 'N/A'):<10}")


def save_results_to_csv(results: Dict[str, Any], reference_id: str, aligned_id: str, 
                        output_path: str):
    """Save results to CSV file."""
    with open(output_path, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Write header
        writer.writerow(['Identity Analysis Results'])
        writer.writerow(['Reference sequence', reference_id])
        writer.writerow(['Aligned sequence', aligned_id])
        writer.writerow([])
        writer.writerow(['Summary'])
        writer.writerow(['Matches', results['matches']])
        writer.writerow(['Mismatches', results['mismatches']])
        writer.writerow(['Gaps', results['gaps']])
        writer.writerow(['Not found', results['not_found']])
        writer.writerow(['Total positions', results['total_positions']])
        writer.writerow(['Identity percentage', f"{results['identity_percentage']:.2f}%"])
        writer.writerow([])
        
        # Write detailed results
        writer.writerow(['Reference_Position', 'Aligned_Position', 'Status', 'Ref_AA', 'Aligned_AA'])
        for res in results['position_results']:
            writer.writerow([
                res['position'],
                res.get('aligned_position', 'N/A'),
                res['status'],
                res.get('ref_aa', 'N/A'),
                res.get('aligned_aa', 'N/A')
            ])
    
    print(f"\nResults saved to: {output_path}")


def main():
    """Main function to run the analysis."""
    parser = argparse.ArgumentParser(
        description='Calculate identity at specific positions between aligned sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python scripts/calculate_alignment_identity.py \\
      --alignment data/processed/mouse/aligned.txt \\
      --positions "[33-35,60,62-63]" \\
      --reference P42209
  
  # With output file
  python scripts/calculate_alignment_identity.py \\
      --alignment data/processed/mouse/aligned.txt \\
      --positions "[33-35,60,62-63,94-96]" \\
      --reference P42209 \\
      --output results.csv
        """
    )
    
    parser.add_argument(
        '--alignment',
        type=str,
        required=True,
        help='Path to CLUSTAL format alignment file'
    )
    
    parser.add_argument(
        '--positions',
        type=str,
        required=True,
        help='Position list in format: "[33-35,60,62-63,94-96]"'
    )
    
    parser.add_argument(
        '--reference',
        type=str,
        required=True,
        help='ID of the reference sequence (positions refer to this sequence)'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Optional: Path to output CSV file'
    )
    
    args = parser.parse_args()
    
    # Check if alignment file exists
    if not Path(args.alignment).exists():
        print(f"Error: Alignment file not found: {args.alignment}")
        return 1
    
    # Parse alignment
    print(f"Parsing alignment from {args.alignment}...")
    try:
        sequences = parse_clustal_alignment(args.alignment)
    except Exception as e:
        print(f"Error parsing alignment file: {e}")
        return 1
    
    if len(sequences) < 2:
        print(f"Error: Need at least 2 sequences, found {len(sequences)}")
        print(f"Found sequences: {list(sequences.keys())}")
        return 1
    
    # Get sequence IDs
    seq_ids = list(sequences.keys())
    seq1_id = seq_ids[0]
    seq2_id = seq_ids[1]
    
    # Determine which sequence is the non-reference (aligned) sequence
    if args.reference == seq1_id:
        aligned_seq_id = seq2_id
    elif args.reference == seq2_id:
        aligned_seq_id = seq1_id
    else:
        print(f"Error: Reference sequence '{args.reference}' not found in alignment")
        print(f"Available sequences: {seq_ids}")
        return 1
    
    print(f"Found sequences: {seq_ids}")
    print(f"Reference sequence ({args.reference}): {len(sequences[args.reference])} characters")
    print(f"Aligned sequence ({aligned_seq_id}): {len(sequences[aligned_seq_id])} characters")
    
    # Parse positions
    try:
        positions = parse_position_list(args.positions)
    except Exception as e:
        print(f"Error parsing position list: {e}")
        return 1
    
    if not positions:
        print("Error: No valid positions found in position list")
        return 1
    
    print(f"\nParsed positions: {positions}")
    print(f"Total positions to check: {len(positions)}")
    
    # Calculate identity
    print(f"\nUsing {args.reference} as reference sequence...")
    try:
        results = calculate_identity_at_positions(
            sequences[seq1_id],
            sequences[seq2_id],
            positions,
            args.reference,
            aligned_seq_id,
            sequences
        )
    except Exception as e:
        print(f"Error calculating identity: {e}")
        return 1
    
    # Print results
    print_results(results, args.reference, aligned_seq_id)
    
    # Save to CSV if requested
    if args.output:
        try:
            save_results_to_csv(results, args.reference, aligned_seq_id, args.output)
        except Exception as e:
            print(f"Error saving results to CSV: {e}")
            return 1
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())


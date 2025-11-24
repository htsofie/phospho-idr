#!/usr/bin/env python3
"""
Script to align mouse/rat proteins with human orthologs using EBI EMBOSS needle API
and map interface residues (IRES) from human positions to mouse/rat positions.

Usage:
    python scripts/align_interfaces.py --input data/processed/mouse/int_interfaces.csv --species mouse
    python scripts/align_interfaces.py --input data/processed/rat/int_interfaces.csv --species rat
"""

import re
import argparse
import logging
import time
import requests
import pandas as pd
from typing import List, Dict, Any, Optional
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# EBI REST API base URL
EBI_BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/emboss_needle"


def format_fasta_sequence(seq_id: str, sequence: str) -> str:
    """
    Format sequence in FASTA format.
    
    Args:
        seq_id: Sequence identifier
        sequence: Sequence data
        
    Returns:
        FASTA-formatted string: ">seq_id\nsequence"
    """
    return f">{seq_id}\n{sequence}"


def submit_emboss_needle_job(asequence: str, bsequence: str, email: str = "user@example.com") -> str:
    """
    Submit EMBOSS needle alignment job to EBI REST API.
    
    Args:
        asequence: FASTA-formatted sequence (protein_id + full_sequence)
        bsequence: FASTA-formatted sequence (human_ortholog_id + human_ortholog_sequence)
        email: Email address (required by API)
        
    Returns:
        Job ID string
        
    Raises:
        Exception: If job submission fails
    """
    url = f"{EBI_BASE_URL}/run"
    data = {
        'asequence': asequence,
        'bsequence': bsequence,
        'email': email
    }
    
    try:
        response = requests.post(url, data=data, timeout=30)
        response.raise_for_status()
        job_id = response.text.strip()
        logger.debug(f"Submitted job, got ID: {job_id}")
        return job_id
    except requests.HTTPError as e:
        # Get the actual error message from the API
        error_detail = ""
        if hasattr(e, 'response') and e.response is not None:
            error_detail = e.response.text
            logger.error(f"API error response: {error_detail}")
        raise Exception(f"Failed to submit alignment job: {e} - API response: {error_detail}")
    except requests.RequestException as e:
        raise Exception(f"Failed to submit alignment job: {e}")


def check_job_status(job_id: str) -> str:
    """
    Check status of EMBOSS needle alignment job.
    
    Args:
        job_id: Job ID from submission
        
    Returns:
        Status string ("RUNNING", "FINISHED", "ERROR", etc.)
        
    Raises:
        Exception: If status check fails
    """
    url = f"{EBI_BASE_URL}/status/{job_id}"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        status = response.text.strip()
        return status
    except requests.RequestException as e:
        raise Exception(f"Failed to check job status: {e}")


def get_alignment_result(job_id: str) -> str:
    """
    Retrieve alignment result from EBI REST API.
    
    Args:
        job_id: Job ID from submission
        
    Returns:
        Alignment text in EMBOSS format
        
    Raises:
        Exception: If result retrieval fails
    """
    url = f"{EBI_BASE_URL}/result/{job_id}/aln"
    
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return response.text
    except requests.RequestException as e:
        raise Exception(f"Failed to retrieve alignment result: {e}")


def parse_emboss_alignment(alignment_text: str) -> Dict[str, str]:
    """
    Parse EMBOSS alignment output format.
    
    Args:
        alignment_text: Alignment text from EMBOSS API
        
    Returns:
        Dictionary mapping sequence IDs to aligned sequences (with gaps)
        Format: {seq1_id: aligned_seq1, seq2_id: aligned_seq2}
        
    Raises:
        Exception: If parsing fails
    """
    sequences = {}
    lines = alignment_text.split('\n')
    
    # First, extract sequence IDs from headers like "# 1: P42209" and "# 2: Q8WYJ6"
    seq_ids = []
    for line in lines:
        # Look for pattern: "# 1: P42209" or "# 2: Q8WYJ6"
        match = re.search(r'#\s*(\d+):\s*([A-Z0-9_]+)', line)
        if match:
            seq_num = int(match.group(1))
            seq_id = match.group(2)
            # Store in order (seq_num is 1-indexed, so subtract 1 for list index)
            while len(seq_ids) < seq_num:
                seq_ids.append(None)
            seq_ids[seq_num - 1] = seq_id
    
    if len(seq_ids) != 2:
        raise Exception(f"Failed to find 2 sequence IDs in headers. Found: {seq_ids}")
    
    # Initialize sequences
    sequences = {seq_ids[0]: [], seq_ids[1]: []}
    
    # Now extract sequences from alternating blocks
    # Format: "P42209             1 -----MDKEYVGFAALPNQLHRKSVKKGFDFTLMVAGESGLGKSTLINSL     45"
    #         "Q8WYJ6             1 MAGGVMDKEYVGFAALPNQLHRKSVKKGFDFTLMVAGESGLGKSTLINSL     50"
    # Sequences alternate, and each line has: ID, start_pos, sequence, end_pos
    
    for line in lines:
        line = line.rstrip()
        
        # Skip comment lines and empty lines
        if line.startswith('#') or not line.strip():
            continue
        
        # Check if this line contains sequence data
        if not any(c in line for c in 'ACDEFGHIKLMNPQRSTVWY-'):
            continue
        
        # Try to match sequence line format: "P42209             1 -----MDKEY...     45"
        # Pattern: sequence_id, whitespace, start_number, sequence, end_number
        match = re.match(r'^([A-Z0-9_]+)\s+\d+\s+([ACDEFGHIKLMNPQRSTVWY\-]+)\s+\d+', line)
        if match:
            seq_id = match.group(1)
            seq_part = match.group(2)
            
            # Add to the appropriate sequence
            if seq_id in sequences:
                sequences[seq_id].append(seq_part)
        else:
            # Try alternative: sequence might be in middle with numbers on both sides
            # Extract just the sequence part (amino acids and gaps)
            parts = line.split()
            if len(parts) >= 3:
                # First part should be ID, last part should be end position
                # Middle parts might contain the sequence
                potential_seq = ''.join(parts[1:-1])  # Everything between ID and end position
                # Check if it looks like a sequence
                if all(c in 'ACDEFGHIKLMNPQRSTVWY-' for c in potential_seq):
                    seq_id = parts[0]
                    if seq_id in sequences:
                        sequences[seq_id].append(potential_seq)
    
    # Join all parts for each sequence
    for seq_id in sequences:
        sequences[seq_id] = ''.join(sequences[seq_id])
    
    # Validate we have both sequences
    if len(sequences) != 2 or any(not seq for seq in sequences.values()):
        raise Exception(f"Failed to parse alignment: expected 2 sequences, found {len(sequences)}. "
                       f"Sequences found: {list(sequences.keys())}. "
                       f"Lengths: {[len(s) for s in sequences.values()]}")
    
    return sequences


def call_ebi_emboss_needle(seq1_id: str, seq1: str, seq2_id: str, seq2: str, 
                           email: str = "user@example.com", max_poll_attempts: int = 60) -> Dict[str, str]:
    """
    Main function to orchestrate EBI API calls for EMBOSS needle alignment.
    
    Args:
        seq1_id: First sequence ID (protein_id)
        seq1: First sequence (full_sequence)
        seq2_id: Second sequence ID (human_ortholog_id)
        seq2: Second sequence (human_ortholog_sequence)
        email: Email address for API
        max_poll_attempts: Maximum number of polling attempts
        
    Returns:
        Dictionary: {seq1_id: aligned_seq1, seq2_id: aligned_seq2}
        
    Raises:
        Exception: If alignment fails
    """
    # Format sequences in FASTA
    asequence = format_fasta_sequence(seq1_id, seq1)
    bsequence = format_fasta_sequence(seq2_id, seq2)
    
    # Submit job
    logger.info(f"Submitting alignment job for {seq1_id} vs {seq2_id}...")
    time.sleep(0.5)  # Rate limiting
    job_id = submit_emboss_needle_job(asequence, bsequence, email)
    logger.debug(f"Job ID: {job_id}")
    
    # Poll for completion
    logger.info(f"Waiting for alignment to complete...")
    for attempt in range(max_poll_attempts):
        time.sleep(2)  # Wait 2 seconds between checks
        status = check_job_status(job_id)
        logger.debug(f"Job status (attempt {attempt + 1}): {status}")
        
        if status == "FINISHED":
            break
        elif status == "ERROR":
            raise Exception(f"Alignment job failed with ERROR status")
        elif status not in ["RUNNING", "PENDING", "QUEUED"]:
            raise Exception(f"Unexpected job status: {status}")
    
    if status != "FINISHED":
        raise Exception(f"Alignment job timed out after {max_poll_attempts} attempts")
    
    # Retrieve result
    logger.info(f"Retrieving alignment result...")
    time.sleep(0.5)  # Rate limiting
    alignment_text = get_alignment_result(job_id)
    
    # Parse alignment
    sequences = parse_emboss_alignment(alignment_text)
    
    # Ensure we have the right sequence IDs
    # EMBOSS may return sequences in different order or with modified IDs
    # Map back to our original IDs
    aligned_seqs = {}
    seq_ids_in_result = list(sequences.keys())
    
    # Try to match by ID
    if seq1_id in sequences:
        aligned_seqs[seq1_id] = sequences[seq1_id]
    elif seq_ids_in_result:
        aligned_seqs[seq1_id] = sequences[seq_ids_in_result[0]]
    
    if seq2_id in sequences:
        aligned_seqs[seq2_id] = sequences[seq2_id]
    elif len(seq_ids_in_result) > 1:
        aligned_seqs[seq2_id] = sequences[seq_ids_in_result[1]]
    elif len(seq_ids_in_result) == 1:
        # Only one sequence found, use it for seq2
        aligned_seqs[seq2_id] = sequences[seq_ids_in_result[0]]
    
    if len(aligned_seqs) != 2:
        # Fallback: use sequences in order
        aligned_seqs = {
            seq1_id: sequences[seq_ids_in_result[0]] if seq_ids_in_result else '',
            seq2_id: sequences[seq_ids_in_result[1]] if len(seq_ids_in_result) > 1 else sequences[seq_ids_in_result[0]] if seq_ids_in_result else ''
        }
    
    logger.info(f"Alignment complete: {len(aligned_seqs[seq1_id])} characters")
    return aligned_seqs


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
                logger.warning(f"Could not parse range '{part}', skipping")
        else:
            # Single position
            try:
                positions.append(int(part))
            except ValueError:
                logger.warning(f"Could not parse position '{part}', skipping")
    
    return sorted(positions)


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


def map_positions_via_alignment(human_positions: List[int], human_aligned: str, 
                                mouse_aligned: str, human_id: str, mouse_id: str) -> Dict[str, Any]:
    """
    Map human IRES positions to mouse/rat positions using alignment.
    
    Args:
        human_positions: List of positions in human sequence (1-indexed)
        human_aligned: Human aligned sequence (with gaps)
        mouse_aligned: Mouse/rat aligned sequence (with gaps)
        human_id: Human sequence ID (for logging)
        mouse_id: Mouse/rat sequence ID (for logging)
        
    Returns:
        Dictionary with:
        - 'mapped_positions': List of mapped positions (None for gaps/unmappable)
        - 'matched_positions': List of positions where amino acids match
        - 'mismatched_positions': List of positions where amino acids don't match
        - 'gap_positions': List of positions that map to gaps
    """
    mapped_positions = []
    matched_positions = []
    mismatched_positions = []
    gap_positions = []
    
    for human_pos in human_positions:
        # Map human position to alignment coordinate
        aligned_pos = map_position_to_alignment(human_aligned, human_pos)
        
        if aligned_pos == -1:
            # Position not found in human sequence
            mapped_positions.append(None)
            gap_positions.append(human_pos)
            continue
        
        # Get amino acids at this aligned position
        if aligned_pos >= len(human_aligned) or aligned_pos >= len(mouse_aligned):
            mapped_positions.append(None)
            gap_positions.append(human_pos)
            continue
        
        human_aa = human_aligned[aligned_pos]
        mouse_aa = mouse_aligned[aligned_pos]
        
        # Check for gaps
        if human_aa == '-' or mouse_aa == '-':
            mapped_positions.append(None)
            gap_positions.append(human_pos)
            continue
        
        # Map to mouse sequence position (1-indexed, non-gap)
        mouse_pos = map_alignment_to_position(mouse_aligned, aligned_pos)
        
        if mouse_pos == -1:
            mapped_positions.append(None)
            gap_positions.append(human_pos)
            continue
        
        mapped_positions.append(mouse_pos)
        
        # Check if amino acids match
        if human_aa == mouse_aa:
            matched_positions.append(mouse_pos)
        else:
            mismatched_positions.append(mouse_pos)
    
    return {
        'mapped_positions': mapped_positions,
        'matched_positions': sorted(matched_positions),
        'mismatched_positions': sorted(mismatched_positions),
        'gap_positions': gap_positions
    }


def calculate_alignment_identity(human_aligned: str, mouse_aligned: str) -> float:
    """
    Calculate percent identity between aligned sequences.
    
    Args:
        human_aligned: Human aligned sequence (with gaps)
        mouse_aligned: Mouse/rat aligned sequence (with gaps)
        
    Returns:
        Percent identity (0-100)
    """
    if len(human_aligned) != len(mouse_aligned):
        logger.warning(f"Alignment length mismatch: {len(human_aligned)} vs {len(mouse_aligned)}")
        return 0.0
    
    matches = 0
    total = 0
    
    for i in range(len(human_aligned)):
        human_aa = human_aligned[i]
        mouse_aa = mouse_aligned[i]
        
        # Skip positions where either sequence has a gap
        if human_aa == '-' or mouse_aa == '-':
            continue
        
        total += 1
        if human_aa == mouse_aa:
            matches += 1
    
    if total == 0:
        return 0.0
    
    return (matches / total) * 100.0


def format_ires_output(positions: List[int]) -> str:
    """
    Convert list of positions to IRES format: [33-35,60,62-63,...]
    
    Args:
        positions: List of positions (may include None for gaps)
        
    Returns:
        Formatted string with ranges for consecutive positions
    """
    # Filter out None values
    valid_positions = [p for p in positions if p is not None]
    
    if not valid_positions:
        return "[]"
    
    # Sort and group consecutive positions
    valid_positions = sorted(set(valid_positions))
    ranges = []
    start = valid_positions[0]
    end = valid_positions[0]
    
    for i in range(1, len(valid_positions)):
        if valid_positions[i] == end + 1:
            # Consecutive, extend range
            end = valid_positions[i]
        else:
            # Gap found, save current range
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{end}")
            start = valid_positions[i]
            end = valid_positions[i]
    
    # Add last range
    if start == end:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}-{end}")
    
    return "[" + ",".join(ranges) + "]"


def process_row_alignment(row: pd.Series, species: str, email: str = "user@example.com") -> Dict[str, Any]:
    """
    Process single row from int_interfaces.csv: align sequences and map IRES positions.
    
    Args:
        row: DataFrame row with protein data
        species: Species name (mouse or rat)
        email: Email for EBI API
        
    Returns:
        Dictionary with new column values:
        - 'aligned_IRES': Mapped positions
        - 'alignment_percent': Percent identity
        - 'matched_IRES': Matching positions
        - 'mismatched_IRES': Mismatching positions
        - 'alignment_error': Error message if any
    """
    result = {
        'aligned_IRES': None,
        'alignment_percent': None,
        'matched_IRES': None,
        'mismatched_IRES': None,
        'alignment_error': None
    }
    
    try:
        # Extract required data
        protein_id = row.get('protein_id')
        full_sequence = row.get('full_sequence')
        human_ortholog_id = row.get('human_ortholog_id')
        human_ortholog_sequence = row.get('human_ortholog_sequence')
        human_phosphoprotein_IRES = row.get('human_phosphoprotein_IRES')
        
        # Validate required data
        if pd.isna(protein_id) or not protein_id:
            result['alignment_error'] = "Missing protein_id"
            return result
        
        if pd.isna(full_sequence) or not full_sequence or full_sequence.strip() == '':
            result['alignment_error'] = "Missing full_sequence"
            return result
        
        if pd.isna(human_ortholog_id) or not human_ortholog_id:
            result['alignment_error'] = "Missing human_ortholog_id"
            return result
        
        if pd.isna(human_ortholog_sequence) or not human_ortholog_sequence or human_ortholog_sequence.strip() == '':
            result['alignment_error'] = "Missing human_ortholog_sequence"
            return result
        
        if pd.isna(human_phosphoprotein_IRES) or not human_phosphoprotein_IRES:
            result['alignment_error'] = "Missing human_phosphoprotein_IRES"
            return result
        
        # Parse IRES positions
        human_positions = parse_position_list(str(human_phosphoprotein_IRES))
        if not human_positions:
            result['alignment_error'] = "No valid positions in human_phosphoprotein_IRES"
            return result
        
        # Call EBI API for alignment
        aligned_seqs = call_ebi_emboss_needle(
            protein_id, full_sequence,
            human_ortholog_id, human_ortholog_sequence,
            email
        )
        
        # Get aligned sequences
        mouse_aligned = aligned_seqs.get(protein_id, '')
        human_aligned = aligned_seqs.get(human_ortholog_id, '')
        
        if not mouse_aligned or not human_aligned:
            result['alignment_error'] = "Failed to extract aligned sequences from API response"
            return result
        
        # Calculate alignment identity
        identity = calculate_alignment_identity(human_aligned, mouse_aligned)
        result['alignment_percent'] = round(identity, 2)
        
        # Map positions
        mapping_result = map_positions_via_alignment(
            human_positions, human_aligned, mouse_aligned,
            human_ortholog_id, protein_id
        )
        
        # Format outputs
        result['aligned_IRES'] = format_ires_output(mapping_result['mapped_positions'])
        result['matched_IRES'] = format_ires_output(mapping_result['matched_positions'])
        result['mismatched_IRES'] = format_ires_output(mapping_result['mismatched_positions'])
        
        logger.info(f"Processed {protein_id}: {len(mapping_result['mapped_positions'])} positions mapped, "
                   f"{len(mapping_result['matched_positions'])} matches, {identity:.2f}% identity")
        
    except Exception as e:
        logger.error(f"Error processing row for {row.get('protein_id', 'unknown')}: {e}")
        result['alignment_error'] = str(e)
    
    return result


def main():
    """Main function to process int_interfaces.csv and add alignment columns."""
    parser = argparse.ArgumentParser(
        description='Align mouse/rat proteins with human orthologs and map IRES positions',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process mouse data
  python scripts/align_interfaces.py --input data/processed/mouse/int_interfaces.csv --species mouse
  
  # Process rat data
  python scripts/align_interfaces.py --input data/processed/rat/int_interfaces.csv --species rat
  
  # With custom email
  python scripts/align_interfaces.py --input data/processed/mouse/int_interfaces.csv --species mouse --email your@email.com
        """
    )
    
    parser.add_argument(
        '--input',
        type=str,
        required=True,
        help='Path to int_interfaces.csv file'
    )
    
    parser.add_argument(
        '--species',
        type=str,
        required=True,
        choices=['mouse', 'rat'],
        help='Species name (mouse or rat)'
    )
    
    parser.add_argument(
        '--email',
        type=str,
        default='user@example.com',
        help='Email address for EBI API (required by API)'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Output file path (default: {species}_aligned_interfaces.csv in same directory)'
    )
    
    args = parser.parse_args()
    
    # Check if input file exists
    input_path = Path(args.input)
    if not input_path.exists():
        logger.error(f"Input file not found: {args.input}")
        return 1
    
    # Determine output path
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.parent / "aligned_interfaces.csv"
    
    # Load input file
    logger.info(f"Loading input file: {args.input}")
    try:
        df = pd.read_csv(args.input)
        logger.info(f"Loaded {len(df)} rows")
    except Exception as e:
        logger.error(f"Failed to load input file: {e}")
        return 1
    
    # Check required columns
    required_cols = ['protein_id', 'full_sequence', 'human_ortholog_id', 
                     'human_ortholog_sequence', 'human_phosphoprotein_IRES']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        logger.error(f"Missing required columns: {missing_cols}")
        logger.error(f"Available columns: {list(df.columns)}")
        return 1
    
    # Initialize new columns
    df['aligned_IRES'] = None
    df['alignment_percent'] = None
    df['matched_IRES'] = None
    df['mismatched_IRES'] = None
    df['alignment_error'] = None
    
    # Count rows with human_phosphoprotein_IRES data
    rows_with_ires = df[
        df['human_phosphoprotein_IRES'].notna() &
        (df['human_phosphoprotein_IRES'] != '')
    ]
    
    logger.info(f"Found {len(rows_with_ires)} rows with human_phosphoprotein_IRES data to process")
    logger.info(f"Total rows in file: {len(df)}")
    
    if len(rows_with_ires) == 0:
        logger.warning("No rows with human_phosphoprotein_IRES data. Saving input file with empty new columns.")
        df.to_csv(output_path, index=False)
        logger.info(f"Saved to: {output_path}")
        return 0
    
    # Process each row - skip if no human_phosphoprotein_IRES
    logger.info("Starting alignment processing...")
    processed_count = 0
    skipped_count = 0
    
    for idx, row in df.iterrows():
        # Skip rows without human_phosphoprotein_IRES
        human_ires = row.get('human_phosphoprotein_IRES')
        if pd.isna(human_ires) or not human_ires or str(human_ires).strip() == '':
            skipped_count += 1
            continue
        
        processed_count += 1
        logger.info(f"Processing row {processed_count}/{len(rows_with_ires)} (index {idx}): {row.get('protein_id', 'unknown')}")
        
        result = process_row_alignment(row, args.species, args.email)
        
        # Update dataframe
        df.at[idx, 'aligned_IRES'] = result['aligned_IRES']
        df.at[idx, 'alignment_percent'] = result['alignment_percent']
        df.at[idx, 'matched_IRES'] = result['matched_IRES']
        df.at[idx, 'mismatched_IRES'] = result['mismatched_IRES']
        df.at[idx, 'alignment_error'] = result['alignment_error']
        
        # Small delay between rows for rate limiting
        time.sleep(0.5)
    
    logger.info(f"Skipped {skipped_count} rows without human_phosphoprotein_IRES data")
    
    # Save results
    logger.info(f"Saving results to: {output_path}")
    try:
        df.to_csv(output_path, index=False)
        logger.info(f"✓ Successfully saved {len(df)} rows")
    except Exception as e:
        logger.error(f"Failed to save output file: {e}")
        return 1
    
    # Print summary
    logger.info("")
    logger.info("=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    logger.info(f"Total rows processed: {len(df)}")
    logger.info(f"Rows with successful alignment: {df['alignment_error'].isna().sum()}")
    logger.info(f"Rows with errors: {df['alignment_error'].notna().sum()}")
    if df['alignment_percent'].notna().any():
        logger.info(f"Average alignment identity: {df['alignment_percent'].mean():.2f}%")
    logger.info(f"Output file: {output_path}")
    logger.info("")
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())

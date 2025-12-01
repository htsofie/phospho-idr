#!/usr/bin/env python3
"""
Compile MoRF data from CV1–CV4 FASTA files into a single CSV.

This script:
1. Parses CV1–CV4 .af files to extract:
   - ID
   - sequence
   - MoRFTrain-based residues (as [start-end,...])
   - database_type (uniprot / disprot / ideal)
2. Optionally maps DisProt/IDEAL IDs to UniProt using the UniProt ID Mapping API.
3. Saves the compiled dataset to a CSV that can be reused for both mouse and rat.

Usage examples:
    # Compile and map IDs, default output: data/morf_interfaces.csv
    python scripts/compile_morf_interfaces.py

    # Custom output and skip UniProt mapping
    python scripts/compile_morf_interfaces.py --output data/custom_morf_interfaces.csv --no-map
"""

import argparse
import logging
import os
from pathlib import Path
from typing import Optional, List, Dict

import pandas as pd
import requests
import time

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

# API rate limiting delay (seconds)
API_DELAY = 0.15


def extract_morf_positions(annotation_line: str, sequence: str) -> str:
    """
    Parse annotation line and convert '1' positions to range format.

    Args:
        annotation_line: String of 0, 1, - characters
        sequence: Protein sequence (for potential future validation)

    Returns:
        Formatted string: [20-30,40,50-60] (1-indexed positions)
    """
    positions: List[int] = []

    # Find all positions where character is '1' (1-indexed)
    for i, char in enumerate(annotation_line, start=1):
        if char == "1":
            positions.append(i)

    if not positions:
        return "[]"

    # Group consecutive positions into ranges
    ranges: List[str] = []
    start = positions[0]
    end = positions[0]

    for i in range(1, len(positions)):
        if positions[i] == end + 1:
            # Consecutive, extend range
            end = positions[i]
        else:
            # Gap found, save current range
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{end}")
            start = positions[i]
            end = positions[i]

    # Add last range
    if start == end:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}-{end}")

    return "[" + ",".join(ranges) + "]"


def classify_database_type(protein_id: str) -> str:
    """
    Classify ID as uniprot, disprot, or ideal.

    Args:
        protein_id: Protein ID string

    Returns:
        'disprot', 'ideal', or 'uniprot'
    """
    if protein_id.startswith("DP"):
        return "disprot"
    if protein_id.startswith("IID"):
        return "ideal"
    return "uniprot"


def parse_morf_fasta(fasta_file: str) -> pd.DataFrame:
    """
    Parse FASTA file and extract MoRF data.

    Format:
        >ID
        sequence
        MoRFTest annotation (skip)
        MoRFTrain annotation (use this - third line after ID)

    Args:
        fasta_file: Path to FASTA file

    Returns:
        DataFrame with columns: ID, sequence, morf_residues, database_type
    """
    logger.info(f"Parsing {fasta_file}...")

    data: List[Dict[str, str]] = []

    with open(fasta_file, "r") as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # Check for ID line
        if line.startswith(">"):
            protein_id = line[1:].strip()  # Remove '>'

            # Get sequence (next line)
            if i + 1 < len(lines):
                sequence = lines[i + 1].strip()
            else:
                logger.warning(f"Missing sequence for {protein_id}")
                i += 1
                continue

            # Get MoRFTrain annotation (third line after ID)
            if i + 3 < len(lines):
                morftrain_annotation = lines[i + 3].strip()
            else:
                logger.warning(f"Missing MoRFTrain annotation for {protein_id}")
                i += 1
                continue

            morf_residues = extract_morf_positions(morftrain_annotation, sequence)
            database_type = classify_database_type(protein_id)

            data.append(
                {
                    "ID": protein_id,
                    "sequence": sequence,
                    "morf_residues": morf_residues,
                    "database_type": database_type,
                }
            )

            i += 4  # Move to next entry
        else:
            i += 1

    df = pd.DataFrame(data)
    logger.info(f"  Extracted {len(df)} entries from {fasta_file}")
    return df


def compile_morf_data(fasta_files: List[str]) -> pd.DataFrame:
    """
    Process CV1–CV4 files and combine into single DataFrame.

    Args:
        fasta_files: List of paths to FASTA files

    Returns:
        Combined DataFrame
    """
    logger.info("Compiling MoRF data from FASTA files...")

    all_dfs: List[pd.DataFrame] = []
    for fasta_file in fasta_files:
        if os.path.exists(fasta_file):
            df = parse_morf_fasta(fasta_file)
            all_dfs.append(df)
        else:
            logger.warning(f"File not found: {fasta_file}")

    if not all_dfs:
        logger.error("No FASTA files could be parsed")
        return pd.DataFrame(columns=["ID", "sequence", "morf_residues", "database_type"])

    combined_df = pd.concat(all_dfs, ignore_index=True)
    logger.info(f"Total entries compiled: {len(combined_df)}")

    return combined_df


def map_id_to_uniprot(protein_id: str, from_db: str, cache: Dict[str, Optional[str]]) -> Optional[str]:
    """
    Map Disprot/Ideal ID to UniProt using UniProt ID Mapping API.

    Args:
        protein_id: Disprot or Ideal ID
        from_db: 'DisProt' or 'IDEAL'
        cache: Dictionary to cache results

    Returns:
        UniProt ID or None if not found
    """
    cache_key = f"{from_db}:{protein_id}"
    if cache_key in cache:
        return cache[cache_key]

    url = "https://rest.uniprot.org/idmapping/run"

    data = {
        "from": from_db,
        "to": "UniProtKB",
        "ids": protein_id,
    }

    try:
        response = requests.post(url, data=data, timeout=30)
        response.raise_for_status()

        job_id = response.json().get("jobId")
        if not job_id:
            logger.error(f"No job ID returned for {protein_id}")
            cache[cache_key] = None
            return None

        status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
        max_attempts = 30
        for attempt in range(max_attempts):
            time.sleep(API_DELAY)
            status_response = requests.get(status_url, timeout=30)
            status_response.raise_for_status()
            status_data = status_response.json()

            job_status = status_data.get("jobStatus")

            # Check if results are directly in the status response
            if "results" in status_data and len(status_data["results"]) > 0:
                result = status_data["results"][0]
                if "to" in result:
                    to_data = result["to"]
                    if isinstance(to_data, dict) and "primaryAccession" in to_data:
                        uniprot_id = to_data["primaryAccession"]
                        cache[cache_key] = uniprot_id
                        return uniprot_id
                    if isinstance(to_data, str):
                        uniprot_id = to_data
                        cache[cache_key] = uniprot_id
                        return uniprot_id

            if job_status == "FINISHED":
                # Fallback to stream endpoint
                results_url = f"https://rest.uniprot.org/idmapping/stream/{job_id}?format=tsv"
                results_response = requests.get(results_url, timeout=30)
                results_response.raise_for_status()

                lines = results_response.text.strip().split("\n")
                if len(lines) > 1:
                    parts = lines[1].split("\t")
                    if len(parts) > 1:
                        uniprot_id = parts[1].strip()
                        if uniprot_id:
                            cache[cache_key] = uniprot_id
                            return uniprot_id

                logger.warning(f"No mapping found for {protein_id}")
                cache[cache_key] = None
                return None
            if job_status == "ERROR":
                error_msg = status_data.get("error", status_data.get("messages", ["Unknown error"]))
                logger.error(f"Job error for {protein_id}: {error_msg}")
                cache[cache_key] = None
                return None
            if job_status in ["RUNNING", "NEW"]:
                continue
            if job_status is None:
                if attempt < max_attempts - 1:
                    time.sleep(1)
                    continue
                logger.warning(f"Job status is None for {protein_id} after {max_attempts} attempts")
                cache[cache_key] = None
                return None

            logger.warning(f"Unexpected job status for {protein_id}: {job_status}")
            cache[cache_key] = None
            return None

        logger.warning(f"Timeout waiting for mapping result for {protein_id}")
        cache[cache_key] = None
        return None

    except requests.exceptions.RequestException as e:
        logger.error(f"Error mapping {protein_id} from {from_db}: {e}")
        cache[cache_key] = None
        return None


def map_all_ids_to_uniprot(morf_df: pd.DataFrame) -> pd.DataFrame:
    """
    Map all Disprot/Ideal IDs to UniProt and add uniprot_id / mapping_error columns.
    """
    logger.info("Mapping Disprot/Ideal IDs to UniProt...")

    morf_df = morf_df.copy()
    morf_df["uniprot_id"] = None
    morf_df["mapping_error"] = None

    cache: Dict[str, Optional[str]] = {}

    to_map = morf_df[morf_df["database_type"].isin(["disprot", "ideal"])]
    logger.info(f"  Found {len(to_map)} IDs to map")

    for idx, row in to_map.iterrows():
        protein_id = row["ID"]
        db_type = row["database_type"]
        from_db = "DisProt" if db_type == "disprot" else "IDEAL"

        logger.info(f"  Mapping {protein_id} ({from_db})...")
        uniprot_id = map_id_to_uniprot(protein_id, from_db, cache)
        time.sleep(API_DELAY)

        if uniprot_id:
            morf_df.at[idx, "uniprot_id"] = uniprot_id
            logger.info(f"    -> {uniprot_id}")
        else:
            morf_df.at[idx, "mapping_error"] = f"Failed to map {protein_id} from {from_db}"
            logger.warning("    -> Mapping failed")

    # For uniprot type, use original ID
    uniprot_mask = morf_df["database_type"] == "uniprot"
    morf_df.loc[uniprot_mask, "uniprot_id"] = morf_df.loc[uniprot_mask, "ID"]

    logger.info(f"  Mapping complete. {len(morf_df[morf_df['uniprot_id'].notna()])} IDs have UniProt mappings")
    return morf_df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compile MoRF data from CV1–CV4 into a reusable CSV."
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to output CSV (default: data/morf_interfaces.csv)",
    )
    parser.add_argument(
        "--no-map",
        action="store_true",
        help="Skip mapping DisProt/IDEAL IDs to UniProt (no uniprot_id column).",
    )

    args = parser.parse_args()

    # Determine project root and MoRF data directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    morf_dir = os.path.join(project_root, "data", "MoRF2")

    fasta_files = [
        os.path.join(morf_dir, "CV1.af"),
        os.path.join(morf_dir, "CV2.af"),
        os.path.join(morf_dir, "CV3.af"),
        os.path.join(morf_dir, "CV4.af"),
    ]

    # Determine output path
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = Path(project_root) / "data" / "morf_interfaces.csv"

    logger.info("=" * 80)
    logger.info("Compile MoRF Interfaces")
    logger.info("=" * 80)
    logger.info(f"Output file: {output_path}")
    logger.info(f"Map DisProt/IDEAL IDs to UniProt: {not args.no_map}")
    logger.info("")

    # Step 1: Compile MoRF data
    morf_df = compile_morf_data(fasta_files)
    if len(morf_df) == 0:
        logger.error("No MoRF data extracted. Exiting.")
        return

    # Step 2: Map IDs (optional)
    if not args.no_map:
        morf_df = map_all_ids_to_uniprot(morf_df)

    # Save compiled MoRF data
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Saving compiled MoRF data to {output_path}...")
    morf_df.to_csv(output_path, index=False)
    logger.info(f"✓ Successfully saved {len(morf_df)} rows")

    logger.info("")
    logger.info("=" * 80)
    logger.info("Summary")
    logger.info("=" * 80)
    logger.info(f"Total MoRF entries: {len(morf_df)}")
    logger.info(f"  - UniProt IDs: {len(morf_df[morf_df['database_type'] == 'uniprot'])}")
    logger.info(f"  - Disprot IDs: {len(morf_df[morf_df['database_type'] == 'disprot'])}")
    logger.info(f"  - Ideal IDs: {len(morf_df[morf_df['database_type'] == 'ideal'])}")
    if "uniprot_id" in morf_df.columns:
        logger.info(f"Successfully mapped IDs: {len(morf_df[morf_df['uniprot_id'].notna()])}")
    logger.info(f"Output file: {output_path}")
    logger.info("")


if __name__ == "__main__":
    main()



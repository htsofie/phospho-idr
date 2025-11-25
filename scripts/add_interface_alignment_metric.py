#!/usr/bin/env python3
"""
Utility to rename alignment columns and compute localized interface alignment percentage.

Usage:
    python scripts/add_interface_alignment_metric.py \
        --input data/processed/mouse/aligned_interfaces.csv \
        --output data/processed/mouse/aligned_interfaces_with_interface_pct.csv
"""

import argparse
import logging
import os
from typing import List

import pandas as pd

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def parse_position_list(position_str: str) -> List[int]:
    """Convert a string like "[33-35,60,62]" into a list of integers."""
    if position_str is None or (isinstance(position_str, float) and pd.isna(position_str)):
        return []

    position_str = str(position_str).strip()
    if not position_str:
        return []

    # Remove enclosing brackets if present
    if position_str.startswith('[') and position_str.endswith(']'):
        position_str = position_str[1:-1]

    if not position_str:
        return []

    positions: List[int] = []
    for part in position_str.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            try:
                start, end = part.split('-', 1)
                start = int(start)
                end = int(end)
                if start <= end:
                    positions.extend(range(start, end + 1))
                else:
                    positions.extend(range(end, start + 1))
            except ValueError:
                logger.debug("Unable to parse range '%s'", part)
        else:
            try:
                positions.append(int(part))
            except ValueError:
                logger.debug("Unable to parse position '%s'", part)

    return positions


def compute_interface_alignment(df: pd.DataFrame, matched_col: str, total_col: str,
                                 output_col: str = 'interface_alignment_percent') -> pd.DataFrame:
    """Compute interface alignment percentage from matched and total interface columns."""
    matched_counts = df[matched_col].apply(lambda x: len(parse_position_list(x)))
    total_counts = df[total_col].apply(lambda x: len(parse_position_list(x)))

    def calc_percent(matched: int, total: int) -> float:
        if total == 0:
            return None
        return round((matched / total) * 100.0, 2)

    df[output_col] = [calc_percent(m, t) for m, t in zip(matched_counts, total_counts)]
    df[f'{matched_col}_count'] = matched_counts
    df[f'{total_col}_count'] = total_counts
    return df


def main():
    parser = argparse.ArgumentParser(
        description='Rename alignment column and compute interface alignment percentage'
    )
    parser.add_argument('--input', required=True, help='Path to aligned_interfaces CSV')
    parser.add_argument('--output', help='Path to save updated CSV (default: overwrite input)')
    parser.add_argument('--matched-col', default='matched_IRES',
                        help='Column containing matched interface residues (default: matched_IRES)')
    parser.add_argument('--total-col', default='human_phosphoprotein_IRES',
                        help='Column containing total interface residues (default: human_phosphoprotein_IRES)')
    parser.add_argument('--rename-only', action='store_true',
                        help='Only rename alignment column without computing interface percent')

    args = parser.parse_args()

    if not os.path.exists(args.input):
        logger.error("Input file not found: %s", args.input)
        return 1

    output_path = args.output if args.output else args.input

    logger.info("Loading %s", args.input)
    try:
        df = pd.read_csv(args.input)
    except Exception as exc:
        logger.error("Failed to read CSV: %s", exc)
        return 1

    required_cols = ['alignment_percent']
    for col in [args.matched_col, args.total_col]:
        if not args.rename_only:
            required_cols.append(col)

    missing = [col for col in required_cols if col not in df.columns]
    if missing:
        logger.error("Missing required columns: %s", missing)
        return 1

    if 'total_alignment' in df.columns:
        logger.warning("Column 'total_alignment' already exists. It will be overwritten.")

    df = df.rename(columns={'alignment_percent': 'total_alignment'})

    if not args.rename_only:
        logger.info("Computing interface alignment percent using '%s' and '%s'",
                    args.matched_col, args.total_col)
        df = compute_interface_alignment(df, args.matched_col, args.total_col)
    else:
        logger.info("Rename-only mode: skipping interface percent computation")

    logger.info("Saving results to %s", output_path)
    try:
        df.to_csv(output_path, index=False)
    except Exception as exc:
        logger.error("Failed to save CSV: %s", exc)
        return 1

    logger.info("Completed. Total rows: %d", len(df))
    if not args.rename_only:
        non_null = df['interface_alignment_percent'].notna().sum()
        logger.info("Rows with interface alignment percent: %d", non_null)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

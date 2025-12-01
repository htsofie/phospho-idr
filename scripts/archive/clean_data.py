#!/usr/bin/env python3
import argparse
import os
import re
from typing import Dict, Optional, Tuple

import pandas as pd
import yaml

# python scripts/clean_data.py --input data/processed/mouse/test_data.csv --output data/processed/mouse/cleaned_test_data.csv --species mouse
# # Original config-based approach
# python scripts/clean_data.py --config configs/mouse.yaml --mode sample --output-format csv
# python scripts/clean_data.py --config configs/rat.yaml --mode full --output-format parquet


def load_config(config_path: str) -> Dict:
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)


def infer_species_from_config_path(config_path: str) -> str:
    base_name = os.path.basename(config_path)
    species = os.path.splitext(base_name)[0]
    return species


def get_processed_paths(species: str, mode: str) -> Tuple[str, Optional[str]]:
    processed_dir = os.path.join('data', 'processed', species)
    base = 'test_data' if mode == 'sample' else 'full_data'
    parquet_path = os.path.join(processed_dir, f'{base}.parquet')
    csv_path = os.path.join(processed_dir, f'{base}.csv')
    # Prefer parquet if present, else csv. Return both for logging.
    if os.path.exists(parquet_path):
        return parquet_path, csv_path if os.path.exists(csv_path) else None
    if os.path.exists(csv_path):
        return csv_path, None
    return parquet_path, csv_path  # non-existent hints for error message


def read_processed(path: str) -> pd.DataFrame:
    if path.endswith('.parquet'):
        return pd.read_parquet(path)
    if path.endswith('.csv'):
        return pd.read_csv(path)
    raise ValueError(f'Unsupported input format for {path}')


# -------------------- Rat-specific logic --------------------

def filter_ids_rat(df: pd.DataFrame, config: Dict) -> pd.DataFrame:
    # Don't remove rows based on missing IDs - keep all rows
    # Only apply localization threshold filtering
    return df


def apply_localization_threshold(df: pd.DataFrame, config: Dict) -> pd.DataFrame:
    loc_cfg = config.get('filters', {}).get('localization', {})
    if not loc_cfg:
        return df
    if loc_cfg.get('type') != 'numerical':
        return df
    column_name = loc_cfg.get('column')
    threshold = loc_cfg.get('threshold')
    if column_name not in df.columns or threshold is None:
        return df
    return df.loc[pd.to_numeric(df[column_name], errors='coerce') >= float(threshold)].copy()


# -------------------- Mouse-specific logic --------------------

def extract_ids_mouse(df: pd.DataFrame, config: Dict) -> pd.DataFrame:
    # Look for the standardized protein_IDs column
    source_col = 'protein_IDs'
    
    if source_col not in df.columns:
        # If not found, skip extraction gracefully
        df['swissprot_id'] = None
        df['trembl_id'] = None
        df['ensembl'] = None
        return df

    def parse_ids(cell: object) -> pd.Series:
        if pd.isna(cell):
            return pd.Series({'swissprot_id': None, 'trembl_id': None, 'ensembl': None})
        # Cell is expected like: "IPI:IPI00118763.1|SWISS-PROT:O54826|TREMBL:Q543D7;Q8VDP9|ENSEMBL:ENSMUSP000000" etc.
        try:
            parts = str(cell).split('|')
        except Exception:
            parts = [str(cell)]
        swissprot_val: Optional[str] = None
        trembl_val: Optional[str] = None
        ensembl_val: Optional[str] = None
        for part in parts:
            part_stripped = str(part).strip()
            if part_stripped.startswith('SWISS-PROT:'):
                swissprot_val = part_stripped.split(':', 1)[1].strip() or None
            elif part_stripped.startswith('TREMBL:'):
                trembl_val = part_stripped.split(':', 1)[1].strip() or None
            elif part_stripped.startswith('ENSEMBL:'):
                ensembl_id = part_stripped.split(':', 1)[1].strip() or None
                # Filter out incomplete ENSEMBL IDs - this appears to be a data quality issue
                # in the original dataset where some ENSEMBL IDs are truncated (e.g., "ENSMUS" instead of "ENSMUSP00000012345")
                # We only keep complete ENSEMBL IDs that match the expected pattern
                if ensembl_id and len(ensembl_id) >= 15 and (ensembl_id.startswith('ENSMUSP') or ensembl_id.startswith('ENSMUSG')):
                    ensembl_val = ensembl_id
                # If ENSEMBL ID is incomplete, we set it to None but keep the row
                # (other identifiers like SWISS-PROT or TREMBL may still be valid)
        
        return pd.Series({'swissprot_id': swissprot_val, 'trembl_id': trembl_val, 'ensembl': ensembl_val})

    ids = df[source_col].apply(parse_ids)
    
    # Don't remove rows based on missing IDs - keep all rows
    # Just extract the IDs and keep all rows
    ids_filtered = ids
    df_filtered = df.drop(columns=[source_col])
    
    # Concatenate ID columns first, then other columns
    df_out = pd.concat([ids_filtered, df_filtered], axis=1)
    
    return df_out


def apply_mouse_categorical_filters(df: pd.DataFrame, config: Dict) -> pd.DataFrame:
    loc_cfg = config.get('filters', {}).get('localization', {})
    if not loc_cfg:
        return df
    if loc_cfg.get('type') != 'categorical':
        return df
    column_name = loc_cfg.get('column')
    exclude_values = set(loc_cfg.get('exclude', []))
    if column_name not in df.columns:
        return df
    return df.loc[~df[column_name].isin(exclude_values)].copy()


# -------------------- Modified sequence processing --------------------

def clean_sequence_motif(motif_seq: str, species: str) -> Tuple[str, Optional[object]]:
    """
    Clean sequence motif(s) and extract phosphorylation site position(s).

    - Supports single motifs or multiple motifs separated by ';'.
    - Each segment is cleaned independently, then rejoined by ';'.
    
    Args:
        motif_seq: Raw sequence motif string, possibly containing ';' separators
        species: 'rat' or 'mouse' to determine cleaning rules
        
    Returns:
        Tuple of (cleaned_motif, motif_position_or_positions)
        motif_position_or_positions is either an int (single motif) or a ';'-joined string of ints (multiple motifs),
        representing the 1-based index of phosphorylation site in each cleaned motif.
    """
    if pd.isna(motif_seq) or not motif_seq:
        return "", None

    raw = str(motif_seq).strip()

    # Split on ';' to support multiple sequence windows
    segments = [seg.strip() for seg in raw.split(';')] if ';' in raw else [raw]

    cleaned_segments: list[str] = []
    positions: list[int] = []

    for seg in segments:
        seq = seg
        if not seq:
            cleaned_segments.append("")
            positions.append(None)  # placeholder to keep alignment
            continue

        if species == 'rat':
            # Rat: site is at position 7; remove underscores, dashes and adjust for underscores before position 7
            underscores_before_7 = seq[:6].count('_')
            seq = seq.replace('_', '')
            seq = seq.replace('-', '')
            motif_pos = 7 - underscores_before_7
        elif species == 'mouse':
            # Mouse: '*' marks site at position 7; remove '*', underscores, and dashes, adjust for underscores before 7
            underscores_before_7 = seq[:6].count('_')
            seq = seq.replace('*', '')
            seq = seq.replace('_', '')
            seq = seq.replace('-', '')
            motif_pos = 7 - underscores_before_7
        else:
            # Generic cleaning
            underscores_before_7 = seq[:6].count('_')
            seq = seq.replace('*', '')
            seq = seq.replace('_', '')
            seq = seq.replace('-', '')
            motif_pos = 7 - underscores_before_7

        cleaned_segments.append(seq)
        positions.append(motif_pos)

    if len(cleaned_segments) == 1:
        return cleaned_segments[0], positions[0]
    # Join multiple segments and their corresponding positions with ';'
    cleaned_joined = ';'.join(cleaned_segments)
    pos_joined = ';'.join(str(p) if p is not None else '' for p in positions)
    return cleaned_joined, pos_joined


def process_sequence_motifs(df: pd.DataFrame, species: str) -> pd.DataFrame:
    """Process sequence motifs for both rat and mouse data.

    - Supports semicolon-separated multiple sequence windows per row. Each is cleaned independently
      and rejoined with ';'. The corresponding positions are also joined by ';'.
    """
    if 'site_motif' not in df.columns:
        print(f"Warning: 'site_motif' column not found in dataset. Available columns: {list(df.columns)}")
        return df
    
    # Apply cleaning to each row
    results = df['site_motif'].apply(lambda x: clean_sequence_motif(x, species))
    
    # Split results into separate columns
    cleaned_motifs, motif_positions = zip(*results)
    
    # Add new columns
    df = df.copy()
    df['cleaned_site_motif'] = cleaned_motifs
    df['motif_position'] = motif_positions
    
    # Reorder columns to put cleaned_site_motif right after site_motif
    cols = list(df.columns)
    motif_idx = cols.index('site_motif')
    
    # Remove the new columns from their current positions
    cols.remove('cleaned_site_motif')
    cols.remove('motif_position')
    
    # Insert them right after site_motif
    cols.insert(motif_idx + 1, 'cleaned_site_motif')
    cols.insert(motif_idx + 2, 'motif_position')
    
    return df[cols]


# -------------------- Standardization --------------------

def standardize_id_columns(df: pd.DataFrame, species: str, config: Dict) -> pd.DataFrame:
    df_out = df.copy()
    if species == 'rat':
        # Rename rat ID columns to standardized names if present
        if 'uniprot_id' in df_out.columns:
            df_out = df_out.rename(columns={'uniprot_id': 'uniprot'})
        if 'ensembl_id' in df_out.columns:
            df_out = df_out.rename(columns={'ensembl_id': 'ensembl'})
    # For mouse, IDs are already in 'swissprot_id', 'trembl_id', and 'ensembl' columns
    return df_out


def detect_species(config: Dict, config_path: str) -> str:
    # Prefer to infer from config filename, else from raw_data path parent
    species = infer_species_from_config_path(config_path)
    if species in {'mouse', 'rat'}:
        return species
    # Fallback to parent folder of raw_data
    raw_path = config.get('raw_data')
    if raw_path:
        candidate = os.path.basename(os.path.dirname(raw_path))
        if candidate in {'mouse', 'rat'}:
            return candidate
    # Default to species token from filename even if unknown
    return species


def write_output(df: pd.DataFrame, species: str, mode: str, output_format: str) -> str:
    output_dir = os.path.join('data', 'processed', species)
    os.makedirs(output_dir, exist_ok=True)
    out_base = f'cleaned_{mode}'
    out_path = os.path.join(output_dir, f'{out_base}.{"parquet" if output_format == "parquet" else "csv"}')
    if output_format == 'parquet':
        try:
            df.to_parquet(out_path, index=False)
        except Exception as exc:
            # Fallback to CSV
            out_path = os.path.join(output_dir, f'{out_base}.csv')
            df.to_csv(out_path, index=False)
    else:
        df.to_csv(out_path, index=False)
    return out_path


def main() -> int:
    parser = argparse.ArgumentParser(description='Clean and filter processed phosphorylation datasets')
    
    # Input options - either config-based or direct file paths
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--config', '-c', help='Path to species YAML configuration file')
    input_group.add_argument('--input', '-i', help='Path to input CSV file')
    
    # Output options
    parser.add_argument('--output', '-o', help='Path to output CSV file (required when using --input)')
    parser.add_argument('--species', '-s', choices=['mouse', 'rat'], help='Species name (required when using --input)')
    
    # Legacy options for config-based processing
    parser.add_argument('--mode', '-m', choices=['sample', 'full'], default='sample', help='Input mode: sample (CSV) or full (Parquet)')
    parser.add_argument('--output-format', '-f', choices=['csv', 'parquet'], default='parquet', help='Output file format')
    
    args = parser.parse_args()

    # Determine processing mode
    if args.config:
        # Config-based processing (legacy mode)
        config = load_config(args.config)
        species = detect_species(config, args.config)
        input_path, _ = get_processed_paths(species, args.mode)
        output_path = None  # Will be determined by write_output
    else:
        # Direct file processing
        if not args.input or not args.output or not args.species:
            parser.error("When using --input, --output and --species are required")
        
        species = args.species
        input_path = args.input
        output_path = args.output
        
        # Load default config for the species
        config_path = f'configs/{species}.yaml'
        if not os.path.exists(config_path):
            parser.error(f"Config file not found: {config_path}")
        config = load_config(config_path)

    # Check if input file exists
    if not os.path.exists(input_path):
        print(f'Input file not found: {input_path}')
        return 1

    # Read the input file
    if input_path.endswith('.parquet'):
        df = pd.read_parquet(input_path)
    elif input_path.endswith('.csv'):
        df = pd.read_csv(input_path)
    else:
        print(f'Unsupported file format: {input_path}')
        return 1

    # Species-specific processing
    if species == 'rat':
        df = apply_localization_threshold(df, config)
        df = filter_ids_rat(df, config)
        df = standardize_id_columns(df, species, config)
        df = process_sequence_motifs(df, species)
    elif species == 'mouse':
        df = apply_mouse_categorical_filters(df, config)
        df = extract_ids_mouse(df, config)
        df = standardize_id_columns(df, species, config)
        df = process_sequence_motifs(df, species)
        # Ensure 'Protein' is first and 'protein_description' is second if present
        if 'Protein' in df.columns:
            ordered_cols = ['Protein']
            if 'protein_description' in df.columns:
                ordered_cols.append('protein_description')
            other_cols = [c for c in df.columns if c not in ordered_cols]
            df = df[ordered_cols + other_cols]
    else:
        # If unknown species, attempt generic filters if present
        df = apply_mouse_categorical_filters(df, config)
        df = apply_localization_threshold(df, config)
        df = filter_ids_rat(df, config)
        df = standardize_id_columns(df, species, config)
        df = process_sequence_motifs(df, species)

    # Save output
    if output_path:
        # Direct output path specified
        saved_to = output_path
        df.to_csv(output_path, index=False)
    else:
        # Use legacy write_output function
        saved_to = write_output(df, species, args.mode, args.output_format)

    # Brief report
    print(f'species: {species}')
    print(f'input: {input_path}')
    print(f'rows: {len(df)}')
    print(f'saved_to: {saved_to}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

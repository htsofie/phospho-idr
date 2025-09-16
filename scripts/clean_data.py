#!/usr/bin/env python3
import argparse
import os
from typing import Dict, Optional, Tuple

import pandas as pd
import yaml

# To run in command line: python scripts/clean_data.py --config configs/rat.yaml --mode sample --output-format csv
# specify: --mode sample or --mode full, --output-format csv or --output-format parquet

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
    required_cols = config.get('filters', {}).get('required_ids', [])
    if not required_cols:
        return df
    present = [c for c in required_cols if c in df.columns]
    if not present:
        return df
    mask = df[present].notna().any(axis=1)
    return df.loc[mask].copy()


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
    
    # Remove rows that have no valid identifiers (all three ID types are None)
    # This ensures we only keep rows with at least one valid protein identifier
    valid_rows = ids.notna().any(axis=1)
    ids_filtered = ids[valid_rows]
    df_filtered = df[valid_rows].drop(columns=[source_col])
    
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
    parser.add_argument('--config', '-c', required=True, help='Path to species YAML configuration file')
    parser.add_argument('--mode', '-m', choices=['sample', 'full'], default='sample', help='Input mode: sample (CSV) or full (Parquet)')
    parser.add_argument('--output-format', '-f', choices=['csv', 'parquet'], default='parquet', help='Output file format')
    args = parser.parse_args()

    config = load_config(args.config)
    species = detect_species(config, args.config)

    input_path, alt_path = get_processed_paths(species, args.mode)
    if not os.path.exists(input_path):
        msg = f'Input not found: {input_path}'
        if alt_path is not None:
            msg += f' and {alt_path}'
        raise FileNotFoundError(msg)

    df = read_processed(input_path)

    # Species-specific processing
    if species == 'rat':
        df = apply_localization_threshold(df, config)
        df = filter_ids_rat(df, config)
        df = standardize_id_columns(df, species, config)
    elif species == 'mouse':
        df = apply_mouse_categorical_filters(df, config)
        df = extract_ids_mouse(df, config)
        df = standardize_id_columns(df, species, config)
    else:
        # If unknown species, attempt generic filters if present
        df = apply_mouse_categorical_filters(df, config)
        df = apply_localization_threshold(df, config)
        df = filter_ids_rat(df, config)
        df = standardize_id_columns(df, species, config)

    saved_to = write_output(df, species, args.mode, args.output_format)

    # Brief report
    print(f'species: {species}')
    print(f'input: {input_path}')
    print(f'rows: {len(df)}')
    print(f'saved_to: {saved_to}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

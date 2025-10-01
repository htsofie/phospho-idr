#!/usr/bin/env python3
import pandas as pd
import yaml
import argparse
import os


def load_config(config_path):
    #Load configuration from YAML file.
    with open(config_path, 'r') as file:
        return yaml.safe_load(file)


def main():
   # Creates command line argument parser, supporting sample or full processing and output format
   # Examples:python scripts/gen_data.py --config configs/rat.yaml --mode sample --sample-size 30 --random-state 42 --output-format csv
   #   Sample: 
   #   Full:   python scripts/gen_data.py --config configs/rat.yaml --mode full --output-format parquet
    parser = argparse.ArgumentParser(description='Generate datasets from phosphorylation data')
    parser.add_argument('--config', '-c', required=True, help='Path to configuration YAML file')
    parser.add_argument('--mode', '-m', choices=['sample', 'full'], default='sample', help='Process a random sample or the full dataset')
    parser.add_argument('--sample-size', '-n', type=int, default=30, help='Number of proteins to sample (when mode=sample)')
    parser.add_argument('--random-state', '-r', type=int, default=42, help='Random state for reproducibility (when mode=sample)')
    parser.add_argument('--output-format', '-f', choices=['csv', 'parquet'], default='csv', help='Output file format')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Load raw data
    raw_data_path = config['raw_data']
    if raw_data_path.endswith('.xlsb'):
        df = pd.read_excel(raw_data_path, engine='pyxlsb')
    else:
        df = pd.read_excel(raw_data_path)

    # Infer species from raw data path (parent directory name)
    species_or_group = os.path.basename(os.path.dirname(raw_data_path))

    # Helper to extract IPI for grouping and output column 'Protein'
    def extract_mouse_ipi(cell: object) -> str:
        # Expected like: "IPI:IPI00121135.5|SWISS-PROT:Q62093|TREMBL:A2AA29;Q99MY4;Q99MY5|ENSEMBL:ENSMUSP..."
        if pd.isna(cell):
            return None
        try:
            parts = str(cell).split('|')
        except Exception:
            parts = [str(cell)]
        for part in parts:
            p = str(part).strip()
            if p.startswith('IPI:'):
                val = p.split(':', 1)[1].strip()
                return val or None
        return None

    # Build a grouping/protein column before sampling
    group_col = '_protein_group'
    if species_or_group == 'mouse':
        source_col = 'Protein IPI (Assigned by Sequest)'
        if source_col in df.columns:
            df[group_col] = df[source_col].apply(extract_mouse_ipi)
        else:
            df[group_col] = None
    elif species_or_group == 'rat':
        # Rat dataset is expected to have a 'Protein' column containing IPI
        df[group_col] = df['Protein'] if 'Protein' in df.columns else None
    else:
        df[group_col] = None
    
    # Determine input frame depending on mode
    if args.mode == 'sample':
        # Sample by protein groups: choose random proteins, include all rows for those proteins
        if group_col in df.columns and df[group_col].notna().any():
            unique_proteins = (
                pd.Series(df[group_col].dropna().unique())
                .sample(n=min(args.sample_size, df[group_col].dropna().nunique()), random_state=args.random_state)
                .tolist()
            )
            df_input = df[df[group_col].isin(unique_proteins)].copy()
        else:
            # Fallback to row sampling if grouping info is unavailable
            df_input = df.sample(n=min(args.sample_size, len(df)), random_state=args.random_state)
        output_basename = 'test_data'
    else:
        df_input = df
        output_basename = 'full_data'
    
    # Filter and rename columns using config
    columns_mapping = config['columns_to_keep']
    available_columns = [col for col in columns_mapping.values() if col in df_input.columns]
    
    if not available_columns:
        print("Error: No valid columns found")
        return 1
    
    # Select only available columns and rename them to standardized names
    df_filtered = df_input[available_columns].copy()
    
    # Create reverse mapping: original_name -> standardized_name
    reverse_mapping = {v: k for k, v in columns_mapping.items() if v in df_input.columns}
    df_filtered = df_filtered.rename(columns=reverse_mapping)

    # Ensure 'Protein' column is present in output
    # - For mouse, derive from standardized 'protein_IDs' (if present) or precomputed group_col
    # - For rat, prefer existing 'Protein' from raw; otherwise use group_col if available
    if 'Protein' not in df_filtered.columns:
        protein_values = None
        if species_or_group == 'mouse':
            if 'protein_IDs' in df_filtered.columns:
                protein_values = df_filtered['protein_IDs'].apply(extract_mouse_ipi)
            elif group_col in df_input.columns:
                protein_values = df_input.loc[df_filtered.index, group_col]
        elif species_or_group == 'rat':
            if 'Protein' in df_input.columns:
                protein_values = df_input.loc[df_filtered.index, 'Protein']
            elif group_col in df_input.columns:
                protein_values = df_input.loc[df_filtered.index, group_col]
        if protein_values is None:
            protein_values = pd.Series([None] * len(df_filtered), index=df_filtered.index)
        df_filtered['Protein'] = protein_values
    
    # Move 'Protein' to be the first column and 'protein_description' second if present
    if 'Protein' in df_filtered.columns:
        # Determine desired order
        ordered_cols = ['Protein']
        if 'protein_description' in df_filtered.columns:
            ordered_cols.append('protein_description')
        other_cols = [c for c in df_filtered.columns if c not in ordered_cols]
        df_filtered = df_filtered[ordered_cols + other_cols]

    # Build output path and write
    output_dir = os.path.join('data', 'processed', species_or_group)
    os.makedirs(output_dir, exist_ok=True)
    
    if args.output_format == 'csv':
        processed_path = os.path.join(output_dir, f"{output_basename}.csv")
        df_filtered.to_csv(processed_path, index=False)
    else:
        processed_path = os.path.join(output_dir, f"{output_basename}.parquet")
        try:
            df_filtered.to_parquet(processed_path, index=False)
        except Exception as exc:
            print(f"Parquet write failed ({exc}). Falling back to CSV.")
            processed_path = os.path.join(output_dir, f"{output_basename}.csv")
            df_filtered.to_csv(processed_path, index=False)
    
    # Print results
    print(df_filtered.columns.tolist())
    print("shape:", df_filtered.shape)
    print(f"saved_to: {processed_path}")
    
    return 0


if __name__ == "__main__":
    exit(main())
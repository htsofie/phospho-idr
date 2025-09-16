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
   # Examples:
   #   Sample: python scripts/gen_data.py --config configs/rat.yaml --mode sample --sample-size 30 --random-state 42 --output-format csv
   #   Full:   python scripts/gen_data.py --config configs/rat.yaml --mode full --output-format parquet
    parser = argparse.ArgumentParser(description='Generate datasets from phosphorylation data')
    parser.add_argument('--config', '-c', required=True, help='Path to configuration YAML file')
    parser.add_argument('--mode', '-m', choices=['sample', 'full'], default='sample', help='Process a random sample or the full dataset')
    parser.add_argument('--sample-size', '-n', type=int, default=30, help='Number of samples to generate (when mode=sample)')
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
    
    # Determine input frame depending on mode
    if args.mode == 'sample':
        df_input = df.sample(n=args.sample_size, random_state=args.random_state)
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
    
    # Build output path and write
    species_or_group = os.path.basename(os.path.dirname(raw_data_path))
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
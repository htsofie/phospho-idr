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
   # Creates command line argument parser, with three arguments: --config (required), --sample-size (optionsal, default 20), --random-state (optional, default 40)
   # Parses command line and stores values in args
   # Example usage: python3 gen_test_data.py --config configs/rat.yaml --sample-size 30 --random-state 42
    parser = argparse.ArgumentParser(description='Generate test datasets from phosphorylation data')
    parser.add_argument('--config', '-c', required=True, help='Path to configuration YAML file')
    parser.add_argument('--sample-size', '-n', type=int, default=20, help='Number of samples to generate')
    parser.add_argument('--random-state', '-r', type=int, default=42, help='Random state for reproducibility')
    
    args = parser.parse_args()
    
    # Load configuration
    config = load_config(args.config)
    
    # Load raw data
    raw_data_path = config['raw_data']
    if raw_data_path.endswith('.xlsb'):
        df = pd.read_excel(raw_data_path, engine='pyxlsb')
    else:
        df = pd.read_excel(raw_data_path)
    
    # Sample data
    test_df = df.sample(n=args.sample_size, random_state=args.random_state)
    
    # Save to processed location
    processed_path = f"../data/processed/{os.path.basename(os.path.dirname(raw_data_path))}/test_data.csv"
    os.makedirs(os.path.dirname(processed_path), exist_ok=True)
    test_df.to_csv(processed_path, index=False)
    
    # Read back and filter by columns
    df_processed = pd.read_csv(processed_path)
    
    # Filter by columns specified in config
    columns_to_keep = list(config['columns_to_keep'].values())
    available_columns = [col for col in columns_to_keep if col in df_processed.columns]
    
    if not available_columns:
        print("Error: No valid columns found")
        return 1
    
    df_filtered = df_processed[available_columns]
    df_filtered.to_csv(processed_path, index=False)
    
    # Print results
    print(df_filtered.columns.tolist())
    print("shape:", df_filtered.shape)
    
    return 0


if __name__ == "__main__":
    exit(main())
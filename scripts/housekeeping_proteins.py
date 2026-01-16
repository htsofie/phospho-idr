#!/usr/bin/env python3
"""
Identify housekeeping proteins (expressed across all tissue types) for mouse and rat.

This script:
1. Filters proteins expressed in all specified tissue types
2. Calculates normalized TSC values for mouse
3. Calculates coefficient of variation (CV) for both species

Usage:
    python scripts/housekeeping_proteins.py --species mouse
    python scripts/housekeeping_proteins.py --species rat
"""

import argparse
import logging
import os
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def filter_housekeeping_proteins(df: pd.DataFrame, tissue_columns: List[str]) -> pd.DataFrame:
    """
    Filter rows where all specified tissue columns have values > 0.
    
    Args:
        df: Input DataFrame
        tissue_columns: List of column names to check
        
    Returns:
        Filtered DataFrame with only housekeeping proteins
    """
    logger.info(f"Filtering for proteins expressed in all {len(tissue_columns)} tissues...")
    
    # Check that all required columns exist
    missing_cols = [col for col in tissue_columns if col not in df.columns]
    if missing_cols:
        logger.error(f"Missing required columns: {missing_cols}")
        logger.error(f"Available columns: {list(df.columns)}")
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Filter rows where all tissue columns are > 0
    mask = True
    for col in tissue_columns:
        # Convert to numeric, replacing non-numeric values with NaN
        df[col] = pd.to_numeric(df[col], errors='coerce')
        mask = mask & (df[col] > 0)
    
    filtered_df = df[mask].copy()
    logger.info(f"Found {len(filtered_df)} housekeeping proteins (from {len(df)} total proteins)")
    
    return filtered_df


def calculate_normalized_tsc(df: pd.DataFrame, tsc_columns: List[str]) -> pd.DataFrame:
    """
    Calculate normalized TSC values for mouse data.
    
    Formula: ((TSC per tissue) / (sum of all TSC)) * 1,000,000
    
    Args:
        df: DataFrame with TSC columns
        tsc_columns: List of TSC column names
        
    Returns:
        DataFrame with added normalized TSC columns
    """
    logger.info("Calculating normalized TSC values...")
    
    df = df.copy()
    
    # Calculate sum of all TSC values for each row
    df['_sum_all_tsc'] = df[tsc_columns].sum(axis=1)
    
    # Calculate normalized TSC for each tissue
    for col in tsc_columns:
        # Extract tissue name from column (e.g., "brain TSCs (Not Phosphorylated)" -> "brain")
        tissue_name = col.split()[0].lower()
        normalized_col = f'normalized_{tissue_name}_TSC'
        
        # Calculate: (TSC_tissue / sum_all_TSC) * 1,000,000
        df[normalized_col] = (df[col] / df['_sum_all_tsc']) * 1_000_000
        
        # Replace inf/NaN with 0 (in case sum is 0)
        df[normalized_col] = df[normalized_col].replace([np.inf, -np.inf, np.nan], 0)
    
    # Remove temporary column
    df = df.drop(columns=['_sum_all_tsc'])
    
    logger.info(f"Added {len(tsc_columns)} normalized TSC columns")
    
    return df


def calculate_coefficient_of_variation(df: pd.DataFrame, value_columns: List[str], cv_column_name: str = 'CV') -> pd.DataFrame:
    """
    Calculate coefficient of variation (CV) across tissue values.
    
    Formula: (std_dev / mean) * 100%
    
    Args:
        df: DataFrame with value columns
        value_columns: List of column names to use for CV calculation
        cv_column_name: Name for the CV column
        
    Returns:
        DataFrame with added CV column
    """
    logger.info(f"Calculating coefficient of variation using {len(value_columns)} columns...")
    
    df = df.copy()
    
    # Calculate mean and std dev across all tissue columns for each row
    means = df[value_columns].mean(axis=1)
    std_devs = df[value_columns].std(axis=1)
    
    # Calculate CV: (std_dev / mean) * 100%
    # Handle division by zero
    cv_values = np.where(means != 0, (std_devs / means) * 100, np.nan)
    
    df[cv_column_name] = cv_values
    
    # Round to 2 decimal places
    df[cv_column_name] = df[cv_column_name].round(2)
    
    logger.info(f"Added {cv_column_name} column")
    logger.info(f"CV statistics: mean={df[cv_column_name].mean():.2f}%, std={df[cv_column_name].std():.2f}%")
    
    return df


def process_mouse_data(input_file: str, output_file: str) -> pd.DataFrame:
    """
    Process mouse data to identify housekeeping proteins.
    
    Args:
        input_file: Path to input Excel file
        output_file: Path to output CSV file
        
    Returns:
        Processed DataFrame
    """
    logger.info("=" * 80)
    logger.info("Processing Mouse Data")
    logger.info("=" * 80)
    
    # Define tissue columns for mouse
    mouse_tissue_columns = [
        'brain TSCs (Not Phosphorylated)',
        'brownfat TSCs (Not Phosphorylated)',
        'heart TSCs (Not Phosphorylated)',
        'kidney TSCs (Not Phosphorylated)',
        'liver TSCs (Not Phosphorylated)',
        'lung TSCs (Not Phosphorylated)',
        'pancreas TSCs (Not Phosphorylated)',
        'spleen TSCs (Not Phosphorylated)',
        'testis TSCs (Not Phosphorylated)'
    ]
    
    # Load Excel file
    logger.info(f"Loading Excel file: {input_file}")
    try:
        # For .xls files, let pandas auto-detect or use xlrd
        df = pd.read_excel(input_file)
        logger.info(f"Loaded {len(df)} rows")
    except Exception as e:
        logger.error(f"Error loading Excel file: {e}")
        raise
    
    # Check for IPI column
    if 'IPI' not in df.columns:
        logger.error(f"IPI column not found. Available columns: {list(df.columns)}")
        raise ValueError("IPI column not found in input file")
    
    # Filter housekeeping proteins
    housekeeping_df = filter_housekeeping_proteins(df, mouse_tissue_columns)
    
    if len(housekeeping_df) == 0:
        logger.warning("No housekeeping proteins found!")
        return pd.DataFrame()
    
    # Select output columns: IPI + all TSC columns
    output_columns = ['IPI'] + mouse_tissue_columns
    housekeeping_df = housekeeping_df[output_columns].copy()
    
    # Calculate normalized TSC
    housekeeping_df = calculate_normalized_tsc(housekeeping_df, mouse_tissue_columns)
    
    # Get normalized column names for CV calculation
    normalized_columns = [f'normalized_{col.split()[0].lower()}_TSC' for col in mouse_tissue_columns]
    
    # Calculate CV using normalized TSC values
    housekeeping_df = calculate_coefficient_of_variation(housekeeping_df, normalized_columns, 'CV')
    
    # Save to CSV
    logger.info(f"Saving results to {output_file}...")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    housekeeping_df.to_csv(output_file, index=False)
    logger.info(f"✓ Successfully saved {len(housekeeping_df)} rows")
    
    return housekeeping_df


def process_rat_data(input_file: str, output_file: str) -> pd.DataFrame:
    """
    Process rat data to identify housekeeping proteins.
    
    Args:
        input_file: Path to input Excel file
        output_file: Path to output CSV file
        
    Returns:
        Processed DataFrame
    """
    logger.info("=" * 80)
    logger.info("Processing Rat Data")
    logger.info("=" * 80)
    
    # Define tissue columns for rat
    rat_tissue_columns = [
        'Brain',
        'Heart',
        'Intestine',
        'Kidney',
        'Liver',
        'Lung',
        'Pancreas',
        'PerirenalFat',
        'Spleen',
        'Stomach',
        'Testicle',
        'Thymus'
    ]
    
    # Load Excel file
    logger.info(f"Loading Excel file: {input_file}")
    try:
        # For .xls files, let pandas auto-detect or use xlrd
        df = pd.read_excel(input_file)
        logger.info(f"Loaded {len(df)} rows")
    except Exception as e:
        logger.error(f"Error loading Excel file: {e}")
        raise
    
    # Check for Protein ID column
    if 'Protein ID' not in df.columns:
        logger.error(f"Protein ID column not found. Available columns: {list(df.columns)}")
        raise ValueError("Protein ID column not found in input file")
    
    # Filter housekeeping proteins
    housekeeping_df = filter_housekeeping_proteins(df, rat_tissue_columns)
    
    if len(housekeeping_df) == 0:
        logger.warning("No housekeeping proteins found!")
        return pd.DataFrame()
    
    # Select output columns: Protein ID + all tissue columns
    output_columns = ['Protein ID'] + rat_tissue_columns
    housekeeping_df = housekeeping_df[output_columns].copy()
    
    # Calculate CV using normalized intensity values (already normalized in source data)
    housekeeping_df = calculate_coefficient_of_variation(housekeeping_df, rat_tissue_columns, 'CV')
    
    # Save to CSV
    logger.info(f"Saving results to {output_file}...")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    housekeeping_df.to_csv(output_file, index=False)
    logger.info(f"✓ Successfully saved {len(housekeeping_df)} rows")
    
    return housekeeping_df


def main():
    """Main workflow function."""
    parser = argparse.ArgumentParser(
        description='Identify housekeeping proteins (expressed across all tissue types)'
    )
    parser.add_argument(
        '--species',
        type=str,
        required=True,
        choices=['mouse', 'rat'],
        help='Species: mouse or rat'
    )
    parser.add_argument(
        '--input',
        type=str,
        default=None,
        help='Path to input Excel file (default: species-specific default path)'
    )
    parser.add_argument(
        '--output',
        type=str,
        default=None,
        help='Path to output CSV file (default: data/processed/{species}/housekeeping_proteins.csv)'
    )
    
    args = parser.parse_args()
    
    # Determine script directory and project root
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    
    # Set default input and output paths based on species
    if args.species == 'mouse':
        default_input = os.path.join(project_root, 'data', 'raw', 'mouse', 'protein_expression.xls')
        default_output = os.path.join(project_root, 'data', 'processed', 'mouse', 'housekeeping_proteins.csv')
    else:  # rat
        default_input = os.path.join(project_root, 'data', 'raw', 'rat', '14rat_protein_expression.xls')
        default_output = os.path.join(project_root, 'data', 'processed', 'rat', 'housekeeping_proteins.csv')
    
    input_file = args.input if args.input else default_input
    output_file = args.output if args.output else default_output
    
    # Validate input file
    if not os.path.exists(input_file):
        logger.error(f"Input file not found: {input_file}")
        return 1
    
    logger.info("=" * 80)
    logger.info("Housekeeping Proteins Identification")
    logger.info("=" * 80)
    logger.info(f"Species: {args.species}")
    logger.info(f"Input file: {input_file}")
    logger.info(f"Output file: {output_file}")
    logger.info("")
    
    try:
        # Process based on species
        if args.species == 'mouse':
            result_df = process_mouse_data(input_file, output_file)
        else:  # rat
            result_df = process_rat_data(input_file, output_file)
        
        # Print summary
        logger.info("")
        logger.info("=" * 80)
        logger.info("Summary")
        logger.info("=" * 80)
        logger.info(f"Species: {args.species}")
        logger.info(f"Housekeeping proteins identified: {len(result_df)}")
        if len(result_df) > 0:
            logger.info(f"Output file: {output_file}")
            if 'CV' in result_df.columns:
                logger.info(f"CV statistics:")
                logger.info(f"  Mean CV: {result_df['CV'].mean():.2f}%")
                logger.info(f"  Std CV: {result_df['CV'].std():.2f}%")
                logger.info(f"  Min CV: {result_df['CV'].min():.2f}%")
                logger.info(f"  Max CV: {result_df['CV'].max():.2f}%")
        logger.info("")
        
        return 0
        
    except Exception as e:
        logger.error(f"Error processing data: {e}", exc_info=True)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())


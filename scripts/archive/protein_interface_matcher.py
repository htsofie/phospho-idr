#!/usr/bin/env python3
"""
Protein Interface Matcher

This script searches for protein ID matches between:
- A disorder CSV file (searches protein_id column)
- An interface CSV file (searches P1 and P2 columns)

When matches are found, it creates a new CSV with interface data
and adds 'P1_ID_match' and 'P2_ID_match' columns indicating which protein IDs matched.
Both columns can have values in the same row if both P1 and P2 match.

Usage:
    python protein_interface_matcher.py <disorder_file> <interface_file> <output_file>

Example:
    python protein_interface_matcher.py ../data/processed/mouse/full_disorder.csv ../data/interface/M_musculus_interfacesHQ.csv mouse_interface_matches.csv
"""

import pandas as pd
import sys
import os
from pathlib import Path

def find_protein_interface_matches(disorder_file, interface_file, output_file):
    """
    Find protein ID matches between disorder and interface files.
    
    Args:
        disorder_file (str): Path to disorder CSV file
        interface_file (str): Path to interface CSV file  
        output_file (str): Path for output CSV file
    """
    
    print("=" * 80)
    print("PROTEIN INTERFACE MATCHER")
    print("=" * 80)
    
    # Check if input files exist
    if not os.path.exists(disorder_file):
        print(f"Error: Disorder file not found: {disorder_file}")
        return False
        
    if not os.path.exists(interface_file):
        print(f"Error: Interface file not found: {interface_file}")
        return False
    
    # Load disorder data
    print(f"Loading disorder data from: {disorder_file}")
    try:
        df_disorder = pd.read_csv(disorder_file)
        print(f"Loaded {len(df_disorder)} rows from disorder file")
    except Exception as e:
        print(f"Error loading disorder file: {e}")
        return False
    
    # Check if protein_id column exists
    if 'protein_id' not in df_disorder.columns:
        print(f"Error: 'protein_id' column not found in disorder file")
        print(f"Available columns: {list(df_disorder.columns)}")
        return False
    
    # Load interface data
    print(f"Loading interface data from: {interface_file}")
    try:
        df_interface = pd.read_csv(interface_file)
        print(f"Loaded {len(df_interface)} rows from interface file")
    except Exception as e:
        print(f"Error loading interface file: {e}")
        return False
    
    # Check if P1 and P2 columns exist
    if 'P1' not in df_interface.columns or 'P2' not in df_interface.columns:
        print(f"Error: 'P1' or 'P2' columns not found in interface file")
        print(f"Available columns: {list(df_interface.columns)}")
        return False
    
    # Get unique protein IDs from disorder file
    unique_protein_ids = set(df_disorder['protein_id'].dropna().unique())
    print(f"Found {len(unique_protein_ids)} unique protein IDs in disorder file")
    
    # Find matches
    print("\nSearching for matches...")
    matches = []
    
    for idx, row in df_interface.iterrows():
        p1_id = row['P1']
        p2_id = row['P2']
        p1_match = None
        p2_match = None
        
        # Check P1 column independently
        if p1_id in unique_protein_ids:
            p1_match = p1_id
            
        # Check P2 column independently
        if p2_id in unique_protein_ids:
            p2_match = p2_id
        
        # If we found at least one match, add to results
        if p1_match or p2_match:
            match_row = row.copy()
            match_row['P1_ID_match'] = p1_match if p1_match else None
            match_row['P2_ID_match'] = p2_match if p2_match else None
            matches.append(match_row)
    
    if not matches:
        print("No matches found between disorder and interface files.")
        return False
    
    # Create results DataFrame
    df_matches = pd.DataFrame(matches)
    print(f"Found {len(df_matches)} interface matches")
    
    # Count matches by column
    p1_matches = df_matches['P1_ID_match'].notna().sum()
    p2_matches = df_matches['P2_ID_match'].notna().sum()
    both_matches = df_matches[(df_matches['P1_ID_match'].notna()) & (df_matches['P2_ID_match'].notna())].shape[0]
    
    print(f"Matches in P1 column: {p1_matches}")
    print(f"Matches in P2 column: {p2_matches}")
    print(f"Rows where both P1 and P2 match: {both_matches}")
    
    # Show unique matched proteins
    p1_unique = df_matches['P1_ID_match'].dropna().nunique()
    p2_unique = df_matches['P2_ID_match'].dropna().nunique()
    all_matched_ids = set(df_matches['P1_ID_match'].dropna().unique()) | set(df_matches['P2_ID_match'].dropna().unique())
    print(f"Unique proteins with P1 match: {p1_unique}")
    print(f"Unique proteins with P2 match: {p2_unique}")
    print(f"Total unique proteins with interface data: {len(all_matched_ids)}")
    
    # Save results
    print(f"\nSaving results to: {output_file}")
    try:
        df_matches.to_csv(output_file, index=False)
        print(f"Successfully saved {len(df_matches)} matches to {output_file}")
    except Exception as e:
        print(f"Error saving results: {e}")
        return False
    
    # Show sample of results
    print(f"\nSample of matches:")
    print(df_matches[['P1', 'P2', 'Source', 'P1_ID_match', 'P2_ID_match']].head(10).to_string(index=False))
    
    return True

def main():
    """Main function to handle command line arguments."""
    
    if len(sys.argv) != 4:
        print("Usage: python protein_interface_matcher.py <disorder_file> <interface_file> <output_file>")
        print("\nExample:")
        print("  python protein_interface_matcher.py ../data/processed/mouse/full_disorder.csv ../data/interface/M_musculus_interfacesHQ.csv mouse_interface_matches.csv")
        print("  python protein_interface_matcher.py ../data/processed/rat/full_disorder.csv ../data/interface/H_sapiens_interfacesHQ.csv rat_interface_matches.csv")
        sys.exit(1)
    
    disorder_file = sys.argv[1]
    interface_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    success = find_protein_interface_matches(disorder_file, interface_file, output_file)
    
    if success:
        print(f"\n✅ Analysis complete! Results saved to: {output_file}")
    else:
        print(f"\n❌ Analysis failed. Please check the error messages above.")
        sys.exit(1)

if __name__ == "__main__":
    main()

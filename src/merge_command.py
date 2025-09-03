#!/usr/bin/env python3
"""
Simple command to merge the two cleaned datasets.
Run this after both enhanced_data_cleaner.py and enhanced_uniprot_fetch.py have been completed.
"""

import pandas as pd

def main():
    print("🔗 Merging Enhanced Data with Sequences")
    print("=" * 40)
    
    # Load both datasets
    print("📖 Loading datasets...")
    df_enhanced = pd.read_parquet("../data/processed/cleaned_phosphomouse_site_data.parquet")
    df_sequences = pd.read_parquet("../data/processed/cleaned_phosphorylation_sites_with_sequences.parquet")
    
    print(f"Enhanced data: {len(df_enhanced):,} sites")
    print(f"Sequences data: {len(df_sequences):,} sites")
    
    # Merge sequences into enhanced data
    print("🔗 Merging sequences...")
    df_enhanced['Full Protein Sequence'] = None
    df_enhanced['Sequence Source'] = None
    
    # Create mapping from sequences data
    sequences_map = {}
    core_columns = ['Gene Symbol', 'Protein Description', 'Sequence', 'Residue', 'Site', 'Localized']
    
    for idx, row in df_sequences.iterrows():
        key = tuple(row[col] for col in core_columns)
        sequences_map[key] = {
            'Full Protein Sequence': row['Full Protein Sequence'],
            'Sequence Source': row['Sequence Source']
        }
    
    # Add sequences to enhanced data
    sequences_added = 0
    for idx, row in df_enhanced.iterrows():
        key = tuple(row[col] for col in core_columns)
        if key in sequences_map:
            df_enhanced.at[idx, 'Full Protein Sequence'] = sequences_map[key]['Full Protein Sequence']
            df_enhanced.at[idx, 'Sequence Source'] = sequences_map[key]['Sequence Source']
            sequences_added += 1
    
    # Save merged data
    output_file = "../data/processed/cleaned_phosphomouse_complete.parquet"
    df_enhanced.to_parquet(output_file, index=False)
    
    # Report results
    print(f"\n✅ Merge complete!")
    print(f"Sequences added: {sequences_added:,}/{len(df_enhanced):,}")
    print(f"Final dataset saved to: {output_file}")
    
    # Show final structure
    sequences_with_data = df_enhanced['Full Protein Sequence'].notna().sum()
    print(f"\n📊 Final Results:")
    print(f"Total sites: {len(df_enhanced):,}")
    print(f"Sites with sequences: {sequences_with_data:,} ({sequences_with_data/len(df_enhanced)*100:.1f}%)")
    print(f"Total columns: {len(df_enhanced.columns)}")
    
    print(f"\n📋 All Columns:")
    for i, col in enumerate(df_enhanced.columns):
        print(f"  {i+1:2d}. {col}")

if __name__ == "__main__":
    main()

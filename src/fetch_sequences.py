#!/usr/bin/env python3
"""
Enhanced script to fetch protein sequences using multiple identifier types.
Tries in order: SWISS-PROT → Ensembl → TREMBL → RefSeq
"""

import pandas as pd
import requests
import time
from Bio import ExPASy
from Bio import SwissProt

def get_sequence_via_swiss_prot(swiss_prot_id):
    """Fetch protein sequence from UniProt using SWISS-PROT ID."""
    if pd.isna(swiss_prot_id) or swiss_prot_id is None:
        return None
    
    try:
        clean_id = str(swiss_prot_id).strip()
        handle = ExPASy.get_sprot_raw(clean_id)
        record = SwissProt.read(handle)
        handle.close()
        return record.sequence
    except Exception as e:
        print(f"  ✗ SWISS-PROT {swiss_prot_id}: {e}")
        return None



def get_sequence_via_ensembl(ensembl_id):
    """Fetch protein sequence using Ensembl ID via Ensembl REST API."""
    if pd.isna(ensembl_id) or ensembl_id is None:
        return None
    
    try:
        clean_id = str(ensembl_id).strip()
        url = f"http://rest.ensembl.org/sequence/id/{clean_id}?content-type=text/x-fasta;type=protein"
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            fasta_text = response.text
            # Extract sequence from FASTA format
            lines = fasta_text.strip().split('\n')
            sequence = ''.join(lines[1:])  # Skip header line
            return sequence
        else:
            print(f"  ✗ Ensembl {ensembl_id}: HTTP {response.status_code}")
            return None
    except Exception as e:
        print(f"  ✗ Ensembl {ensembl_id}: {e}")
        return None

def get_sequence_via_refseq(refseq_id):
    """Fetch protein sequence using RefSeq ID via UniProt ID mapping."""
    if pd.isna(refseq_id) or refseq_id is None:
        return None
    
    try:
        clean_id = str(refseq_id).strip()
        # Map RefSeq to UniProt accession
        url = "https://www.uniprot.org/uploadlists/"
        params = {
            "from": "P_REFSEQ_AC",
            "to": "ACC",
            "format": "tab",
            "query": clean_id
        }
        response = requests.post(url, data=params, timeout=10)
        
        if response.status_code == 200 and response.text.strip():
            lines = response.text.strip().split('\n')
            if len(lines) > 1:  # Has data beyond header
                uniprot_id = lines[1].split('\t')[1]
                if uniprot_id and uniprot_id != 'null':
                    # Now get sequence using UniProt ID
                    return get_sequence_via_swiss_prot(uniprot_id)
        
        print(f"  ✗ RefSeq {refseq_id}: No UniProt mapping found")
        return None
    except Exception as e:
        print(f"  ✗ RefSeq {refseq_id}: {e}")
        return None

def get_sequence_via_trembl(trembl_id):
    """Fetch protein sequence using TREMBL ID via UniProt."""
    if pd.isna(trembl_id) or trembl_id is None:
        return None
    
    try:
        clean_id = str(trembl_id).strip()
        # TREMBL IDs are already UniProt IDs, so use directly
        return get_sequence_via_swiss_prot(clean_id)
    except Exception as e:
        print(f"  ✗ TREMBL {trembl_id}: {e}")
        return None

def get_sequence_comprehensive(protein_data):
    """
    Try to get sequence using multiple identifier types in order:
    SWISS-PROT → Ensembl → TREMBL → RefSeq
    """
    swiss_prot_id = protein_data['Protein Identifier (SWISS-PROT)']
    ensembl_id = protein_data['Protein Identifier (Ensembl)']
    refseq_id = protein_data['Protein Identifier (RefSeq)']
    trembl_id = protein_data['Protein Identifier (TREMBL)']
    
    # Try SWISS-PROT first
    if swiss_prot_id:
        sequence = get_sequence_via_swiss_prot(swiss_prot_id)
        if sequence:
            return sequence, 'SWISS-PROT'
    
    # Try Ensembl
    if ensembl_id:
        sequence = get_sequence_via_ensembl(ensembl_id)
        if sequence:
            return sequence, 'Ensembl'
    
    # Try TREMBL
    if trembl_id:
        sequence = get_sequence_via_trembl(trembl_id)
        if sequence:
            return sequence, 'TREMBL'
    
    # Try RefSeq last
    if refseq_id:
        sequence = get_sequence_via_refseq(refseq_id)
        if sequence:
            return sequence, 'RefSeq'
    
    return None, None

def main():
    print("🧬 Enhanced Protein Sequence Fetching")
    print("=" * 50)
    print("Trying identifiers in order: SWISS-PROT → Ensembl → TREMBL → RefSeq")
    print()
    
    # Load the enhanced cleaned data
    print("📖 Loading enhanced cleaned phosphorylation data...")
    df = pd.read_parquet("../data/processed/cleaned_phosphomouse_site_data.parquet")
    print(f"Loaded {len(df):,} phosphorylation sites")
    
    # Get unique proteins (combine all identifier types)
    print("🔍 Identifying unique proteins...")
    unique_proteins = set()
    
    for col in ['Protein Identifier (SWISS-PROT)', 'Protein Identifier (Ensembl)', 
                'Protein Identifier (RefSeq)', 'Protein Identifier (TREMBL)']:
        unique_proteins.update(df[col].dropna().unique())
    
    print(f"Found {len(unique_proteins):,} unique protein identifiers across all types")
    
    # Fetch sequences
    print("🔍 Fetching protein sequences...")
    protein_sequences = {}
    source_stats = {'SWISS-PROT': 0, 'Ensembl': 0, 'RefSeq': 0, 'TREMBL': 0, 'Failed': 0}
    
    # Create a mapping of all identifiers to their protein data
    protein_mapping = {}
    for idx, row in df.iterrows():
        for col in ['Protein Identifier (SWISS-PROT)', 'Protein Identifier (Ensembl)', 
                    'Protein Identifier (RefSeq)', 'Protein Identifier (TREMBL)']:
            identifier = row[col]
            if pd.notna(identifier):
                protein_mapping[identifier] = row
    
    processed = 0
    for identifier in unique_proteins:
        if identifier in protein_mapping:
            protein_data = protein_mapping[identifier]
            processed += 1
            
            print(f"Processing {processed}/{len(unique_proteins)}: {identifier}")
            
            sequence, source = get_sequence_comprehensive(protein_data)
            if sequence:
                protein_sequences[identifier] = sequence
                source_stats[source] += 1
                print(f"  ✓ Success via {source}: {len(sequence)} amino acids")
            else:
                source_stats['Failed'] += 1
                print(f"  ✗ Failed with all methods")
            
            # Be respectful to the servers
            time.sleep(0.5)
    
    # Add sequences to dataframe
    print("\n📝 Adding sequences to dataframe...")
    df['Full Protein Sequence'] = None
    df['Sequence Source'] = None
    
    for idx, row in df.iterrows():
        # Try to find a sequence for this row
        for col in ['Protein Identifier (SWISS-PROT)', 'Protein Identifier (Ensembl)', 
                    'Protein Identifier (RefSeq)', 'Protein Identifier (TREMBL)']:
            identifier = row[col]
            if pd.notna(identifier) and identifier in protein_sequences:
                df.at[idx, 'Full Protein Sequence'] = protein_sequences[identifier]
                # Find which source was used
                sequence, source = get_sequence_comprehensive(row)
                if sequence:
                    df.at[idx, 'Sequence Source'] = source
                break
    
    # Save enhanced data
    output_file = "../data/processed/cleaned_phosphorylation_sites_with_sequences.parquet"
    df.to_parquet(output_file, index=False)
    
    # Report results
    sequences_added = df['Full Protein Sequence'].notna().sum()
    print(f"\n✅ Enhanced sequence fetching complete!")
    print(f"Sequences retrieved: {len(protein_sequences):,}/{len(unique_proteins):,}")
    print(f"Sites with sequences: {sequences_added:,}/{len(df):,} ({sequences_added/len(df)*100:.1f}%)")
    print(f"Enhanced data saved to: {output_file}")
    
    print(f"\n📊 Success by Source:")
    for source, count in source_stats.items():
        if source != 'Failed':
            print(f"  {source:12s}: {count:,} sequences")
    print(f"  {'Failed':12s}: {source_stats['Failed']:,} sequences")

if __name__ == "__main__":
    main()

"""
✅ Enhanced sequence fetching complete!
Sequences retrieved: 6,852/6,940
Sites with sequences: 5,373/5,523 (97.3%)

📊 Success by Source:
  SWISS-PROT  : 5,809 sequences
  Ensembl     : 660 sequences
  RefSeq      : 0 sequences
  TREMBL      : 383 sequences
  Failed      : 88 sequences
"""
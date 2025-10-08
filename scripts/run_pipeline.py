#!/usr/bin/env python3
"""
Pipeline wrapper script to run the complete phosphorylation site mapping pipeline:
clean_data.py -> paper_blast.py -> total_blast.py -> align_to_full_seq.py

Usage:
    python scripts/run_pipeline.py --input data/processed/mouse/test_data.csv --species mouse
"""

import argparse
import subprocess
import sys
import os
import pandas as pd


def run_command(cmd: list, description: str) -> tuple:
    """Run a command and return success status and output path."""
    print(f"\n{'='*80}")
    print(f"STEP: {description}")
    print(f"{'='*80}")
    print(f"Running: {' '.join(cmd)}\n")
    
    result = subprocess.run(cmd, capture_output=False, text=True)
    
    if result.returncode != 0:
        print(f"\n❌ Error: {description} failed with exit code {result.returncode}")
        return False, None
    
    print(f"\n✓ {description} completed successfully")
    return True, None


def remove_dashes_from_sequences(input_file: str, output_file: str, manual_review_only: bool = False) -> bool:
    """Remove dashes from cleaned_site_motif column. Can target all rows or only manual_review==TRUE rows."""
    step_name = "Recleaning manual review rows" if manual_review_only else "Removing dashes from sequences"
    print(f"\n{'='*80}")
    print(f"STEP: {step_name}")
    print(f"{'='*80}")
    
    try:
        # Load the data
        print(f"Loading data from {input_file}")
        df = pd.read_csv(input_file)
        print(f"Loaded {len(df)} rows")
        
        # Check if cleaned_site_motif column exists
        if 'cleaned_site_motif' not in df.columns:
            print("❌ Error: 'cleaned_site_motif' column not found in input file")
            return False
        
        # Determine which rows to process
        if manual_review_only:
            # Check if manual_review column exists
            if 'manual_review' not in df.columns:
                print("❌ Error: 'manual_review' column not found in input file")
                return False
            
            # Filter to manual_review == TRUE rows
            rows_to_process = df[df['manual_review'] == True]
            print(f"Found {len(rows_to_process)} rows with manual_review == TRUE")
            
            if len(rows_to_process) == 0:
                print("No manual review rows found. Nothing to process.")
                return True
        else:
            # Process all rows
            rows_to_process = df
            print(f"Processing all {len(rows_to_process)} rows")
        
        # Count rows with dashes in cleaned_site_motif
        dash_count = 0
        for idx, row in rows_to_process.iterrows():
            if pd.notna(row['cleaned_site_motif']) and '-' in str(row['cleaned_site_motif']):
                dash_count += 1
        
        print(f"Found {dash_count} rows with dashes in cleaned_site_motif")
        
        if dash_count == 0:
            print("No dashes found. Nothing to process.")
            return True
        
        # Remove dashes from the sequences
        cleaned_count = 0
        for idx, row in rows_to_process.iterrows():
            if pd.notna(row['cleaned_site_motif']) and '-' in str(row['cleaned_site_motif']):
                original_motif = str(row['cleaned_site_motif'])
                # Remove all dashes from the sequence
                cleaned_motif = original_motif.replace('-', '')
                df.at[idx, 'cleaned_site_motif'] = cleaned_motif
                cleaned_count += 1
                if manual_review_only:  # Only show details for manual review processing
                    print(f"  Row {idx}: '{original_motif}' → '{cleaned_motif}'")
        
        print(f"Processed {cleaned_count} sequences")
        
        # Save the updated data
        print(f"Saving updated data to {output_file}")
        df.to_csv(output_file, index=False)
        
        print(f"✓ Dash removal completed successfully")
        return True
        
    except Exception as e:
        print(f"❌ Error during dash removal: {e}")
        return False


def reclean_manual_review_rows(input_file: str, output_file: str) -> bool:
    """Reclean cleaned_site_motif column by removing dashes from manual_review==TRUE rows."""
    return remove_dashes_from_sequences(input_file, output_file, manual_review_only=True)


def main():
    parser = argparse.ArgumentParser(
        description='Run the complete phosphorylation site mapping pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process mouse data (from raw/uncleaned data)
  python scripts/run_pipeline.py --input data/processed/mouse/test_data.csv --species mouse
  
  # Process rat data (from raw/uncleaned data)
  python scripts/run_pipeline.py --input data/processed/rat/full_data.csv --species rat
  
  # Skip cleaning step if data is already cleaned
  python scripts/run_pipeline.py --input data/processed/mouse/cleaned_test_data_grouped_processed.csv --species mouse --skip-cleaning
  
  # Use custom BLAST database
  python scripts/run_pipeline.py --input data.csv --species mouse --blast-db data/blast_dbs/custom_mouse.fasta
  
  # Reclean manual review rows from existing dataset (remove dashes and reblast/realign)
  python scripts/run_pipeline.py --input data/processed/mouse/cleaned_full_data_aligned.csv --species mouse --reclean-manual-review
        """
    )
    
    parser.add_argument('--input', '-i', required=True,
                       help='Path to input CSV file (raw data, cleaned data, or existing processed dataset)')
    parser.add_argument('--species', '-s', required=True, choices=['mouse', 'rat'],
                       help='Species name (mouse or rat)')
    parser.add_argument('--skip-cleaning', action='store_true',
                       help='Skip data cleaning step (use if input is already cleaned)')
    parser.add_argument('--blast-db', 
                       help='Path to BLAST database FASTA file (optional, uses default paper DB)')
    parser.add_argument('--output-dir',
                       help='Output directory (default: data/processed/{species}/)')
    parser.add_argument('--keep-intermediate', action='store_true',
                       help='Keep intermediate output files (default: only keep final aligned output)')
    parser.add_argument('--reclean-manual-review', action='store_true',
                       help='Reclean cleaned_site_motif column (remove dashes) and reblast/realign only manual_review==TRUE rows from the specified input dataset')
    
    args = parser.parse_args()
    
    # Validate input file exists
    if not os.path.exists(args.input):
        print(f"❌ Error: Input file not found: {args.input}")
        sys.exit(1)
    
    # Determine output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = os.path.join('data', 'processed', args.species)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Get base filename for naming intermediate files
    input_basename = os.path.splitext(os.path.basename(args.input))[0]
    
    # Define output paths for each stage
    cleaned_output = os.path.join(output_dir, f"cleaned_{input_basename}.csv")
    paper_blast_output = os.path.join(output_dir, f"cleaned_{input_basename}_paper_blast.csv")
    total_blast_output = os.path.join(output_dir, f"cleaned_{input_basename}_total_blast.csv")
    aligned_output = os.path.join(output_dir, f"cleaned_{input_basename}_aligned.csv")
    
    # Handle recleaning workflow - process only manual_review rows from specified dataset
    if args.reclean_manual_review:
        print("\n" + "="*80)
        print("RECLEANING MANUAL REVIEW ROWS WORKFLOW")
        print("="*80)
        print(f"Processing dataset: {args.input}")
        print(f"Species:           {args.species}")
        print(f"Output directory:  {output_dir}")
        print(f"\nThis will:")
        print(f"  1. Reclean cleaned_site_motif column (remove dashes)")
        print(f"  2. Reblast only manual_review==TRUE rows")
        print(f"  3. Realign only manual_review==TRUE rows")
        print(f"  4. Preserve all other rows unchanged")
        
        # Step 1: Reclean the data (only manual_review rows with dashes)
        recleaned_output = os.path.join(output_dir, f"recleaned_{input_basename}.csv")
        success = reclean_manual_review_rows(args.input, recleaned_output)
        if not success:
            sys.exit(1)
        
        # Step 2: Reblast only manual_review rows
        reblast_output = os.path.join(output_dir, f"recleaned_{input_basename}_reblast.csv")
        cmd_reblast = [
            sys.executable, 'scripts/total_blast.py',
            '--input', recleaned_output,
            '--species', args.species,
            '--output', reblast_output
        ]
        if args.blast_db:
            cmd_reblast.extend(['--database', args.blast_db])
        
        success, _ = run_command(cmd_reblast, "Re-BLAST manual review rows")
        if not success:
            sys.exit(1)
        
        # Step 3: Realign only manual_review rows
        realigned_output = os.path.join(output_dir, f"recleaned_{input_basename}_realigned.csv")
        cmd_realign = [
            sys.executable, 'scripts/align_to_full_seq.py',
            '--input', reblast_output,
            '--species', args.species,
            '--output', realigned_output
        ]
        
        success, _ = run_command(cmd_realign, "Re-align manual review rows")
        if not success:
            sys.exit(1)
        
        # Final success message for recleaning workflow
        print(f"\n{'='*80}")
        print("✓ RECLEANING WORKFLOW COMPLETED SUCCESSFULLY")
        print(f"{'='*80}")
        print(f"\nFinal output: {realigned_output}")
        print(f"\n✓ Recleaned sequences, reblasted, and realigned manual review rows.")
        print(f"✓ All other rows preserved unchanged.")
        return
    
    print("\n" + "="*80)
    print("PHOSPHORYLATION SITE MAPPING PIPELINE")
    print("="*80)
    print(f"Input file:    {args.input}")
    print(f"Species:       {args.species}")
    print(f"Output dir:    {output_dir}")
    print(f"\nPipeline stages:")
    
    # Track which stages to run
    stage_num = 1
    intermediate_files = []
    
    # Determine the starting input file
    current_input = args.input
    
    # Stage 0 (optional): clean_data.py
    if not args.skip_cleaning:
        print(f"  {stage_num}. Data Cleaning → {cleaned_output}")
        stage_num += 1
    else:
        print(f"  (Skipping data cleaning)")
        cleaned_output = None
    
    # Stage 0.5: Dash removal (always run after cleaning)
    dash_removed_output = None
    if not args.skip_cleaning:
        dash_removed_output = os.path.join(output_dir, f"dash_removed_{input_basename}.csv")
        print(f"  {stage_num}. Dash Removal → {dash_removed_output}")
        stage_num += 1
    
    print(f"  {stage_num}. Paper BLAST → {paper_blast_output}")
    stage_num += 1
    print(f"  {stage_num}. Total BLAST → {total_blast_output}")
    stage_num += 1
    print(f"  {stage_num}. Alignment   → {aligned_output}")
    
    # Stage 0: clean_data.py (optional)
    if not args.skip_cleaning:
        cmd0 = [
            sys.executable, 'scripts/clean_data.py',
            '--input', current_input,
            '--species', args.species,
            '--output', cleaned_output
        ]
        
        success, _ = run_command(cmd0, "Data Cleaning")
        if not success:
            sys.exit(1)
        
        current_input = cleaned_output
        intermediate_files.append(cleaned_output)
        
        # Stage 0.5: Remove dashes from sequences
        success = remove_dashes_from_sequences(current_input, dash_removed_output, manual_review_only=False)
        if not success:
            sys.exit(1)
        
        current_input = dash_removed_output
        intermediate_files.append(dash_removed_output)
    
    # Stage 1: paper_blast.py
    cmd1 = [
        sys.executable, 'scripts/paper_blast.py',
        '--input', current_input,
        '--species', args.species,
        '--output', paper_blast_output
    ]
    if args.blast_db:
        cmd1.extend(['--database', args.blast_db])
    
    success, _ = run_command(cmd1, "Paper BLAST")
    if not success:
        sys.exit(1)
    
    intermediate_files.append(paper_blast_output)
    
    # Stage 2: total_blast.py
    cmd2 = [
        sys.executable, 'scripts/total_blast.py',
        '--input', paper_blast_output,
        '--species', args.species,
        '--output', total_blast_output
    ]
    
    success, _ = run_command(cmd2, "Total BLAST (manual review rows)")
    if not success:
        sys.exit(1)
    
    intermediate_files.append(total_blast_output)
    
    # Stage 3: align_to_full_seq.py
    cmd3 = [
        sys.executable, 'scripts/align_to_full_seq.py',
        '--input', total_blast_output,
        '--species', args.species,
        '--output', aligned_output
    ]
    
    success, _ = run_command(cmd3, "Sequence Alignment")
    if not success:
        sys.exit(1)
    
    # Clean up intermediate files if not requested to keep
    if not args.keep_intermediate:
        print(f"\n{'='*80}")
        print("Cleaning up intermediate files...")
        print(f"{'='*80}")
        for intermediate_file in intermediate_files:
            try:
                if os.path.exists(intermediate_file):
                    os.remove(intermediate_file)
                    print(f"  Removed: {intermediate_file}")
            except Exception as e:
                print(f"  Warning: Could not remove {intermediate_file}: {e}")
    
    # Final success message
    print(f"\n{'='*80}")
    print("✓ PIPELINE COMPLETED SUCCESSFULLY")
    print(f"{'='*80}")
    print(f"\nFinal output: {aligned_output}")
    
    if args.keep_intermediate:
        print(f"\nIntermediate files:")
        for intermediate_file in intermediate_files:
            print(f"  - {intermediate_file}")
    
    print(f"\nYou can now analyze the results in: {aligned_output}")


if __name__ == "__main__":
    main()


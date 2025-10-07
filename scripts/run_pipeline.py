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
        """
    )
    
    parser.add_argument('--input', '-i', required=True,
                       help='Path to input CSV file (raw or cleaned data)')
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


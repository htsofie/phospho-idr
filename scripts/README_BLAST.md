# BLAST Phosphosite Mapping

This script maps phosphosite sequences for protein groups flagged for manual review using BLAST against SwissProt/TrEMBL databases.

## Features

- **Automatic Database Download**: Downloads and filters SwissProt/TrEMBL databases by species
- **Species Filtering**: Only searches proteins from the target species
- **Group-based Mapping**: Maps all phosphosite sequences in a group to the same protein
- **Quality Scoring**: Selects proteins with best alignment across all sites
- **Position Tracking**: Calculates position differences and alignment quality

## Usage

### Basic Usage

```bash
python scripts/blast_phosphosite_mapper.py \
  --input data/processed/rat/cleaned_test_data_grouped_processed_aligned_remapped.csv \
  --species rat \
  --output results/rat_blast_mapped.csv
```

### Advanced Options

```bash
python scripts/blast_phosphosite_mapper.py \
  --input data/processed/mouse/cleaned_test_data_grouped_processed_aligned_remapped.csv \
  --species mouse \
  --output results/mouse_blast_mapped.csv \
  --db-type trembl \
  --min-identity 85.0
```

## Parameters

- `--input, -i`: Input CSV/Parquet file from remapping script (required)
- `--species, -s`: Species name - mouse, rat, or human (required)
- `--output, -o`: Output CSV file path (required)
- `--db-type`: Database type - swissprot (default) or trembl
- `--min-identity`: Minimum identity threshold for alignments (default: 80.0)

## How It Works

1. **Load Data**: Reads the remapped CSV/Parquet file
2. **Filter Groups**: Identifies groups where `manual_review_flag = True`
3. **Download Database**: Downloads SwissProt/TrEMBL and filters by species
4. **Create BLAST DB**: Builds local BLAST database from filtered sequences
5. **BLAST Sequences**: Searches each phosphosite sequence against the database
6. **Score Proteins**: Ranks candidate proteins by alignment quality and hit count
7. **Align Groups**: Attempts to align all sites in a group to the same protein
8. **Select Best**: Chooses protein with best overall alignment across all sites

## Output

The script creates a CSV file with the following columns:

- `original_protein_id`: Original protein ID from input data
- `mapped_protein_id`: UniProt ID of the best matching protein
- `protein_sequence`: Full protein sequence
- `all_sites_aligned`: Whether all phosphosite sites aligned successfully
- `aligned_sites`: Number of sites that aligned
- `total_sites`: Total number of sites in the group
- `avg_identity`: Average sequence identity across all alignments
- `avg_position_diff`: Average position difference from expected positions

## Database Management

- Databases are downloaded to `data/blast_dbs/` directory
- SwissProt databases are filtered by species before use
- BLAST databases are created automatically and reused
- First run will take longer due to database download

## Requirements

- NCBI BLAST+ installed system-wide
- Biopython (in requirements.txt)
- Internet connection for database download
- Python 3.8+

## Logging

The script provides detailed logging including:
- Database download progress
- BLAST search results
- Alignment quality metrics
- Success/failure status for each group

## Troubleshooting

- **No flagged groups**: Ensure input file has `manual_review_flag` column with `True` values
- **Database download fails**: Check internet connection and UniProt server availability
- **BLAST errors**: Verify NCBI BLAST+ is installed and accessible
- **Memory issues**: Large databases may require significant RAM; consider using smaller species-specific subsets
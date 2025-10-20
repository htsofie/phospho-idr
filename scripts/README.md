# Scripts Documentation

This directory contains Python scripts for processing phosphorylation data across multiple species (mouse, rat). The scripts implement a comprehensive pipeline for data cleaning, BLAST analysis, and sequence alignment.

## Complete Workflow

The data processing pipeline follows this sequence:

```
gen_data.py → clean_data.py → paper_blast.py → total_blast.py → align_to_full_seq.py
```

### Quick Start (Automated Pipeline)

For full dataset processing, use the provided shell scripts:

```bash
# Process mouse data
bash scripts/mouse_blast.sh
# Or to get detailed output as script runs:
bash -x scripts/mouse_blast.sh
# Process rat data  
bash scripts/rat_blast.sh
```

## Script Descriptions

### 1. Data Generation (`gen_data.py`)

Generates datasets from raw phosphorylation data with configurable sampling and output formats.

**Purpose**: Convert raw Excel files to processed CSV/Parquet format with optional sampling.

**Usage**:
```bash
# Generate sample dataset (30 proteins)
python scripts/gen_data.py --config configs/rat.yaml --mode sample --sample-size 30 --output-format csv

# Generate full dataset
python scripts/gen_data.py --config configs/mouse.yaml --mode full --output-format parquet
```

**Parameters**:
- `--config, -c`: Path to YAML configuration file (required)
- `--mode, -m`: Processing mode - 'sample' or 'full' (default: sample)
- `--sample-size, -n`: Number of proteins to sample (default: 30)
- `--random-state, -r`: Random seed for reproducibility (default: 42)
- `--output-format, -f`: Output format - 'csv' or 'parquet' (default: csv)

### 2. Data Cleaning (`clean_data.py`)

Cleans and filters phosphorylation data based on species-specific criteria.

**Purpose**: Remove ambiguous sites, filter by localization probability, and standardize data format.

**Usage**:
```bash
# Clean mouse data
python scripts/clean_data.py --input data/processed/mouse/full_data.csv --species mouse --output data/processed/mouse/cleaned_full_data.csv

# Clean rat data
python scripts/clean_data.py --input data/processed/rat/test_data.csv --species rat --output data/processed/rat/cleaned_test_data.csv
```

**Parameters**:
- `--input, -i`: Input CSV/Parquet file (required)
- `--species, -s`: Species name - 'mouse' or 'rat' (required)
- `--output, -o`: Output file path (required)

**Filtering Criteria**:
- **Mouse**: Removes sites with ambiguous localization ('Amb')
- **Rat**: Filters sites with localization probability < 0.75

### 3. Paper BLAST Analysis (`paper_blast.py`)

Performs BLAST searches against species-specific UniProtKB databases from literature papers.

**Purpose**: Map phosphosites to protein sequences using curated paper databases.

**Usage**:
```bash
# BLAST mouse data against paper database
python scripts/paper_blast.py --input data/processed/mouse/cleaned_full_data.csv --species mouse --output data/processed/mouse/full_paper_blast.csv

# BLAST rat data against paper database
python scripts/paper_blast.py --input data/processed/rat/cleaned_full_data.csv --species rat --output data/processed/rat/full_paper_blast.csv
```

**Parameters**:
- `--input, -i`: Input cleaned data file (required)
- `--species, -s`: Species name - 'mouse' or 'rat' (required)
- `--output, -o`: Output file path (required)

**Features**:
- Uses species-specific BLAST databases from papers
- Handles protein groups with multiple identifiers
- Extracts top BLAST matches with alignment details

### 4. Total BLAST Analysis (`total_blast.py`)

Performs comprehensive BLAST searches against complete UniProtKB species databases.

**Purpose**: Map phosphosites using full UniProtKB databases for comprehensive coverage.

**Usage**:
```bash
# Total BLAST for mouse data
python scripts/total_blast.py --input data/processed/mouse/full_paper_blast.csv --species mouse --output data/processed/mouse/full_total_blast.csv

# Total BLAST for rat data
python scripts/total_blast.py --input data/processed/rat/full_paper_blast.csv --species rat --output data/processed/rat/full_total_blast.csv
```

**Parameters**:
- `--input, -i`: Input paper BLAST results file (required)
- `--species, -s`: Species name - 'mouse' or 'rat' (required)
- `--output, -o`: Output file path (required)

**Features**:
- Filters for manually reviewed sites (`manual_review == TRUE`)
- Uses complete UniProtKB databases
- Provides comprehensive sequence mapping

### 5. Sequence Alignment (`align_to_full_seq.py`)

Aligns phosphosite sequences to full protein sequences and extracts context.

**Purpose**: Map phosphosites to their positions in full protein sequences and extract surrounding context.

**Usage**:
```bash
# Align mouse sequences
python scripts/align_to_full_seq.py --input data/processed/mouse/full_total_blast.csv --species mouse

# Align rat sequences
python scripts/align_to_full_seq.py --input data/processed/rat/full_total_blast.csv --species rat
```

**Parameters**:
- `--input, -i`: Input BLAST results file (required)
- `--species, -s`: Species name - 'mouse' or 'rat' (required)

**Output**:
- Creates aligned sequences with phosphosite context
- Generates `*_aligned.csv` files with full sequence information

## Shell Scripts

### Automated Pipeline Scripts

**`mouse_blast.sh`**: Complete mouse data processing pipeline
```bash
bash scripts/mouse_blast.sh
```

**`rat_blast.sh`**: Complete rat data processing pipeline
```bash
bash scripts/rat_blast.sh
```

These scripts automatically run the complete workflow:
1. Clean the data
2. Run paper BLAST analysis
3. Run total BLAST analysis
4. Align sequences to full proteins

## Configuration Files

### Species-Specific Configurations

**`configs/mouse.yaml`**: Mouse data configuration
- Raw data: `Phosphomouse_phosphorylation_sites.xlsb`
- Tissue-specific columns (Brain, Heart, Liver, etc.)
- Localization filter: Excludes ambiguous sites

**`configs/rat.yaml`**: Rat data configuration
- Raw data: `14rat_data.xls`
- Tissue-specific columns (Cortex, Brainstem, etc.)
- Localization filter: Probability threshold ≥ 0.75

## Data Requirements

### Input Data
- Raw Excel files (.xls, .xlsb) in `data/raw/{species}/`
- BLAST databases in `data/blast_dbs/`
- Processed data in `data/processed/{species}/`

### Dependencies
- Python 3.8+
- BioPython
- pandas
- pyxlsb (for .xlsb files)
- BLAST+ command line tools

## Output Files

The pipeline generates several intermediate and final output files:

1. **`*_data.csv`**: Initial processed data
2. **`cleaned_*_data.csv`**: Cleaned and filtered data
3. **`*_paper_blast.csv`**: Paper BLAST results
4. **`*_total_blast.csv`**: Total BLAST results
5. **`*_aligned.csv`**: Final aligned sequences

## Troubleshooting

### Common Issues

1. **BLAST not found**: Ensure BLAST+ is installed and in PATH
2. **Database missing**: Check that BLAST databases exist in `data/blast_dbs/`
3. **Memory issues**: Use sample mode for large datasets
4. **File format errors**: Ensure input files are in correct format (CSV/Parquet)

### Logging

All scripts provide detailed logging output. Check console output for processing status and any error messages.

## Future Improvements

- Nextflow pipeline integration for better workflow management
- Parallel processing for large datasets
- Additional species support (human)
- Enhanced error handling and recovery
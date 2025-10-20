# Phosphorylation Data Analysis Project

A comprehensive data analysis pipeline for processing and analyzing phosphorylation datasets across multiple species (mouse, rat). This repository contains scripts, notebooks, and data for cleaning, processing, BLAST analysis, and sequence alignment of phosphorylation site data.

## Overview

This project processes phosphorylation site data from multiple species, performing comprehensive sequence analysis including:
- Data cleaning and filtering
- BLAST analysis against species-specific databases
- Sequence alignment and context extraction
- Multi-species comparative analysis

## Quick Start

### Prerequisites
- Python 3.8+
- Virtual environment (recommended)

### Installation
```bash
# Clone the repository
git clone https://github.com/htsofie/phospho-idr
cd phospho_root

# Create and activate virtual environment
python3 -m venv lab_env
source lab_env/bin/activate  # On Windows: lab_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage

#### Quick Start (Automated Pipeline)
```bash
# Activate virtual environment
source lab_env/bin/activate

# Process mouse data (complete pipeline)
bash scripts/mouse_blast.sh

# Process rat data (complete pipeline)
bash scripts/rat_blast.sh
```

#### Manual Step-by-Step Processing
```bash
# 1. Generate data from raw files
python scripts/gen_data.py --config configs/mouse.yaml --mode full --output-format csv

# 2. Clean the data
python scripts/clean_data.py --input data/processed/mouse/full_data.csv --species mouse --output data/processed/mouse/cleaned_full_data.csv

# 3. Run paper BLAST analysis
python scripts/paper_blast.py --input data/processed/mouse/cleaned_full_data.csv --species mouse --output data/processed/mouse/full_paper_blast.csv

# 4. Run total BLAST analysis
python scripts/total_blast.py --input data/processed/mouse/full_paper_blast.csv --species mouse --output data/processed/mouse/full_total_blast.csv

# 5. Align sequences to full proteins
python scripts/align_to_full_seq.py --input data/processed/mouse/full_total_blast.csv --species mouse
```
### To Push/Pull Changes to/from GitHub
```bash
# To push changes
git status #checks status of what is new
git add . #adds all modified files
git commit -m "Add message here on new changes"
git push origin master # use master not main unless i want to rename branch

# To pull updates
cd /home/htsofie/Desktop/phospho_root && git fetch --all --prune | cat && git pull --rebase --autostash | cat
# OR better..
cd git pull origin master
```
## Project Structure

```
phospho_root/
├── data/                             # Data storage organized by species
│   ├── raw/                          # Original data files
│   │   ├── rat/                      # Rat phosphorylation data
│   │   │   └── 14rat_data.xls        # Rat phosphorylation sites dataset
│   │   └── mouse/                    # Mouse phosphorylation data
│   │       └── Phosphomouse_phosphorylation_sites.xlsb
│   ├── processed/                    # Processed and cleaned data
│   │   ├── rat/                      # Rat processed data
│   │   │   ├── full_data.csv         # Full rat dataset
│   │   │   ├── cleaned_full_data.csv # Cleaned rat data
│   │   │   ├── full_paper_blast.csv  # Paper BLAST results
│   │   │   ├── full_total_blast.csv  # Total BLAST results
│   │   │   └── full_total_blast_aligned.csv # Final aligned sequences
│   │   └── mouse/                    # Mouse processed data
│   │       ├── full_data.csv         # Full mouse dataset
│   │       ├── cleaned_full_data.csv # Cleaned mouse data
│   │       ├── full_paper_blast.csv  # Paper BLAST results
│   │       ├── full_total_blast.csv  # Total BLAST results
│   │       └── full_total_blast_aligned.csv # Final aligned sequences
│   └── blast_dbs/                    # BLAST databases
│       ├── mouse_full_blast.*        # Mouse full UniProtKB database
│       ├── mouse_paper_blast.*       # Mouse paper database
│       ├── rat_full_blast.*          # Rat full UniProtKB database
│       └── rat_paper_blast.*         # Rat paper database
│
├── scripts/                          # Python scripts for data processing
│   ├── gen_data.py                   # Generate datasets from raw data
│   ├── clean_data.py                 # Clean and filter data
│   ├── paper_blast.py                # BLAST against paper databases
│   ├── total_blast.py                # BLAST against full UniProtKB
│   ├── align_to_full_seq.py          # Align sequences to full proteins
│   ├── mouse_blast.sh                # Mouse processing pipeline
│   ├── rat_blast.sh                  # Rat processing pipeline
│   ├── README.md                     # Scripts documentation
│   └── old/                          # Archived scripts
│
├── notebooks/                        # Jupyter notebooks and analysis
│   ├── alignment_check.ipynb         # Alignment validation
│   ├── percent_phosphosites.ipynb    # Phosphosite analysis
│   └── README.md                     # Notebooks documentation
│
├── configs/                          # Configuration files
│   ├── mouse.yaml                    # Mouse data configuration
│   └── rat.yaml                      # Rat data configuration
│
├── docs/                             # Documentation
│   ├── general_decisions.md          # Project decisions
│   └── rat_data.md                   # Rat data documentation
│
├── lab_env/                         # Python virtual environment
├── requirements.txt                 # Python package dependencies
├── lab_env_requirements.txt         # Lab environment requirements
└── README.md                        # This file
```

## Key Features

### Complete Data Processing Pipeline
1. **Data Generation**: Convert raw Excel files (.xls, .xlsb) to processed CSV/Parquet format
2. **Data Cleaning**: Filter ambiguous sites and apply species-specific quality filters
3. **BLAST Analysis**: Map phosphosites to protein sequences using two-tier approach:
   - Paper databases (curated from literature)
   - Full UniProtKB databases (comprehensive coverage)
4. **Sequence Alignment**: Align phosphosite sequences to full protein sequences
5. **Context Extraction**: Extract surrounding sequence context for analysis

### Multi-Species Support
- **Mouse**: Phosphomouse phosphorylation sites dataset (~35,965 sites)
- **Rat**: 14rat phosphorylation sites dataset (~23,015 sites)
- **Extensible**: Ready for additional species (human, etc.)

### Analysis Tools
- **Automated Pipelines**: Shell scripts for complete data processing
- **Interactive Analysis**: Jupyter notebooks for data exploration
- **BLAST Integration**: Comprehensive sequence mapping and alignment
- **Quality Control**: Built-in filtering and validation steps

## Data Sources

### Mouse Data
- **Source**: Phosphomouse phosphorylation sites
- **Format**: Excel (.xlsb) → Parquet
- **Size**: ~35,965 phosphorylation sites
- **Features**: Protein identifiers, localization probabilities, tissue-specific data

### Rat Data
- **Source**: 14rat phosphorylation sites
- **Format**: Excel (.xls) → CSV/Parquet
- **Size**: ~23,015 phosphorylation sites
- **Features**: UniProt IDs, tissue distribution, kinase motifs

## Data Processing Workflow

### Step-by-Step Pipeline

1. **Data Generation** (`gen_data.py`)
   - Converts raw Excel files to processed format
   - Supports both sample and full dataset generation
   - Configurable output formats (CSV/Parquet)

2. **Data Cleaning** (`clean_data.py`)
   - Applies species-specific quality filters
   - Removes ambiguous phosphorylation sites
   - Standardizes data format

3. **Paper BLAST** (`paper_blast.py`)
   - Maps phosphosites using curated paper databases
   - Handles protein groups with multiple identifiers
   - Extracts top BLAST matches

4. **Total BLAST** (`total_blast.py`)
   - Comprehensive mapping using full UniProtKB databases
   - Filters for manually reviewed sites
   - Provides complete sequence coverage

5. **Sequence Alignment** (`align_to_full_seq.py`)
   - Aligns phosphosite sequences to full proteins
   - Extracts surrounding sequence context
   - Generates final aligned datasets

### Scripts Overview

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `gen_data.py` | Generate datasets from raw data | Raw Excel files | Processed CSV/Parquet |
| `clean_data.py` | Clean and filter data | Processed data | Cleaned data |
| `paper_blast.py` | BLAST against paper databases | Cleaned data | Paper BLAST results |
| `total_blast.py` | BLAST against full UniProtKB | Paper BLAST results | Total BLAST results |
| `align_to_full_seq.py` | Align sequences to full proteins | BLAST results | Aligned sequences |

## Usage Examples

### Complete Pipeline (Automated)
```bash
# Process mouse data
bash scripts/mouse_blast.sh

# Process rat data
bash scripts/rat_blast.sh
```

### Individual Script Usage
```bash
# Generate sample dataset
python scripts/gen_data.py --config configs/rat.yaml --mode sample --sample-size 100

# Clean data with custom parameters
python scripts/clean_data.py --input data/processed/mouse/full_data.csv --species mouse --output data/processed/mouse/cleaned_data.csv

# Run BLAST analysis
python scripts/paper_blast.py --input data/processed/mouse/cleaned_data.csv --species mouse --output data/processed/mouse/paper_blast.csv
```

## Dependencies

### Core Requirements
- **Python 3.8+**: Required Python version
- **pandas**: Data manipulation and analysis
- **numpy**: Numerical computing
- **pyxlsb**: Excel .xlsb file reading
- **BioPython**: Biological sequence analysis
- **PyYAML**: Configuration file parsing

### BLAST Requirements
- **BLAST+**: Command-line BLAST tools
- **NCBI BLAST databases**: Species-specific protein databases

### Jupyter Ecosystem
- **jupyter**: Interactive computing
- **jupyterlab**: Modern Jupyter interface
- **ipython**: Enhanced Python shell

### Installation
```bash
# Install Python dependencies
pip install -r requirements.txt

# Install BLAST+ (system-dependent)
# On macOS: brew install blast
# On Ubuntu: sudo apt-get install ncbi-blast+
```

See `requirements.txt` and `lab_env_requirements.txt` for complete dependency lists.

## Configuration

### Species-Specific Settings

Each species has its own configuration file in `configs/`:

- **`mouse.yaml`**: Mouse-specific column mappings and filters
- **`rat.yaml`**: Rat-specific column mappings and filters

### BLAST Databases

Required BLAST databases in `data/blast_dbs/`:
- Species-specific paper databases (curated from literature)
- Full UniProtKB species databases (comprehensive coverage)

## Output Files

The pipeline generates several intermediate and final output files:

1. **`*_data.csv`**: Initial processed data from raw files
2. **`cleaned_*_data.csv`**: Cleaned and filtered data
3. **`*_paper_blast.csv`**: BLAST results against paper databases
4. **`*_total_blast.csv`**: BLAST results against full UniProtKB
5. **`*_aligned.csv`**: Final aligned sequences with context

## Troubleshooting

### Common Issues

1. **BLAST not found**: Ensure BLAST+ is installed and in PATH
2. **Database missing**: Check that BLAST databases exist in `data/blast_dbs/`
3. **Memory issues**: Use sample mode for large datasets
4. **File format errors**: Ensure input files are in correct format (CSV/Parquet)

### Logging

All scripts provide detailed logging output. Check console output for processing status and any error messages.

## Contributing

This project is designed to work across Mac and Linux systems. When contributing:
- Ensure compatibility with both operating systems
- Test with both mouse and rat datasets
- Update documentation for any new features
- Follow the existing code structure and naming conventions

## Future Improvements

- Nextflow pipeline integration for better workflow management
- Parallel processing for large datasets
- Additional species support (human)
- Enhanced error handling and recovery
- Web interface for data exploration
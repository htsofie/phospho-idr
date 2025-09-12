# Phosphorylation Data Analysis Project

A comprehensive data analysis project for processing and analyzing phosphorylation datasets across multiple species (mouse, rat, human). This repository contains scripts, notebooks, and data for cleaning, processing, and analyzing phosphorylation site data.

## Quick Start

### Prerequisites
- Python 3.8+
- Virtual environment (recommended)

### Installation
```bash
# Clone the repository
git clone <repository-url>
cd phospho_root

# Create and activate virtual environment
python3 -m venv lab_env
source lab_env/bin/activate  # On Windows: lab_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

### Basic Usage
```bash
# Activate virtual environment
source lab_env/bin/activate

```

##  Project Structure

```
phospho_root/
├── data/                             # Data storage organized by species
│   ├── raw/                          # Original data files
│   │   ├── human/                    # Human phosphorylation data (empty)
│   │   ├── rat/                      # Rat phosphorylation data
│   │   │   └── 14rat_data.xls        # Rat phosphorylation sites dataset
│   │   └── mouse/                    # Mouse phosphorylation data
│   │       └── Phosphomouse_phosphorylation_sites.xlsb
│   └── processed/                    # Processed and cleaned data
│       ├── human/                    # Human processed data (empty)
│       ├── rat/                      # Rat processed data
│       │   ├── test_data.csv         # Sample rat data
│       │   └── test_data_cleaned.csv # Cleaned rat data
│       └── mouse/                    # Mouse processed data
│           ├── Phosphomouse_phosphorylation_sites.parquet
│           └── cleaned_phosphomouse_site_data.parquet
│
├── scripts/                          # Python scripts for data processing
│   ├── convert_to_parquet.py         # Convert Excel files to Parquet format
│   ├── enhanced_clean_data.py        # Clean and process phosphorylation data
│   ├── fetch_sequences.py            # Fetch protein sequences from UniProt
│   ├── gen_test_data.py              # Generate test datasets
│   ├── merge_command.py              # Merge datasets with sequences
│   └── old/                          # Archived scripts
│       ├── enhanced_clean.py         # Previous version of data cleaning
│       ├── get_full_seq.py           # Previous sequence fetching
│       ├── merge_data.py             # Previous merge script
│       └── README.md                 # Documentation for archived scripts
│
├── notebooks/                        # Jupyter notebooks and analysis scripts
│   ├── read_data.py                  # Data reading and exploration
│   ├── search_data.py                # UniProt ID search functionality
│   └── README.md                     # Notebooks documentation
│
├── outputs/                          # Analysis results and reports
│   └── cleaned_data_results.md       # Data cleaning summary report
│
├── ab_env/                          # Python virtual environment
│   └── ...                          # Virtual environment files
│
├── requirements.txt                 # Python package dependencies
└── README.md                        # This file
```

## Key Features

### Data Processing Pipeline
1. **Data Conversion**: Convert Excel files (.xls, .xlsb) to efficient Parquet format
2. **Data Cleaning**: Extract protein identifiers, filter localized sites
3. **Sequence Fetching**: Retrieve protein sequences from UniProt database
4. **Data Merging**: Combine phosphorylation data with protein sequences

### Multi-Species Support
- **Mouse**: Phosphomouse phosphorylation sites dataset
- **Rat**: 14rat phosphorylation sites dataset  
- **Human**: Ready for human phosphorylation data

### Analysis Tools
- **UniProt ID Search**: Find specific proteins in datasets
- **Data Exploration**: Jupyter notebooks for interactive analysis
- **Test Data Generation**: Create sample datasets for testing

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

##  Scripts Overview

| Script | Purpose | Input | Output |
|--------|---------|-------|--------|
| `convert_to_parquet.py` | Convert Excel to Parquet | `.xls/.xlsb` | `.parquet` |
| `enhanced_clean_data.py` | Clean and process data | Raw parquet | Cleaned parquet |
| `fetch_sequences.py` | Fetch protein sequences | Cleaned data | Data with sequences |
| `merge_command.py` | Merge datasets | Multiple parquet | Combined dataset |
| `gen_test_data.py` | Generate test data | Full dataset | Sample dataset |
| `search_data.py` | Search UniProt IDs | Dataset + ID | Search results |

## Usage Examples

### Search for a specific UniProt ID
```python
# In notebooks/search_data.py
target_id = "Q9R1N3"  # Example UniProt ID
result = search_uniprot_id("data/raw/rat/14rat_data.xls", target_id)
```

### Generate test dataset
```python
# In scripts/gen_test_data.py
df = pd.read_excel("../data/raw/rat/14rat_data.xls")
sample = df.sample(n=20, random_state=42)
sample.to_csv("../data/processed/rat/test_data.csv", index=False)
```

## Dependencies

### Core Requirements
- **pandas**: Data manipulation and analysis
- **numpy**: Numerical computing
- **xlrd**: Excel file reading
- **requests**: HTTP requests for API calls

### Jupyter Ecosystem
- **jupyter**: Interactive computing
- **jupyterlab**: Modern Jupyter interface
- **ipython**: Enhanced Python shell

See `requirements.txt` for complete dependency list.


**Note**: This project is designed to work across Mac and Linux systems. Ensure you're using the correct virtual environment for your operating system.
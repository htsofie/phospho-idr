# Phosphorylation Data Analysis Project

A comprehensive data analysis project for processing and analyzing phosphorylation datasets across multiple species (mouse, rat, human). This repository contains scripts, notebooks, and data for cleaning, processing, and analyzing phosphorylation site data.

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
```bash
# Activate virtual environment
source lab_env/bin/activate

```
### To Push/Pull Changes to/from GitHub
```bash
git status #checks status of what is new
git add . #adds all modified files
git commit -m "Add message here on new changes"
git push origin master # use master not main unless i want to rename branch
```
##  Project Structure

```
phospho_root/
├── data/                             # Data storage organized by species
│   ├── raw/                          # Original data files
│   │   ├── rat/                      # Rat phosphorylation data
│   │   │   ├── 14rat_data.xls        # Rat phosphorylation sites dataset
│   │   │   └── .~lock.14rat_data.xls# # Excel lock file
│   │   └── mouse/                    # Mouse phosphorylation data
│   │       └── Phosphomouse_phosphorylation_sites.xlsb
│   └── processed/                    # Processed and cleaned data
│       ├── rat/                      # Rat processed data
│       │   ├── test_data.csv         # Sample rat data
│       │   ├── test_data_cleaned.csv # Cleaned rat data
│       │   └── .~lock.test_data.csv# # Excel lock file
│       └── mouse/                    # Mouse processed data
│           ├── Phosphomouse_phosphorylation_sites.parquet
│           └── cleaned_phosphomouse_site_data.parquet
│
├── scripts/                          # Python scripts for data processing
│   ├── convert_to_parquet.py         # Convert Excel files to Parquet format
│   ├── gen_test_data.py              # Generate test datasets
│   ├── README.md                     # Scripts documentation
│   └── old/                          # Archived scripts (empty)
│
├── notebooks/                        # Jupyter notebooks and analysis scripts
│   ├── search_data.py                # UniProt ID search functionality
│   └── README.md                     # Notebooks documentation
│
├── configs/                          # Configuration files
│   └── rat.yaml                      # Rat data processing configuration
│
├── lab_env/                         # Python virtual environment
│   ├── bin/                         # Virtual environment executables
│   ├── lib/                         # Python packages and libraries
│   ├── etc/                         # Configuration files
│   └── share/                       # Shared data and applications
│
├── .git/                            # Git repository metadata
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
df = pd.read_excel("data/raw/rat/14rat_data.xls")
sample = df.sample(n=20, random_state=42)
sample.to_csv("data/processed/rat/test_data.csv", index=False)
```

### Convert Excel to Parquet
```python
# In scripts/convert_to_parquet.py
convert_excel_to_parquet("data/raw/mouse/Phosphomouse_phosphorylation_sites.xlsb", 
                        "data/processed/mouse/Phosphomouse_phosphorylation_sites.parquet")
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
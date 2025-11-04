# Phospho Analysis Project - Environment Setup

This project supports both Linux and macOS environments. Follow the appropriate setup instructions below.

## Quick Setup

### For Linux (Ubuntu/Debian):
```bash
./setup_environment.sh
```

### For macOS:
```bash
./setup_macos.sh
```

## Manual Setup

If you prefer to set up manually or the scripts don't work:

### Linux Setup:
1. Ensure Python 3.8+ is installed
2. Install packages: `pip3 install --break-system-packages -r requirements.txt`
3. Register kernel: `python3 -m ipykernel install --user --name=phospho_env --display-name="Phospho Analysis Environment"`

### macOS Setup:
1. Install Homebrew if not already installed
2. Install Python 3: `brew install python3`
3. Install packages: `pip3 install --user -r requirements.txt`
4. Register kernel: `python3 -m ipykernel install --user --name=phospho_env --display-name="Phospho Analysis Environment"`

## Using Jupyter Notebooks

After setup, you can start Jupyter in several ways:

1. **Jupyter Lab** (recommended): `jupyter lab`
2. **Jupyter Notebook**: `jupyter notebook`
3. **VS Code**: Open `.ipynb` files and select the "Phospho Analysis Environment" kernel

## Environment Details

- **Python Version**: 3.8+ (tested on 3.12)
- **Key Packages**: pandas, numpy, matplotlib, seaborn, biopython, jupyter
- **Kernel Name**: "Phospho Analysis Environment"
- **Cross-Platform**: Linux (Ubuntu/Debian) and macOS

## Troubleshooting

### Linux Issues:
- If you get "externally-managed-environment" error, use `--break-system-packages` flag
- If pip is not available, install: `sudo apt install python3-pip`

### macOS Issues:
- If Homebrew is not installed, install it first
- If Python 3 is not available, install via Homebrew: `brew install python3`
- If you get permission errors, use `--user` flag with pip

### General Issues:
- Ensure Python 3.8+ is installed
- Check that all packages imported successfully
- Verify kernel is registered: `jupyter kernelspec list`

## Project Structure

```
phospho_root/
├── notebooks/           # Jupyter notebooks for analysis
├── scripts/            # Python scripts for data processing
├── data/               # Data files (raw and processed)
├── configs/            # Configuration files
├── requirements.txt    # Python package requirements
├── setup_environment.sh  # Linux setup script
├── setup_macos.sh     # macOS setup script
└── README.md          # This file
```

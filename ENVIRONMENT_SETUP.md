# Phospho Analysis Project - Environment Setup

This project supports both Linux and macOS environments with a single unified setup script.

## Quick Setup

### For All Platforms (Linux & macOS):
```bash
bash setup.sh
```

That's it! The script will:
1. Check for Python 3.8+
2. Create a virtual environment (`lab_env`)
3. Install all required packages
4. Register the Jupyter kernel

## Manual Setup (Alternative)

If you prefer to set up manually:

```bash
# Create virtual environment
python3 -m venv lab_env

# Activate virtual environment
source lab_env/bin/activate  # On Windows: lab_env\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Register Jupyter kernel
python -m ipykernel install --user --name=lab_env --display-name="Python (lab_env)"
```

## Using the Environment

After setup:

```bash
# Activate the environment
source lab_env/bin/activate

# Start Jupyter Lab
jupyter lab

# Or start Jupyter Notebook
jupyter notebook
```

In your notebooks, select the **"Python (lab_env)"** kernel.

## Environment Details

- **Python Version**: 3.8+ (tested on 3.13)
- **Virtual Environment**: `lab_env/` (created in project root)
- **Key Packages**: pandas, numpy, matplotlib, seaborn, biopython, jupyter
- **Kernel Name**: "Python (lab_env)"
- **Cross-Platform**: Linux and macOS

## Troubleshooting

### Python Not Found:
- **Linux**: `sudo apt install python3`
- **macOS**: `brew install python3`

### Virtual Environment Issues:
- If `lab_env` is corrupted, delete it and rerun `bash setup.sh`
- Make sure you're using `source lab_env/bin/activate` (not just `lab_env/bin/activate`)

### Package Installation Issues:
- Ensure you have internet connection
- Try upgrading pip first: `lab_env/bin/pip install --upgrade pip`

### Jupyter Kernel Not Found:
- Verify kernel is registered: `jupyter kernelspec list`
- Re-register: `lab_env/bin/python -m ipykernel install --user --name=lab_env`

## Project Structure

```
phospho_root/
├── notebooks/           # Jupyter notebooks for analysis
├── scripts/            # Python scripts for data processing
├── data/               # Data files (raw and processed)
├── configs/            # Configuration files
├── requirements.txt    # Python package requirements
├── setup.sh           # Unified setup script (Mac & Linux)
└── README.md          # Main project documentation
```

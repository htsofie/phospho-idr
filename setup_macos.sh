#!/bin/bash

# macOS-specific setup script for Phospho Analysis Project
# Run this script on macOS to set up the environment

set -e  # Exit on any error

echo "Setting up Phospho Analysis Environment on macOS..."

# Check if Homebrew is installed
if ! command -v brew &> /dev/null; then
    echo "Homebrew not found. Installing Homebrew..."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
fi

# Install Python 3 if not already installed
if ! command -v python3 &> /dev/null; then
    echo "Installing Python 3 via Homebrew..."
    brew install python3
fi

# Install pip if not available
if ! command -v pip3 &> /dev/null; then
    echo "Installing pip..."
    python3 -m ensurepip --upgrade
fi

# Install packages
echo "Installing Python packages..."
pip3 install --user -r requirements.txt

# Register the kernel
echo "Registering Jupyter kernel..."
python3 -m ipykernel install --user --name=phospho_env --display-name="Phospho Analysis Environment"

# Test installation
echo "Testing installation..."
python3 -c "
import sys
try:
    import jupyter, ipykernel, pandas, matplotlib, seaborn, numpy
    import Bio
    from Bio import __version__ as bio_version
    print('✓ All packages imported successfully!')
    print(f'✓ Python version: {sys.version}')
    print(f'✓ Biopython version: {bio_version}')
    print('✓ Environment setup complete!')
except ImportError as e:
    print(f'✗ Import error: {e}')
    sys.exit(1)
"

echo ""
echo "Setup complete! You can now:"
echo "1. Start Jupyter Lab: jupyter lab"
echo "2. Start Jupyter Notebook: jupyter notebook"
echo "3. Use the 'Phospho Analysis Environment' kernel in your notebooks"

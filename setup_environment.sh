#!/bin/bash

# Cross-platform environment setup script for Phospho Analysis Project
# Works on both Linux and macOS

set -e  # Exit on any error

echo "Setting up Phospho Analysis Environment..."

# Detect operating system
if [[ "$OSTYPE" == "darwin"* ]]; then
    OS="macos"
    PYTHON_CMD="python3"
    PIP_CMD="pip3"
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
    OS="linux"
    PYTHON_CMD="python3"
    PIP_CMD="pip3"
else
    echo "Unsupported operating system: $OSTYPE"
    exit 1
fi

echo "Detected OS: $OS"

# Check if Python 3 is available
if ! command -v $PYTHON_CMD &> /dev/null; then
    echo "Error: Python 3 is not installed or not in PATH"
    exit 1
fi

# Check Python version
PYTHON_VERSION=$($PYTHON_CMD --version 2>&1 | cut -d' ' -f2 | cut -d'.' -f1,2)
echo "Python version: $PYTHON_VERSION"

# Install packages with appropriate flags for each OS
if [[ "$OS" == "macos" ]]; then
    echo "Installing packages for macOS..."
    $PIP_CMD install --user -r requirements.txt
elif [[ "$OS" == "linux" ]]; then
    echo "Installing packages for Linux..."
    $PIP_CMD install --break-system-packages -r requirements.txt
fi

# Register the kernel
echo "Registering Jupyter kernel..."
$PYTHON_CMD -m ipykernel install --user --name=phospho_env --display-name="Phospho Analysis Environment"

# Test installation
echo "Testing installation..."
$PYTHON_CMD -c "
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
echo ""
echo "To activate this environment in the future, just run:"
echo "  source setup_environment.sh"

#!/bin/bash

# Simple, portable environment setup script for Phospho Analysis Project
# Works on both macOS and Linux
# Usage: bash setup.sh

set -e  # Exit on any error

echo "=========================================="
echo "Phospho Analysis Project - Environment Setup"
echo "=========================================="
echo ""

# Check if Python 3 is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Error: Python 3 is not installed"
    echo "   Please install Python 3.8+ first"
    exit 1
fi

PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
echo "✓ Found Python: $PYTHON_VERSION"
echo ""

# Create virtual environment if it doesn't exist
if [ ! -d "lab_env" ]; then
    echo "Creating virtual environment 'lab_env'..."
    python3 -m venv lab_env
    echo "✓ Virtual environment created"
else
    echo "✓ Virtual environment 'lab_env' already exists"
fi
echo ""

# Use virtual environment's Python and pip directly
VENV_PYTHON="lab_env/bin/python"
VENV_PIP="lab_env/bin/pip"

# Upgrade pip
echo "Upgrading pip..."
$VENV_PIP install --upgrade pip --quiet

# Install requirements
echo "Installing Python packages from requirements.txt..."
$VENV_PIP install -r requirements.txt

echo ""
echo "Registering Jupyter kernel..."
$VENV_PYTHON -m ipykernel install --user --name=lab_env --display-name="Python (lab_env)"

echo ""
echo "=========================================="
echo "✓ Setup complete!"
echo "=========================================="
echo ""
echo "To use the environment:"
echo "  1. Activate: source lab_env/bin/activate"
echo "  2. Start Jupyter: jupyter lab"
echo "  3. Select kernel: 'Python (lab_env)' in your notebooks"
echo ""
echo "To deactivate: deactivate"
echo ""


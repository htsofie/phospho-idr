#!/usr/bin/env python3
"""
Script to convert Phosphomouse phosphorylation sites from Excel (.xlsb) to Parquet format.
"""

import pandas as pd
import os

def convert_xlsb_to_parquet(input_file, output_file=None):
    """
    Convert an Excel .xlsb file to parquet format.
    
    Args:
        input_file (str): Path to the input .xlsb file
        output_file (str): Path for the output .parquet file (optional)
    
    Returns:
        str: Path to the created parquet file
    """
    
    # Generate output filename if not provided
    if output_file is None:
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}.parquet"
    
    try:
        # Read the Excel .xlsb file using pyxlsb engine
        df = pd.read_excel(input_file, engine='pyxlsb')
        
        # Convert to parquet
        df.to_parquet(output_file, index=False)
        
        print(f"✅ Converted: {input_file} → {output_file}")
        print(f"   Data: {df.shape[0]} rows × {df.shape[1]} columns")
        
        return output_file
        
    except Exception as e:
        print(f"❌ Error: {e}")
        return None

def main():
    """Main function to run the conversion."""
    
    # Input file path
    input_file = "Phosphomouse_phosphorylation_sites.xlsb"
    
    # Check if input file exists
    if not os.path.exists(input_file):
        print(f"❌ File not found: {input_file}")
        return
    
    # Convert the file
    output_file = convert_xlsb_to_parquet(input_file)
    
    if output_file:
        print(f"✅ Conversion completed successfully!")
    else:
        print("❌ Conversion failed!")

if __name__ == "__main__":
    main()

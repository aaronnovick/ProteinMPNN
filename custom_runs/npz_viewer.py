#!/usr/bin/env python3
"""
NPZ File Viewer

This script opens and displays the contents of .npz files.
It shows detailed information about the data structure, shapes, and basic statistics.
"""

import numpy as np
import argparse
import os
import sys
from pathlib import Path

def view_npz_file(npz_path):
    """
    Open and display the contents of an NPZ file.
    
    Args:
        npz_path (str): Path to the NPZ file
    """
    try:
        # Load the NPZ file
        data = np.load(npz_path)
        
        print(f"\n{'='*60}")
        print(f"NPZ File: {npz_path}")
        print(f"{'='*60}")
        
        # Get file size
        file_size = os.path.getsize(npz_path)
        print(f"File size: {file_size:,} bytes ({file_size/1024:.2f} KB)")
        
        # List all arrays in the NPZ file
        print(f"\nArrays in file:")
        print(f"{'Array Name':<20} {'Shape':<15} {'Data Type':<12} {'Size':<10}")
        print(f"{'-'*20} {'-'*15} {'-'*12} {'-'*10}")
        
        total_elements = 0
        for key in data.keys():
            array = data[key]
            shape_str = str(array.shape)
            dtype_str = str(array.dtype)
            size = array.size
            total_elements += size
            
            print(f"{key:<20} {shape_str:<15} {dtype_str:<12} {size:,}")
        
        print(f"\nTotal elements: {total_elements:,}")
        
        # Display detailed information for each array
        print(f"\n{'='*60}")
        print("DETAILED ARRAY INFORMATION")
        print(f"{'='*60}")
        
        for key in data.keys():
            array = data[key]
            print(f"\nArray: '{key}'")
            print(f"{'-'*40}")
            print(f"Shape: {array.shape}")
            print(f"Data type: {array.dtype}")
            print(f"Size: {array.size:,} elements")
            print(f"Memory usage: {array.nbytes:,} bytes ({array.nbytes/1024:.2f} KB)")
            
            # Basic statistics for numeric arrays
            if np.issubdtype(array.dtype, np.number):
                print(f"Min: {np.min(array):.6f}")
                print(f"Max: {np.max(array):.6f}")
                print(f"Mean: {np.mean(array):.6f}")
                print(f"Std: {np.std(array):.6f}")
                print(f"Median: {np.median(array):.6f}")
                
                # Show first few and last few elements
                if array.size <= 10:
                    print(f"All values: {array}")
                else:
                    print(f"First 5 values: {array.flatten()[:5]}")
                    print(f"Last 5 values: {array.flatten()[-5:]}")
            else:
                # For non-numeric arrays, show unique values
                unique_vals = np.unique(array)
                if len(unique_vals) <= 10:
                    print(f"Unique values: {unique_vals}")
                else:
                    print(f"First 5 unique values: {unique_vals[:5]}")
                    print(f"Total unique values: {len(unique_vals)}")
        
        # Close the file
        data.close()
        
    except FileNotFoundError:
        print(f"Error: File '{npz_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading NPZ file: {e}")
        sys.exit(1)

def main():
    """Main function to handle command line arguments and execute the viewer."""
    parser = argparse.ArgumentParser(
        description="View the contents of NPZ files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 npz_viewer.py file.npz
  python3 npz_viewer.py /path/to/scores.npz
  python3 npz_viewer.py *.npz
        """
    )
    
    parser.add_argument(
        'npz_files',
        nargs='+',
        help='NPZ file(s) to view'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Show more detailed information'
    )
    
    args = parser.parse_args()
    
    # Process each NPZ file
    for npz_file in args.npz_files:
        # Handle glob patterns
        if '*' in npz_file or '?' in npz_file:
            import glob
            matching_files = glob.glob(npz_file)
            if not matching_files:
                print(f"No files found matching pattern: {npz_file}")
                continue
            for file in matching_files:
                if file.endswith('.npz'):
                    view_npz_file(file)
        else:
            if not npz_file.endswith('.npz'):
                print(f"Warning: File '{npz_file}' doesn't have .npz extension")
            view_npz_file(npz_file)

if __name__ == "__main__":
    main() 
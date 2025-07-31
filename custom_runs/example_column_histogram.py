#!/usr/bin/env python3
"""
Example script showing how to use the column histogram analysis.

This script demonstrates how to run the column histogram analysis
with different parameters and layouts.
"""

import subprocess
import sys
import os

def run_column_histogram_analysis():
    """Run the column histogram analysis with example parameters."""
    
    # Example 1: Generate only column layout with 20 bins
    print("Example 1: Column layout with 20 bins")
    cmd1 = [
        "python3", "custom_runs/column_histogram_analysis.py",
        "--sample_dirs",
        "outputs/my_variants/5UOI/sample_scores/5UOI_sample_1_score_only_from_fasta",
        "outputs/my_variants/5UOI/sample_scores/5UOI_sample_2_score_only_from_fasta",
        "outputs/my_variants/5UOI/sample_scores/5UOI_sample_3_score_only_from_fasta",
        "--output_dir", "outputs/my_variants/5UOI/example_analysis_1",
        "--bins", "20",
        "--layout", "column"
    ]
    
    try:
        subprocess.run(cmd1, check=True)
        print("✓ Example 1 completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"✗ Example 1 failed: {e}")
    
    # Example 2: Generate only grid layout with mean scores
    print("\nExample 2: Grid layout with mean scores")
    cmd2 = [
        "python3", "custom_runs/column_histogram_analysis.py",
        "--sample_dirs",
        "outputs/my_variants/5UOI/sample_scores/5UOI_sample_1_score_only_from_fasta",
        "outputs/my_variants/5UOI/sample_scores/5UOI_sample_2_score_only_from_fasta",
        "--output_dir", "outputs/my_variants/5UOI/example_analysis_2",
        "--bins", "15",
        "--layout", "grid",
        "--use_mean_scores"
    ]
    
    try:
        subprocess.run(cmd2, check=True)
        print("✓ Example 2 completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"✗ Example 2 failed: {e}")
    
    # Example 3: Generate both layouts with all samples
    print("\nExample 3: Both layouts with all samples")
    cmd3 = [
        "python3", "custom_runs/column_histogram_analysis.py",
        "--sample_dirs",
        "outputs/my_variants/5UOI/sample_scores/5UOI_sample_1_score_only_from_fasta",
        "outputs/my_variants/5UOI/sample_scores/5UOI_sample_2_score_only_from_fasta",
        "outputs/my_variants/5UOI/sample_scores/5UOI_sample_3_score_only_from_fasta",
        "outputs/my_variants/5UOI/sample_scores/5UOI_sample_4_score_only_from_fasta",
        "outputs/my_variants/5UOI/sample_scores/5UOI_sample_5_score_only_from_fasta",
        "--output_dir", "outputs/my_variants/5UOI/example_analysis_3",
        "--bins", "25",
        "--layout", "both"
    ]
    
    try:
        subprocess.run(cmd3, check=True)
        print("✓ Example 3 completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"✗ Example 3 failed: {e}")

def print_usage_examples():
    """Print usage examples for the column histogram analysis."""
    
    print("Column Histogram Analysis - Usage Examples")
    print("=" * 50)
    
    print("\n1. Basic usage (column layout only):")
    print("python3 custom_runs/column_histogram_analysis.py \\")
    print("  --sample_dirs sample_dir1 sample_dir2 \\")
    print("  --output_dir output_directory \\")
    print("  --layout column")
    
    print("\n2. Grid layout with mean scores:")
    print("python3 custom_runs/column_histogram_analysis.py \\")
    print("  --sample_dirs sample_dir1 sample_dir2 \\")
    print("  --output_dir output_directory \\")
    print("  --layout grid \\")
    print("  --use_mean_scores")
    
    print("\n3. Both layouts with custom bins:")
    print("python3 custom_runs/column_histogram_analysis.py \\")
    print("  --sample_dirs sample_dir1 sample_dir2 sample_dir3 \\")
    print("  --output_dir output_directory \\")
    print("  --layout both \\")
    print("  --bins 40")
    
    print("\n4. Run with default settings (all 5UOI samples):")
    print("python3 custom_runs/column_histogram_analysis.py")

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--help":
        print_usage_examples()
    else:
        print("Running column histogram analysis examples...")
        run_column_histogram_analysis()
        print("\nExamples completed! Check the output directories for results.") 
#!/usr/bin/env python3
"""
Overlayed Histogram Analysis for ProteinMPNN Variants

This script creates overlayed histograms of mean score distributions for all discovered sample variants
using NPZ files as input. Each NPZ file contains multiple scores for a single variant,
and this script calculates the mean score for each variant, then displays the distribution
of these mean scores as overlayed histograms for easy comparison.

The script automatically discovers all sample directories in sample_variant_scores,
so no manual configuration is needed for new samples.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import argparse
from pathlib import Path

def load_mean_scores_from_npz_files(npz_dir, exclude_files=None):
    """
    Load mean scores from NPZ files in a directory.
    Each NPZ file contains multiple scores for the same sequence variant.
    This function calculates the mean score for each variant.
    Excludes specified files (default: 5UOI_pdb.npz).
    
    Args:
        npz_dir (str): Directory containing NPZ files
        exclude_files (list): List of filenames to exclude (default: ['5UOI_pdb.npz'])
        
    Returns:
        list: List of mean scores, one per NPZ file
    """
    if exclude_files is None:
        exclude_files = ['5UOI_pdb.npz']
        
    mean_scores = []
    npz_files = glob.glob(os.path.join(npz_dir, "*.npz"))
    
    for npz_file in sorted(npz_files):
        # Skip excluded files
        if os.path.basename(npz_file) in exclude_files:
            print(f"Skipping excluded file: {npz_file}")
            continue
            
        try:
            data = np.load(npz_file)
            if 'score' in data:
                scores = data['score']
                # Calculate mean score for this variant
                mean_score = np.mean(scores)
                mean_scores.append(mean_score)
            data.close()
        except Exception as e:
            print(f"Warning: Could not load {npz_file}: {e}")
    
    return mean_scores

def create_overlayed_histograms(sample_data_dict, output_path, title="Score Distribution Comparison (Overlayed Histograms)"):
    """
    Create overlayed histograms for multiple samples.
    
    Args:
        sample_data_dict (dict): Dictionary with sample names as keys and score lists as values
        output_path (str): Path to save the plot
        title (str): Plot title
    """
    plt.figure(figsize=(12, 8))
    
    # Define colors for each sample
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    # Create overlayed histograms
    for i, (sample_name, scores) in enumerate(sample_data_dict.items()):
        if not scores:  # Skip empty samples
            continue
            
        color = colors[i % len(colors)]
        alpha = 0.6
        
        # Create histogram
        plt.hist(scores, bins=30, alpha=alpha, label=sample_name, 
                color=color, edgecolor='black', linewidth=0.5)
        
        # Add vertical line for mean
        mean_score = np.mean(scores)
        plt.axvline(mean_score, color=color, linestyle='--', linewidth=2, 
                   alpha=0.8, label=f'{sample_name} Mean: {mean_score:.4f}')
    
    plt.xlabel('Score', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    
    # Adjust layout to accommodate legend
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Overlayed histogram saved to: {output_path}")

def print_statistics(sample_data_dict):
    """
    Print statistical summary for all samples.
    
    Args:
        sample_data_dict (dict): Dictionary with sample data
    """
    print("\n" + "="*60)
    print("STATISTICAL SUMMARY")
    print("="*60)
    
    for sample_name, scores in sample_data_dict.items():
        if not scores:
            continue
            
        mean_score = np.mean(scores)
        std_score = np.std(scores)
        median_score = np.median(scores)
        min_score = np.min(scores)
        max_score = np.max(scores)
        
        print(f"\n{sample_name}:")
        print(f"  Count: {len(scores):,}")
        print(f"  Mean: {mean_score:.6f}")
        print(f"  Std: {std_score:.6f}")
        print(f"  Median: {median_score:.6f}")
        print(f"  Min: {min_score:.6f}")
        print(f"  Max: {max_score:.6f}")

def main():
    """Main function to handle command line arguments and execute the analysis."""
    parser = argparse.ArgumentParser(
        description='Create overlayed histograms of mean score distributions from NPZ files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 overlayed_histogram_analysis.py --output_dir output
  python3 overlayed_histogram_analysis.py --sample_dirs dir1 dir2 dir3 dir4 dir5 --output_dir output
  python3 overlayed_histogram_analysis.py --npz_patterns "sample_1/*.npz" "sample_2/*.npz" --output_dir output
        """
    )
    
    parser.add_argument(
        '--sample_dirs',
        nargs='+',
        help='Directories containing NPZ files for each sample'
    )
    
    parser.add_argument(
        '--npz_patterns',
        nargs='+',
        help='Glob patterns for NPZ files (alternative to sample_dirs)'
    )
    
    parser.add_argument(
        '--output_dir',
        default='../outputs/my_variants/5UOI/overlayed_histogram_analysis',
        help='Output directory for plots (default: ../outputs/my_variants/5UOI/overlayed_histogram_analysis)'
    )
    
    parser.add_argument(
        '--title',
        default='Score Distribution Comparison (Overlayed Histograms)',
        help='Title for the overlayed histogram'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    sample_data = {}
    
    if args.sample_dirs:
        # Load data from sample directories
        for i, sample_dir in enumerate(args.sample_dirs):
            sample_name = f"Sample_{i+1}"  # Use numbered names instead of directory basename
            scores = load_mean_scores_from_npz_files(sample_dir)
            sample_data[sample_name] = scores
            print(f"Loaded {len(scores)} mean scores from {sample_name}")
    
    elif args.npz_patterns:
        # Load data from NPZ patterns
        for i, pattern in enumerate(args.npz_patterns):
            sample_name = f"Sample_{i+1}"
            npz_files = glob.glob(pattern)
            
            mean_scores = []
            for npz_file in sorted(npz_files):
                try:
                    data = np.load(npz_file)
                    if 'score' in data:
                        scores = data['score']
                        # Calculate mean score for this variant
                        mean_score = np.mean(scores)
                        mean_scores.append(mean_score)
                    data.close()
                except Exception as e:
                    print(f"Warning: Could not load {npz_file}: {e}")
            
            sample_data[sample_name] = mean_scores
            print(f"Loaded {len(mean_scores)} mean scores from {sample_name} ({len(npz_files)} files)")
    
    else:
        # Default: Automatically discover all sample directories in sample_variant_scores
        base_dir = "../outputs/my_variants/5UOI/sample_variant_scores/5score"
        sample_dirs = []
        
        if os.path.exists(base_dir):
            # Get all subdirectories that contain 'score_only' folders
            for item in sorted(os.listdir(base_dir)):
                item_path = os.path.join(base_dir, item)
                score_only_path = os.path.join(item_path, "score_only")
                
                # Check if this is a directory and contains a score_only subdirectory
                if os.path.isdir(item_path) and os.path.isdir(score_only_path):
                    # Check if it's a numeric sample directory
                    if item.isdigit():
                        sample_dirs.append(score_only_path)
                        print(f"Found sample directory: {item}")
                    else:
                        print(f"Skipping non-numeric directory: {item}")
            
            print(f"Discovered {len(sample_dirs)} sample directories")
        else:
            print(f"Warning: Base directory {base_dir} does not exist")
            sample_dirs = []
        
        # Process all discovered sample directories
        for i, sample_dir in enumerate(sample_dirs):
            if os.path.exists(sample_dir):
                # Extract sample number from the path for better naming
                sample_number = os.path.basename(os.path.dirname(sample_dir))
                sample_name = f"Sample_{sample_number}"
                scores = load_mean_scores_from_npz_files(sample_dir)
                sample_data[sample_name] = scores
                print(f"Loaded {len(scores)} mean scores from {sample_name}")
            else:
                print(f"Warning: Sample directory does not exist: {sample_dir}")
        
        # Sort sample_data by sample number to ensure proper ordering in legend
        sorted_sample_data = {}
        for sample_name in sorted(sample_data.keys(), key=lambda x: int(x.split('_')[1])):
            sorted_sample_data[sample_name] = sample_data[sample_name]
        
        sample_data = sorted_sample_data
    
    if not sample_data:
        print("Error: No data loaded. Please check your input directories or patterns.")
        return
    
    # Create overlayed histogram
    output_path = os.path.join(args.output_dir, "overlayed_histogram.png")
    create_overlayed_histograms(sample_data, output_path, args.title)
    
    # Print statistics
    print_statistics(sample_data)

if __name__ == "__main__":
    main() 
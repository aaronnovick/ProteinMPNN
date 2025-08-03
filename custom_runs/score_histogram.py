#!/usr/bin/env python3
"""
Score Histogram Generator for ProteinMPNN Variants

This script generates histogram plots of scores for sequence variants.
It can read from CSV summary files or individual NPZ files.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from pathlib import Path
import argparse

def load_scores_from_csv(csv_path):
    """
    Load scores from CSV summary file.
    
    Args:
        csv_path (str): Path to CSV file containing score summaries
        
    Returns:
        list: List of mean scores for each variant
    """
    df = pd.read_csv(csv_path)
    return df['Mean Score'].tolist()

def load_mean_scores_from_npz_files(npz_dir):
    """
    Load mean scores from NPZ files in a directory.
    Each NPZ file contains multiple scores for the same sequence variant.
    This function calculates the mean score for each variant.
    
    Args:
        npz_dir (str): Directory containing NPZ files
        
    Returns:
        list: List of mean scores, one per NPZ file
    """
    mean_scores = []
    npz_files = glob.glob(os.path.join(npz_dir, "*.npz"))
    
    for npz_file in npz_files:
        data = np.load(npz_file)
        scores = data['score']
        # Calculate mean score for this variant
        mean_score = np.mean(scores)
        mean_scores.append(mean_score)
    
    return mean_scores

def create_histogram(scores, output_path, title="Score Distribution", 
                    xlabel="Score", ylabel="Frequency", bins=30, 
                    color='skyblue', edgecolor='black', alpha=0.7):
    """
    Create and save a histogram plot of scores.
    
    Args:
        scores (list): List of scores to plot
        output_path (str): Path to save the plot
        title (str): Plot title
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        bins (int): Number of histogram bins
        color (str): Histogram color
        edgecolor (str): Histogram edge color
        alpha (float): Transparency level
    """
    plt.figure(figsize=(10, 6))
    
    # Create histogram
    plt.hist(scores, bins=bins, color=color, edgecolor=edgecolor, alpha=alpha)
    
    # Add statistics
    mean_score = np.mean(scores)
    std_score = np.std(scores)
    median_score = np.median(scores)
    
    # Add vertical lines for mean and median
    plt.axvline(mean_score, color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {mean_score:.4f}')
    plt.axvline(median_score, color='green', linestyle='--', linewidth=2, 
                label=f'Median: {median_score:.4f}')
    
    # Add text box with statistics
    stats_text = f'Mean: {mean_score:.4f}\nStd: {std_score:.4f}\nMedian: {median_score:.4f}\nN: {len(scores)}'
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Customize plot
    plt.title(title, fontsize=14, fontweight='bold')
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Histogram saved to: {output_path}")
    print(f"Statistics: Mean={mean_score:.4f}, Std={std_score:.4f}, Median={median_score:.4f}, N={len(scores)}")

def main():
    parser = argparse.ArgumentParser(description='Generate histogram plots of sequence variant scores')
    parser.add_argument('--input', required=True, 
                       help='Input path: CSV file or directory containing NPZ files')
    parser.add_argument('--output', required=True,
                       help='Output path for the histogram plot (e.g., histogram.png)')
    parser.add_argument('--title', default='Score Distribution',
                       help='Title for the histogram plot')
    parser.add_argument('--bins', type=int, default=30,
                       help='Number of histogram bins (default: 30)')
    parser.add_argument('--data_type', choices=['csv', 'npz'], default='auto',
                       help='Data type: csv for summary file, npz for individual files (default: auto-detect)')
    
    args = parser.parse_args()
    
    # Auto-detect data type if not specified
    if args.data_type == 'auto':
        if args.input.endswith('.csv'):
            args.data_type = 'csv'
        elif os.path.isdir(args.input):
            args.data_type = 'npz'
        else:
            raise ValueError("Could not auto-detect data type. Please specify --data_type")
    
    # Load scores based on data type
    if args.data_type == 'csv':
        scores = load_scores_from_csv(args.input)
        print(f"Loaded {len(scores)} mean scores from CSV file")
    elif args.data_type == 'npz':
        scores = load_mean_scores_from_npz_files(args.input)
        print(f"Loaded {len(scores)} mean scores from NPZ files")
    
    # Create histogram
    create_histogram(scores, args.output, title=args.title, bins=args.bins)

if __name__ == "__main__":
    # Example usage without command line arguments
    if len(os.sys.argv) == 1:
        # Default example: plot scores from NPZ files in sample 1
        npz_dir = "../outputs/my_variants/5UOI/sample_variant_scores/5/score_only"
        output_path = "../outputs/my_variants/5UOI/sample_variant_scores/5/score_only/score_histogram.png"
        
        if os.path.exists(npz_dir):
            print("Creating histogram from NPZ files...")
            scores = load_mean_scores_from_npz_files(npz_dir)
            create_histogram(scores, output_path, 
                           title="5UOI Sample 5 Variants - Score Distribution",
                           xlabel="Score", ylabel="Frequency")
        else:
            print(f"NPZ directory not found: {npz_dir}")
            print("Please run with command line arguments or ensure the default directory exists.")
    else:
        main() 
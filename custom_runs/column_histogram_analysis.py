#!/usr/bin/env python3
"""
Column Histogram Analysis for ProteinMPNN Variants

This script generates histograms in a column layout with shared bins
for better comparison across multiple samples.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from pathlib import Path
import argparse

def load_sample_data(sample_dir):
    """
    Load all score data from a sample directory.
    
    Args:
        sample_dir (str): Path to sample directory
        
    Returns:
        dict: Dictionary containing scores and metadata
    """
    csv_path = os.path.join(sample_dir, "score_only", "score_summary.csv")
    npz_dir = os.path.join(sample_dir, "score_only")
    
    data = {}
    
    # Load CSV summary
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
        data['mean_scores'] = df['Mean Score'].tolist()
        data['std_scores'] = df['Std Dev'].tolist()
        data['sequences'] = df['Sequence'].tolist()
    
    # Load individual NPZ scores
    if os.path.exists(npz_dir):
        all_scores = []
        npz_files = glob.glob(os.path.join(npz_dir, "*.npz"))
        
        for npz_file in npz_files:
            file_data = np.load(npz_file)
            scores = file_data['score']
            all_scores.extend(scores)
        
        data['individual_scores'] = all_scores
    
    return data

def create_column_histograms(sample_data_dict, output_path, bins=30, 
                           figsize=(10, 12), use_individual_scores=True):
    """
    Create histograms in a column layout with shared bins.
    
    Args:
        sample_data_dict (dict): Dictionary with sample names as keys and data as values
        output_path (str): Path to save the plot
        bins (int): Number of histogram bins
        figsize (tuple): Figure size (width, height)
        use_individual_scores (bool): Whether to use individual scores or mean scores
    """
    # Collect all scores to determine shared bin range
    all_scores = []
    sample_names = []
    sample_scores = {}
    
    for sample_name, data in sample_data_dict.items():
        if use_individual_scores and 'individual_scores' in data:
            scores = data['individual_scores']
        elif 'mean_scores' in data:
            scores = data['mean_scores']
        else:
            continue
            
        all_scores.extend(scores)
        sample_names.append(sample_name)
        sample_scores[sample_name] = scores
    
    if not all_scores:
        print("No score data found!")
        return
    
    # Determine shared bin edges
    min_score = min(all_scores)
    max_score = max(all_scores)
    bin_edges = np.linspace(min_score, max_score, bins + 1)
    
    # Create subplots
    n_samples = len(sample_names)
    fig, axes = plt.subplots(n_samples, 1, figsize=figsize, sharex=True)
    
    # If only one sample, make axes iterable
    if n_samples == 1:
        axes = [axes]
    
    # Create histograms for each sample
    for i, sample_name in enumerate(sample_names):
        scores = sample_scores[sample_name]
        
        # Create histogram
        n, bins_hist, patches = axes[i].hist(scores, bins=bin_edges, 
                                            color='skyblue', edgecolor='black', 
                                            alpha=0.7, density=True)
        
        # Add statistics
        mean_score = np.mean(scores)
        std_score = np.std(scores)
        median_score = np.median(scores)
        
        # Add vertical lines for mean and median
        axes[i].axvline(mean_score, color='red', linestyle='--', linewidth=2, 
                       label=f'Mean: {mean_score:.4f}')
        axes[i].axvline(median_score, color='green', linestyle='--', linewidth=2, 
                       label=f'Median: {median_score:.4f}')
        
        # Add text box with statistics
        stats_text = f'{sample_name}\nMean: {mean_score:.4f}\nStd: {std_score:.4f}\nN: {len(scores)}'
        axes[i].text(0.02, 0.95, stats_text, transform=axes[i].transAxes, 
                    verticalalignment='top', bbox=dict(boxstyle='round', 
                    facecolor='white', alpha=0.8), fontsize=10)
        
        # Customize subplot
        axes[i].set_ylabel('Density', fontsize=10)
        axes[i].legend(fontsize=8)
        axes[i].grid(True, alpha=0.3)
        
        # Only show legend for the first subplot to avoid clutter
        if i > 0:
            axes[i].legend().remove()
    
    # Set common x-axis label
    axes[-1].set_xlabel('Score', fontsize=12)
    
    # Add overall title
    fig.suptitle('Score Distribution Comparison (Column Layout)', fontsize=14, fontweight='bold')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Column histogram plot saved to: {output_path}")

def create_shared_bin_histograms(sample_data_dict, output_path, bins=30, 
                                figsize=(12, 8), use_individual_scores=True):
    """
    Create histograms with shared bins in a grid layout.
    
    Args:
        sample_data_dict (dict): Dictionary with sample names as keys and data as values
        output_path (str): Path to save the plot
        bins (int): Number of histogram bins
        figsize (tuple): Figure size (width, height)
        use_individual_scores (bool): Whether to use individual scores or mean scores
    """
    # Collect all scores to determine shared bin range
    all_scores = []
    sample_names = []
    sample_scores = {}
    
    for sample_name, data in sample_data_dict.items():
        if use_individual_scores and 'individual_scores' in data:
            scores = data['individual_scores']
        elif 'mean_scores' in data:
            scores = data['mean_scores']
        else:
            continue
            
        all_scores.extend(scores)
        sample_names.append(sample_name)
        sample_scores[sample_name] = scores
    
    if not all_scores:
        print("No score data found!")
        return
    
    # Determine shared bin edges
    min_score = min(all_scores)
    max_score = max(all_scores)
    bin_edges = np.linspace(min_score, max_score, bins + 1)
    
    # Calculate grid dimensions
    n_samples = len(sample_names)
    n_cols = min(3, n_samples)  # Max 3 columns
    n_rows = (n_samples + n_cols - 1) // n_cols
    
    # Create subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, sharex=True, sharey=True)
    
    # Make axes iterable if only one subplot
    if n_samples == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    
    # Create histograms for each sample
    for i, sample_name in enumerate(sample_names):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]
        
        scores = sample_scores[sample_name]
        
        # Create histogram
        n, bins_hist, patches = ax.hist(scores, bins=bin_edges, 
                                       color='skyblue', edgecolor='black', 
                                       alpha=0.7, density=True)
        
        # Add statistics
        mean_score = np.mean(scores)
        std_score = np.std(scores)
        median_score = np.median(scores)
        
        # Add vertical lines for mean and median
        ax.axvline(mean_score, color='red', linestyle='--', linewidth=2, 
                  label=f'Mean: {mean_score:.4f}')
        ax.axvline(median_score, color='green', linestyle='--', linewidth=2, 
                  label=f'Median: {median_score:.4f}')
        
        # Add title with statistics
        title = f'{sample_name}\nMean: {mean_score:.4f}, Std: {std_score:.4f}, N: {len(scores)}'
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    # Hide empty subplots
    for i in range(n_samples, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)
    
    # Set common labels
    fig.text(0.5, 0.02, 'Score', ha='center', fontsize=12)
    fig.text(0.02, 0.5, 'Density', va='center', rotation='vertical', fontsize=12)
    
    # Add overall title
    fig.suptitle('Score Distribution Comparison (Grid Layout)', fontsize=14, fontweight='bold')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Grid histogram plot saved to: {output_path}")

def statistical_analysis(sample_data_dict, use_individual_scores=True):
    """
    Perform statistical analysis on sample data.
    
    Args:
        sample_data_dict (dict): Dictionary with sample data
        use_individual_scores (bool): Whether to use individual scores or mean scores
        
    Returns:
        dict: Statistical summary
    """
    stats_summary = {}
    
    for sample_name, data in sample_data_dict.items():
        if use_individual_scores and 'individual_scores' in data:
            scores = data['individual_scores']
        elif 'mean_scores' in data:
            scores = data['mean_scores']
        else:
            continue
            
        stats_summary[sample_name] = {
            'mean': np.mean(scores),
            'std': np.std(scores),
            'median': np.median(scores),
            'min': np.min(scores),
            'max': np.max(scores),
            'count': len(scores)
        }
    
    return stats_summary

def main():
    parser = argparse.ArgumentParser(description='Column histogram analysis of sequence variant scores')
    parser.add_argument('--sample_dirs', nargs='+', required=True,
                       help='Directories containing sample data')
    parser.add_argument('--output_dir', required=True,
                       help='Output directory for plots')
    parser.add_argument('--bins', type=int, default=30,
                       help='Number of histogram bins (default: 30)')
    parser.add_argument('--use_mean_scores', action='store_true',
                       help='Use mean scores instead of individual scores')
    parser.add_argument('--layout', choices=['column', 'grid', 'both'], default='both',
                       help='Layout type for histograms')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data from all samples
    sample_data = {}
    for sample_dir in args.sample_dirs:
        sample_name = os.path.basename(sample_dir)
        sample_data[sample_name] = load_sample_data(sample_dir)
    
    use_individual = not args.use_mean_scores
    
    # Generate plots based on layout choice
    if args.layout in ['column', 'both']:
        output_path = os.path.join(args.output_dir, "column_histograms.png")
        create_column_histograms(sample_data, output_path, bins=args.bins, 
                               use_individual_scores=use_individual)
    
    if args.layout in ['grid', 'both']:
        output_path = os.path.join(args.output_dir, "grid_histograms.png")
        create_shared_bin_histograms(sample_data, output_path, bins=args.bins, 
                                   use_individual_scores=use_individual)
    
    # Print statistical summary
    stats_summary = statistical_analysis(sample_data, use_individual_scores=use_individual)
    print("\nStatistical Summary:")
    print("=" * 50)
    for sample_name, stats in stats_summary.items():
        print(f"\n{sample_name}:")
        for stat_name, value in stats.items():
            print(f"  {stat_name}: {value:.4f}")

if __name__ == "__main__":
    # Example usage without command line arguments
    if len(os.sys.argv) == 1:
        # Default example: analyze all 5UOI samples
        sample_dirs = [
            "outputs/my_variants/5UOI/sample_scores/5UOI_sample_1_score_only_from_fasta",
            "outputs/my_variants/5UOI/sample_scores/5UOI_sample_2_score_only_from_fasta",
            "outputs/my_variants/5UOI/sample_scores/5UOI_sample_3_score_only_from_fasta",
            "outputs/my_variants/5UOI/sample_scores/5UOI_sample_4_score_only_from_fasta",
            "outputs/my_variants/5UOI/sample_scores/5UOI_sample_5_score_only_from_fasta"
        ]
        
        output_dir = "outputs/my_variants/5UOI/column_histogram_analysis"
        
        # Check which sample directories exist
        existing_dirs = [d for d in sample_dirs if os.path.exists(d)]
        
        if existing_dirs:
            print(f"Analyzing {len(existing_dirs)} sample directories...")
            
            # Load data
            sample_data = {}
            for sample_dir in existing_dirs:
                sample_name = os.path.basename(sample_dir)
                sample_data[sample_name] = load_sample_data(sample_dir)
            
            # Create output directory
            os.makedirs(output_dir, exist_ok=True)
            
            # Generate column histograms
            output_path = os.path.join(output_dir, "column_histograms.png")
            create_column_histograms(sample_data, output_path)
            
            # Generate grid histograms
            output_path = os.path.join(output_dir, "grid_histograms.png")
            create_shared_bin_histograms(sample_data, output_path)
            
            # Print statistics
            stats_summary = statistical_analysis(sample_data)
            print("\nStatistical Summary:")
            print("=" * 50)
            for sample_name, stats in stats_summary.items():
                print(f"\n{sample_name}:")
                for stat_name, value in stats.items():
                    print(f"  {stat_name}: {value:.4f}")
        else:
            print("No sample directories found. Please run with command line arguments.")
    else:
        main() 
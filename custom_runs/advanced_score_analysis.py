#!/usr/bin/env python3
"""
Advanced Score Analysis for ProteinMPNN Variants

This script provides comprehensive analysis and visualization of sequence variant scores,
including histograms, box plots, and statistical comparisons across multiple samples.
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from pathlib import Path
import argparse
import seaborn as sns
from scipy import stats

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

def create_comparison_plot(sample_data_dict, output_path, plot_type='histogram'):
    """
    Create comparison plots for multiple samples.
    
    Args:
        sample_data_dict (dict): Dictionary with sample names as keys and data as values
        output_path (str): Path to save the plot
        plot_type (str): Type of plot ('histogram', 'boxplot', 'violin')
    """
    plt.figure(figsize=(12, 8))
    
    if plot_type == 'histogram':
        for sample_name, data in sample_data_dict.items():
            if 'individual_scores' in data:
                plt.hist(data['individual_scores'], alpha=0.6, label=sample_name, bins=30)
            elif 'mean_scores' in data:
                plt.hist(data['mean_scores'], alpha=0.6, label=sample_name, bins=30)
        
        plt.xlabel('Score')
        plt.ylabel('Frequency')
        plt.title('Score Distribution Comparison')
        plt.legend()
        
    elif plot_type == 'boxplot':
        plot_data = []
        labels = []
        
        for sample_name, data in sample_data_dict.items():
            if 'individual_scores' in data:
                plot_data.append(data['individual_scores'])
                labels.append(sample_name)
            elif 'mean_scores' in data:
                plot_data.append(data['mean_scores'])
                labels.append(sample_name)
        
        plt.boxplot(plot_data, labels=labels)
        plt.ylabel('Score')
        plt.title('Score Distribution Comparison (Box Plot)')
        plt.xticks(rotation=45)
        
    elif plot_type == 'violin':
        plot_data = []
        labels = []
        
        for sample_name, data in sample_data_dict.items():
            if 'individual_scores' in data:
                plot_data.append(data['individual_scores'])
                labels.append(sample_name)
            elif 'mean_scores' in data:
                plot_data.append(data['mean_scores'])
                labels.append(sample_name)
        
        plt.violinplot(plot_data, labels=labels)
        plt.ylabel('Score')
        plt.title('Score Distribution Comparison (Violin Plot)')
        plt.xticks(rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Comparison plot saved to: {output_path}")

def statistical_analysis(sample_data_dict):
    """
    Perform statistical analysis on sample data.
    
    Args:
        sample_data_dict (dict): Dictionary with sample data
        
    Returns:
        dict: Statistical summary
    """
    stats_summary = {}
    
    for sample_name, data in sample_data_dict.items():
        if 'individual_scores' in data:
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

def create_detailed_histogram(scores, output_path, title="Score Distribution", 
                            xlabel="Score", ylabel="Frequency", bins=30):
    """
    Create a detailed histogram with additional statistical information.
    
    Args:
        scores (list): List of scores to plot
        output_path (str): Path to save the plot
        title (str): Plot title
        xlabel (str): X-axis label
        ylabel (str): Y-axis label
        bins (int): Number of histogram bins
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), 
                                   gridspec_kw={'height_ratios': [3, 1]})
    
    # Main histogram
    n, bins_edges, patches = ax1.hist(scores, bins=bins, color='skyblue', 
                                      edgecolor='black', alpha=0.7)
    
    # Statistics
    mean_score = np.mean(scores)
    std_score = np.std(scores)
    median_score = np.median(scores)
    
    # Add vertical lines
    ax1.axvline(mean_score, color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {mean_score:.4f}')
    ax1.axvline(median_score, color='green', linestyle='--', linewidth=2, 
                label=f'Median: {median_score:.4f}')
    
    # Add text box
    stats_text = f'Mean: {mean_score:.4f}\nStd: {std_score:.4f}\nMedian: {median_score:.4f}\nN: {len(scores)}'
    ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax1.set_title(title, fontsize=14, fontweight='bold')
    ax1.set_ylabel(ylabel, fontsize=12)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Cumulative distribution
    sorted_scores = np.sort(scores)
    cumulative = np.arange(1, len(sorted_scores) + 1) / len(sorted_scores)
    ax2.plot(sorted_scores, cumulative, 'b-', linewidth=2)
    ax2.set_xlabel(xlabel, fontsize=12)
    ax2.set_ylabel('Cumulative Probability', fontsize=12)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Detailed histogram saved to: {output_path}")
    print(f"Statistics: Mean={mean_score:.4f}, Std={std_score:.4f}, Median={median_score:.4f}, N={len(scores)}")

def main():
    parser = argparse.ArgumentParser(description='Advanced analysis of sequence variant scores')
    parser.add_argument('--sample_dirs', nargs='+', required=True,
                       help='Directories containing sample data')
    parser.add_argument('--output_dir', required=True,
                       help='Output directory for plots')
    parser.add_argument('--plot_types', nargs='+', default=['histogram', 'boxplot'],
                       choices=['histogram', 'boxplot', 'violin', 'detailed'],
                       help='Types of plots to generate')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data from all samples
    sample_data = {}
    for sample_dir in args.sample_dirs:
        sample_name = os.path.basename(sample_dir)
        sample_data[sample_name] = load_sample_data(sample_dir)
    
    # Generate plots
    for plot_type in args.plot_types:
        if plot_type == 'detailed':
            # Create detailed histogram for each sample
            for sample_name, data in sample_data.items():
                if 'individual_scores' in data:
                    scores = data['individual_scores']
                elif 'mean_scores' in data:
                    scores = data['mean_scores']
                else:
                    continue
                
                output_path = os.path.join(args.output_dir, f"{sample_name}_detailed_histogram.png")
                create_detailed_histogram(scores, output_path, 
                                       title=f"{sample_name} - Score Distribution")
        else:
            # Create comparison plots
            output_path = os.path.join(args.output_dir, f"comparison_{plot_type}.png")
            create_comparison_plot(sample_data, output_path, plot_type)
    
    # Print statistical summary
    stats_summary = statistical_analysis(sample_data)
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
        
        output_dir = "outputs/my_variants/5UOI/score_analysis"
        
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
            
            # Generate detailed histograms for each sample
            for sample_name, data in sample_data.items():
                if 'individual_scores' in data:
                    scores = data['individual_scores']
                elif 'mean_scores' in data:
                    scores = data['mean_scores']
                else:
                    continue
                
                output_path = os.path.join(output_dir, f"{sample_name}_detailed_histogram.png")
                create_detailed_histogram(scores, output_path, 
                                       title=f"{sample_name} - Score Distribution")
            
            # Generate comparison plots
            if len(sample_data) > 1:
                for plot_type in ['histogram', 'boxplot']:
                    output_path = os.path.join(output_dir, f"comparison_{plot_type}.png")
                    create_comparison_plot(sample_data, output_path, plot_type)
            
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
#!/usr/bin/env python3
"""
Coefficient of Variation vs. Number of Scores Analysis for ProteinMPNN Variants

This script creates a plot showing CV vs. number of scores to evaluate the consistency
of each calculation, as required by the assignment.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import argparse

def load_variant_cv_data(npz_dir, exclude_files=None):
    """
    Load coefficient of variation data for each variant in a directory.
    
    Args:
        npz_dir (str): Directory containing NPZ files
        exclude_files (list): List of filenames to exclude (default: ['5UOI_pdb.npz'])
        
    Returns:
        list: List of tuples (cv, n_scores) for each variant
    """
    if exclude_files is None:
        exclude_files = ['5UOI_pdb.npz']
        
    variant_data = []
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
                n_scores = len(scores)
                
                if n_scores > 1:  # Need at least 2 scores to calculate CV
                    # Calculate coefficient of variation
                    mean_score = np.mean(scores)
                    std_score = np.std(scores)
                    
                    if mean_score != 0:  # Avoid division by zero
                        cv = std_score / abs(mean_score)
                        variant_data.append((cv, n_scores))
                
            data.close()
        except Exception as e:
            print(f"Warning: Could not load {npz_file}: {e}")
    
    return variant_data

def create_cv_vs_nscores_plot(sample_data_dict, output_path, title="CV vs. Number of Scores"):
    """
    Create a plot showing coefficient of variation vs. number of scores for all samples.
    This addresses the assignment requirement to evaluate CV vs. number of scores
    to determine the consistency of each calculation.
    
    Args:
        sample_data_dict (dict): Dictionary with sample names as keys and variant data as values
        output_path (str): Path to save the plot
        title (str): Plot title
    """
    plt.figure(figsize=(14, 10))
    
    # Define colors for each sample
    colors = plt.cm.tab20(np.linspace(0, 1, len(sample_data_dict)))
    
    # Create scatter plot for each sample
    for i, (sample_name, variant_data) in enumerate(sample_data_dict.items()):
        if not variant_data:  # Skip empty samples
            continue
            
        # Extract CV and n_scores data
        cvs = [data[0] for data in variant_data]
        n_scores = [data[1] for data in variant_data]
        
        color = colors[i % len(colors)]
        alpha = 0.7
        
        # Create scatter plot
        plt.scatter(n_scores, cvs, alpha=alpha, label=sample_name, 
                   color=color, s=30, edgecolors='black', linewidth=0.5)
        
        # Add trend line for this sample
        if len(n_scores) > 1:
            z = np.polyfit(n_scores, cvs, 1)
            p = np.poly1d(z)
            plt.plot(n_scores, p(n_scores), color=color, alpha=0.5, linestyle='--', linewidth=1)
    
    plt.xlabel('Number of Scores per Variant', fontsize=12)
    plt.ylabel('Coefficient of Variation (CV)', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    plt.grid(True, alpha=0.3)
    
    # Use log scale for better visualization if CV values vary widely
    if any(cv > 10 for sample_data in sample_data_dict.values() for cv, _ in sample_data):
        plt.yscale('log')
        plt.ylabel('Coefficient of Variation (CV) - Log Scale', fontsize=12)
    
    # Adjust layout to accommodate legend
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"CV vs. Number of Scores plot saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(
        description="Analyze coefficient of variation vs. number of scores for ProteinMPNN variants",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 cv_vs_nscores_analysis.py --output_dir output
  python3 cv_vs_nscores_analysis.py --sample_dirs dir1 dir2 dir3 --output_dir output
  python3 cv_vs_nscores_analysis.py --npz_patterns "sample_1/*.npz" "sample_2/*.npz" --output_dir output
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
        default='../outputs/my_variants/5UOI/cv_vs_nscores_analysis',
        help='Output directory for plots (default: ../outputs/my_variants/5UOI/cv_vs_nscores_analysis)'
    )
    
    parser.add_argument(
        '--title',
        default='CV vs. Number of Scores Analysis',
        help='Title for the plot'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    sample_data = {}
    
    if args.sample_dirs:
        # Load data from sample directories
        for i, sample_dir in enumerate(args.sample_dirs):
            sample_name = f"Sample_{i+1}"
            variant_data = load_variant_cv_data(sample_dir)
            sample_data[sample_name] = variant_data
            print(f"Loaded CV data for {len(variant_data)} variants from {sample_name}")
    
    elif args.npz_patterns:
        # Load data from NPZ patterns
        for i, pattern in enumerate(args.npz_patterns):
            sample_name = f"Sample_{i+1}"
            npz_files = glob.glob(pattern)
            
            variant_data = []
            for npz_file in sorted(npz_files):
                try:
                    data = np.load(npz_file)
                    if 'score' in data:
                        scores = data['score']
                        n_scores = len(scores)
                        
                        if n_scores > 1:
                            mean_score = np.mean(scores)
                            std_score = np.std(scores)
                            
                            if mean_score != 0:
                                cv = std_score / abs(mean_score)
                                variant_data.append((cv, n_scores))
                    
                    data.close()
                except Exception as e:
                    print(f"Warning: Could not load {npz_file}: {e}")
            
            sample_data[sample_name] = variant_data
            print(f"Loaded CV data for {len(variant_data)} variants from {sample_name} ({len(npz_files)} files)")
    
    else:
        # Default: Automatically discover all sample directories in sample_variant_scores
        base_dir = "../outputs/my_variants/5UOI/sample_variant_scores"
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
                variant_data = load_variant_cv_data(sample_dir)
                sample_data[sample_name] = variant_data
                print(f"Loaded CV data for {len(variant_data)} variants from {sample_name}")
            else:
                print(f"Warning: Sample directory does not exist: {sample_dir}")
        
        # Sort sample_data by sample number to ensure proper ordering
        sorted_sample_data = {}
        for sample_name in sorted(sample_data.keys(), key=lambda x: int(x.split('_')[1])):
            sorted_sample_data[sample_name] = sample_data[sample_name]
        
        sample_data = sorted_sample_data
    
    if not sample_data:
        print("Error: No data loaded. Please check your input directories or patterns.")
        return
    
    # Create CV vs. Number of Scores plot
    cv_plot_path = os.path.join(args.output_dir, "cv_vs_nscores.png")
    create_cv_vs_nscores_plot(sample_data, cv_plot_path, args.title)
    
    print(f"\nAnalysis complete! CV vs. Number of Scores plot saved to: {cv_plot_path}")

if __name__ == "__main__":
    main()

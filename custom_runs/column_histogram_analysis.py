#!/usr/bin/env python3
"""
Column Histogram Analysis for ProteinMPNN Variants

This script generates histograms in a column layout with shared bins
for better comparison across multiple samples. It uses mean scores from NPZ files,
where each NPZ file contains multiple scores for a single variant, and the script
calculates the mean score for each variant.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import glob
import argparse

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

def create_column_histograms(sample_data_dict, output_path, bins=30, 
                           figsize=(10, 12), use_individual_scores=True, samples_per_page=25):
    """
    Create histograms in a column layout with shared bins.
    If more than samples_per_page samples, creates multiple pages.
    
    Args:
        sample_data_dict (dict): Dictionary with sample names as keys and data as values
        output_path (str): Path to save the plot (base name, will append page numbers)
        bins (int): Number of histogram bins
        figsize (tuple): Figure size (width, height)
        use_individual_scores (bool): Whether to use individual scores or mean scores
        samples_per_page (int): Maximum number of samples per page
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
    
    # Split samples into pages
    total_samples = len(sample_names)
    num_pages = (total_samples + samples_per_page - 1) // samples_per_page
    
    # Open PDF file for writing
    with PdfPages(output_path) as pdf:
        for page in range(num_pages):
            start_idx = page * samples_per_page
            end_idx = min((page + 1) * samples_per_page, total_samples)
            page_samples = sample_names[start_idx:end_idx]
            
            # Create subplots for this page
            n_page_samples = len(page_samples)
            fig, axes = plt.subplots(n_page_samples, 1, figsize=figsize, sharex=True)
            
            # If only one sample, make axes iterable
            if n_page_samples == 1:
                axes = [axes]
            
            # Create histograms for each sample on this page
            for i, sample_name in enumerate(page_samples):
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
                axes[i].set_ylabel('Frequency', fontsize=10)
                axes[i].legend(fontsize=8)
                axes[i].grid(True, alpha=0.3)
            
            # Set common x-axis label
            axes[-1].set_xlabel('Score', fontsize=12)
            
            # Add overall title with page info
            if num_pages > 1:
                title = f'Score Distribution Comparison (Column Layout) - Page {page + 1} of {num_pages}'
            else:
                title = 'Score Distribution Comparison (Column Layout)'
            fig.suptitle(title, fontsize=14, fontweight='bold')
            
            # Adjust layout
            plt.tight_layout()
            
            # Save page to PDF
            pdf.savefig(fig, dpi=300, bbox_inches='tight')
            print(f"Column histogram page {page + 1} added to: {output_path}")
            
            plt.close()

def create_shared_bin_histograms(sample_data_dict, output_path, bins=30, 
                                figsize=(12, 8), use_individual_scores=True, samples_per_page=25):
    """
    Create histograms with shared bins in a grid layout.
    If more than samples_per_page samples, creates multiple pages.
    
    Args:
        sample_data_dict (dict): Dictionary with sample names as keys and data as values
        output_path (str): Path to save the plot (base name, will append page numbers)
        bins (int): Number of histogram bins
        figsize (tuple): Figure size (width, height)
        use_individual_scores (bool): Whether to use individual scores or mean scores
        samples_per_page (int): Maximum number of samples per page
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
    
    # Split samples into pages
    total_samples = len(sample_names)
    num_pages = (total_samples + samples_per_page - 1) // samples_per_page
    
    # Open PDF file for writing
    with PdfPages(output_path) as pdf:
        for page in range(num_pages):
            start_idx = page * samples_per_page
            end_idx = min((page + 1) * samples_per_page, total_samples)
            page_samples = sample_names[start_idx:end_idx]
            
            # Calculate grid dimensions for this page
            n_page_samples = len(page_samples)
            n_cols = min(3, n_page_samples)  # Max 3 columns
            n_rows = (n_page_samples + n_cols - 1) // n_cols
            
            # Create subplots for this page
            fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize, sharex=True, sharey=True)
            
            # Make axes iterable if only one subplot
            if n_page_samples == 1:
                axes = [axes]
            elif n_rows == 1:
                axes = axes.reshape(1, -1)
            elif n_cols == 1:
                axes = axes.reshape(-1, 1)
            
            # Create histograms for each sample on this page
            for i, sample_name in enumerate(page_samples):
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
            for i in range(n_page_samples, n_rows * n_cols):
                row = i // n_cols
                col = i % n_cols
                axes[row, col].set_visible(False)
            
            # Set common labels
            fig.text(0.5, 0.02, 'Score', ha='center', fontsize=12)
            fig.text(0.02, 0.5, 'Frequency', va='center', rotation='vertical', fontsize=12)
            
            # Add overall title with page info
            if num_pages > 1:
                title = f'Score Distribution Comparison (Grid Layout) - Page {page + 1} of {num_pages}'
            else:
                title = 'Score Distribution Comparison (Grid Layout)'
            fig.suptitle(title, fontsize=14, fontweight='bold')
            
            # Adjust layout
            plt.tight_layout()
            
            # Save page to PDF
            pdf.savefig(fig, dpi=300, bbox_inches='tight')
            print(f"Grid histogram page {page + 1} added to: {output_path}")
            
            plt.close()

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
    parser.add_argument('--sample_dirs', nargs='+', required=False,
                       help='Directories containing sample data (each should be a sample directory or its score_only subdir)')
    parser.add_argument('--base_dir', default='../outputs/my_variants/5UOI/sample_variant_scores',
                       help='Base directory to auto-discover sample score directories (default: ../outputs/my_variants/5UOI/sample_variant_scores)')
    parser.add_argument('--output_dir', default='../outputs/my_variants/5UOI/column_histogram_analysis',
                       help='Output directory for plots (default: ../outputs/my_variants/5UOI/column_histogram_analysis)')
    parser.add_argument('--bins', type=int, default=30,
                       help='Number of histogram bins (default: 30)')
    parser.add_argument('--use_mean_scores', action='store_true', default=True,
                       help='Use mean scores instead of individual scores (default: True)')
    parser.add_argument('--layout', choices=['column', 'grid', 'both'], default='both',
                       help='Layout type for histograms')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Helper to resolve a provided directory to its score_only directory and sample name
    def resolve_score_only_dir(path):
        score_only_dir = path
        if os.path.basename(score_only_dir) != 'score_only':
            candidate = os.path.join(score_only_dir, 'score_only')
            if os.path.isdir(candidate):
                score_only_dir = candidate
            else:
                return None, None
        sample_number = os.path.basename(os.path.dirname(score_only_dir))
        if not sample_number.isdigit():
            return None, None
        sample_name = f"Sample_{sample_number}"
        return score_only_dir, sample_name
    
    # Load data from all samples (either provided or auto-discovered)
    sample_data = {}
    
    if args.sample_dirs:
        for input_dir in args.sample_dirs:
            score_only_dir, sample_name = resolve_score_only_dir(input_dir)
            if score_only_dir is None or sample_name is None:
                print(f"Skipping {input_dir}: could not find a valid score_only directory with numeric sample name")
                continue
            mean_scores = load_mean_scores_from_npz_files(score_only_dir)
            sample_data[sample_name] = {'mean_scores': mean_scores}
            print(f"Loaded {len(mean_scores)} mean scores from {sample_name}")
    else:
        base_dir = args.base_dir
        if os.path.exists(base_dir):
            for item in sorted(os.listdir(base_dir)):
                item_path = os.path.join(base_dir, item)
                score_only_path = os.path.join(item_path, 'score_only')
                if os.path.isdir(item_path) and os.path.isdir(score_only_path) and item.isdigit():
                    mean_scores = load_mean_scores_from_npz_files(score_only_path)
                    sample_name = f"Sample_{item}"
                    sample_data[sample_name] = {'mean_scores': mean_scores}
                    print(f"Found sample directory: {item}; loaded {len(mean_scores)} mean scores")
        else:
            print(f"Warning: Base directory {base_dir} does not exist")
    
    # Ensure consistent ordering by sample number
    if sample_data:
        sample_data = {k: sample_data[k] for k in sorted(sample_data.keys(), key=lambda x: int(x.split('_')[1]))}
    
    # Prepare for pagination if more than 10 samples
    samples_per_page = 10
    total_samples = len(sample_data)
    if total_samples > samples_per_page:
        print(f"Found {total_samples} samples, will create multiple pages with max {samples_per_page} histograms per page")
    
    use_individual = not args.use_mean_scores
    
    # Generate plots based on layout choice
    if args.layout in ['column', 'both']:
        output_path = os.path.join(args.output_dir, "column_histograms.pdf")
        create_column_histograms(sample_data, output_path, bins=args.bins, 
                               use_individual_scores=use_individual, samples_per_page=samples_per_page)
    
    if args.layout in ['grid', 'both']:
        output_path = os.path.join(args.output_dir, "grid_histograms.pdf")
        create_shared_bin_histograms(sample_data, output_path, bins=args.bins, 
                                   use_individual_scores=use_individual, samples_per_page=samples_per_page)
    
    # Print statistical summary
    stats_summary = statistical_analysis(sample_data, use_individual_scores=use_individual)
    print("\nStatistical Summary:")
    print("=" * 50)
    for sample_name, stats in stats_summary.items():
        print(f"\n{sample_name}:")
        for stat_name, value in stats.items():
            print(f"  {stat_name}: {value:.4f}")

if __name__ == "__main__":
    main()
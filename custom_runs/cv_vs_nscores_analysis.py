#!/usr/bin/env python3
"""
Coefficient of Variation vs. Number of Scores Analysis for ProteinMPNN Variants

This script creates a plot showing CV vs. number of scores to evaluate the consistency
of each calculation, as required by the assignment.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects
import os
import glob
import argparse
import re

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

def create_cv_vs_nscores_plot(sample_data_dict, output_path, title="Coefficient of Variation vs. Number of Scores"):
    """
    Create a violin plot showing coefficient of variation vs. number of scores for all samples.
    This addresses the assignment requirement to evaluate CV vs. number of scores
    to determine the consistency of each calculation.
    
    Args:
        sample_data_dict (dict): Dictionary with sample names as keys and variant data as values
        output_path (str): Path to save the plot
        title (str): Plot title
    """
    plt.figure(figsize=(16, 12))
    
    # Group data by number of scores for violin plot
    score_groups = {}
    for sample_name, variant_data in sample_data_dict.items():
        if not variant_data:
            continue
        
        # Extract score count from sample name (e.g., "Sample_1_3scores" -> 3)
        if sample_name.endswith('scores'):
            score_count = int(sample_name.split('_')[-1].replace('scores', ''))
            if score_count not in score_groups:
                score_groups[score_count] = []
            score_groups[score_count].extend([data[0] for data in variant_data])  # Extract CV values
    
    if not score_groups:
        print("Warning: No valid score groups found for violin plot")
        return
    
    # Prepare data for violin plot
    score_counts = sorted(score_groups.keys())
    cv_data = [score_groups[count] for count in score_counts]
    
    # Create categorical positions for even spacing
    positions = list(range(len(score_counts)))
    
    # Create violin plot
    violin_parts = plt.violinplot(cv_data, positions=positions, showmeans=True, showmedians=True)
    
    # Customize violin plot appearance
    violin_parts['bodies'][0].set_facecolor('lightblue')
    violin_parts['bodies'][0].set_alpha(0.7)
    
    # Make mean and median lines more prominent and ensure they're above data points
    violin_parts['cmeans'].set_color('white')
    violin_parts['cmeans'].set_linewidth(2)
    violin_parts['cmeans'].set_zorder(20)
    violin_parts['cmeans'].set_path_effects([plt.matplotlib.patheffects.withStroke(linewidth=2, foreground='white')])
    
    violin_parts['cmedians'].set_color('black')
    violin_parts['cmedians'].set_linewidth(2)
    violin_parts['cmedians'].set_zorder(20)
    violin_parts['cmedians'].set_path_effects([plt.matplotlib.patheffects.withStroke(linewidth=2, foreground='white')])
    # Color each violin differently
    # Custom darker color palette (works well with white mean line)
    custom_colors = [
        "blue",  # blue
        "darkred",  # dark red
        "darkgreen",  # dark green
        "darkviolet",  # dark violet
    ]

    for i, (body, color) in enumerate(zip(violin_parts['bodies'], custom_colors)):
        body.set_facecolor(color)
        body.set_alpha(0.7)
        body.set_edgecolor("black")
        body.set_linewidth(0.5)

    

    # Add 10% of data points (chosen randomly) to each violin
    np.random.seed(42)  # For reproducible random selection
    for i, (score_count, cvs) in enumerate(zip(score_counts, cv_data)):
        # Randomly select 10% of the data points
        n_points = len(cvs)
        n_selected = max(1, n_points // 10)  # At least 1 point, 10% of total
        selected_indices = np.random.choice(n_points, size=n_selected, replace=False)
        selected_cvs = [cvs[j] for j in selected_indices]
        
        # Add jitter to x-position to avoid overlap
        jitter = np.random.normal(0, 0.1, len(selected_cvs))
        plt.scatter(positions[i] + jitter, selected_cvs, alpha=0.2, s=25, color='black', zorder=5, edgecolors='none')
    
    # Set x-axis ticks to show actual score count values
    plt.xticks(positions, score_counts)
    
    plt.xlabel('Number of Scores per Variant', fontsize=14)
    plt.ylabel('Coefficient of Variation (CV)', fontsize=14)
    plt.title(title, fontsize=16, fontweight='bold')
    
    # Create custom legend
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    import matplotlib.patheffects as pe

    legend_elements = [
        Patch(facecolor='darkgray', alpha=0.7, label='CV Distribution'),
        Line2D([0], [0], color='white', linewidth=2.5,
            path_effects=[pe.withStroke(linewidth=3.5, foreground='black')],
            label='Mean CV'),
        Line2D([0], [0], color='black', linewidth=2,
            path_effects=[pe.withStroke(linewidth=3)],
            label='Median CV')
]

    plt.legend(handles=legend_elements, loc='upper right', fontsize=12)

    
    plt.grid(True, alpha=0.3)
    
    # Use log scale for better visualization if CV values vary widely
    all_cvs = []
    for cvs in cv_data:
        all_cvs.extend(cvs)
    
    if any(cv > 10 for cv in all_cvs):
        plt.yscale('log')
        plt.ylabel('Coefficient of Variation (CV) - Log Scale', fontsize=14)
    
    # Add some statistical insights
    # Calculate overall trends
    if len(all_cvs) > 1:
        # Group by number of scores and calculate statistics
        score_stats = {}
        for score_count, cvs in zip(score_counts, cv_data):
            if cvs:
                score_stats[score_count] = {
                    'mean': np.mean(cvs),
                    'median': np.median(cvs),
                    'std': np.std(cvs),
                    'count': len(cvs)
                }
        
        # Add text box with insights
        insight_text = "Statistical Summary:\n"
        for score_count in sorted(score_stats.keys()):
            stats = score_stats[score_count]
            insight_text += f"{score_count} scores:\n"
            insight_text += f"  Mean CV: {stats['mean']:.4f}\n"
            insight_text += f"  Median CV: {stats['median']:.4f}\n"
            insight_text += f"  Std CV: {stats['std']:.4f}\n"
            insight_text += f"  N variants: {stats['count']}\n\n"
        
        plt.text(1.05, 0.98, insight_text, transform=plt.gca().transAxes, 
                fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Adjust layout
    plt.tight_layout()
    
    # Save plot
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"CV vs. Number of Scores violin plot saved to: {output_path}")

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
        default='Coefficient of Variation vs. Number of Scores',
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
        # Default: Automatically discover all score directories and sample directories
        base_dir = "../outputs/my_variants/5UOI/sample_variant_scores"
        sample_data = {}
        
        if os.path.exists(base_dir):
            # Look for any score directories (e.g., 3score, 5score, 10score, 15score, etc.)
            score_dirs = []
            for item in sorted(os.listdir(base_dir)):
                item_path = os.path.join(base_dir, item)
                if os.path.isdir(item_path) and 'score' in item.lower():
                    score_dirs.append(item_path)
                    print(f"Found score directory: {item}")
            
            print(f"Discovered {len(score_dirs)} score directories")
            
            # Process each score directory
            for score_dir in score_dirs:
                # Extract score count from directory name (more flexible parsing)
                dir_name = os.path.basename(score_dir)
                
                # Try to extract the number before "score"
                score_match = re.search(r'(\d+)', dir_name)
                if score_match:
                    score_count = int(score_match.group(1))
                else:
                    print(f"Warning: Could not extract score count from directory name: {dir_name}")
                    continue
                
                print(f"\nProcessing {score_count} scores directory...")
                
                # Get all sample subdirectories
                sample_dirs = []
                for item in sorted(os.listdir(score_dir)):
                    item_path = os.path.join(score_dir, item)
                    score_only_path = os.path.join(item_path, "score_only")
                    
                    # Check if this is a directory and contains a score_only subdirectory
                    if os.path.isdir(item_path) and os.path.isdir(score_only_path):
                        # Check if it's a numeric sample directory
                        if item.isdigit():
                            sample_dirs.append(score_only_path)
                        else:
                            print(f"  Skipping non-numeric directory: {item}")
                
                print(f"  Found {len(sample_dirs)} sample directories with {score_count} scores")
                
                # Process all discovered sample directories for this score count
                for sample_dir in sample_dirs:
                    if os.path.exists(sample_dir):
                        # Extract sample number from the path
                        sample_number = os.path.basename(os.path.dirname(sample_dir))
                        sample_name = f"Sample_{sample_number}_{score_count}scores"
                        
                        # Load CV data for this sample
                        variant_data = load_variant_cv_data(sample_dir)
                        
                        if variant_data:
                            # All variants should have the same number of scores
                            for cv, n_scores in variant_data:
                                if n_scores != score_count:
                                    print(f"  Warning: Expected {score_count} scores but got {n_scores} for sample {sample_number}")
                            
                            sample_data[sample_name] = variant_data
                            print(f"    Sample {sample_number}: {len(variant_data)} variants")
                        else:
                            print(f"    Sample {sample_number}: No valid variants found")
                    else:
                        print(f"  Warning: Sample directory does not exist: {sample_dir}")
        else:
            print(f"Warning: Base directory {base_dir} does not exist")
        
        # Sort sample_data by sample number and score count
        sorted_sample_data = {}
        for sample_name in sorted(sample_data.keys(), key=lambda x: (int(x.split('_')[1]), int(x.split('_')[2].replace('scores', '')))):
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

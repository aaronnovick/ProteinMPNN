#!/usr/bin/env python3
"""
Simple script to calculate z-scores for ProteinMPNN scores from .npz files.

This script processes .npz files and calculates z-scores for normalization.
Z-scores are calculated as: (score - mean) / std

Usage examples:
    # Process all .npz files in a directory
    python simple_score_zscore.py --input_dir "outputs/my_variants/5UOI/sample_variant_scores/15score/9/score_only/"
    
    # Process specific files
    python simple_score_zscore.py --files "file1.npz" "file2.npz" "file3.npz"
"""

import numpy as np
import os
import glob
import argparse
from typing import List, Union, Tuple, Dict
import pandas as pd


def load_scores_from_files(input_paths: List[str], labels: List[str] = None) -> tuple:
    """
    Load scores from .npz files.
    
    Args:
        input_paths: List of file paths or directory paths (can be multiple libraries)
        labels: Optional list of labels corresponding to each provided input path. If not provided,
                directory basenames or 'file' will be used.
        
    Returns:
        Tuple of (scores, file_names, sequences, groups)
    """
    all_scores = []
    file_names = []
    sequences = []
    groups = []
    
    # Expand directories to individual files
    file_paths: List[str] = []
    file_to_group: Dict[str, str] = {}

    labels = labels or []
    if labels and len(labels) != len(input_paths):
        raise ValueError("If providing --labels, the number of labels must match the number of input paths")

    for idx, path in enumerate(input_paths):
        label = None
        if labels and idx < len(labels):
            label = labels[idx]
        else:
            label = os.path.basename(os.path.normpath(path)) if os.path.isdir(path) else "file"

        if os.path.isdir(path):
            expanded = glob.glob(os.path.join(path, "**", "*.npz"), recursive=True)
            for p in expanded:
                file_paths.append(p)
                file_to_group[p] = label
        else:
            file_paths.append(path)
            file_to_group[path] = label
    
    print(f"Processing {len(file_paths)} files...")
    
    for file_path in file_paths:
        try:
            data = np.load(file_path)
            
            # Extract score (prefer 'score' over 'global_score')
            if 'score' in data:
                score_data = data['score']
            elif 'global_score' in data:
                score_data = data['global_score']
            else:
                print(f"Warning: No score data in {file_path}")
                continue
            
            # Handle both scalar and array scores
            if score_data.ndim == 0:  # Scalar
                score = score_data.item()
            else:  # Array - take mean
                score = np.mean(score_data)
            
            all_scores.append(score)
            file_names.append(os.path.basename(file_path))
            groups.append(file_to_group.get(file_path, ""))
            
            # Extract sequence if available
            if 'seq_str' in data:
                seq_data = data['seq_str']
                if isinstance(seq_data, np.ndarray) and seq_data.ndim == 0:
                    sequences.append(seq_data.item())
                else:
                    sequences.append(str(seq_data))
            else:
                sequences.append("")
                
        except Exception as e:
            print(f"Error loading {file_path}: {e}")
            continue
    
    return np.array(all_scores), file_names, sequences, groups


def calculate_z_scores(scores: np.ndarray) -> np.ndarray:
    """Calculate z-scores for the given scores."""
    mean_score = np.mean(scores)
    std_score = np.std(scores, ddof=1)  # Sample standard deviation
    
    if std_score == 0:
        print("Warning: Standard deviation is 0, all z-scores will be 0")
        return np.zeros_like(scores)
    
    z_scores = (scores - mean_score) / std_score
    return z_scores


def print_summary(scores: np.ndarray, z_scores: np.ndarray, file_names: List[str], groups: List[str] = None):
    """Print summary statistics."""
    print("\n" + "="*60)
    print("Z-SCORE NORMALIZATION SUMMARY")
    print("="*60)
    print(f"Total sequences: {len(scores)}")
    print(f"\nOriginal scores:")
    print(f"  Mean: {np.mean(scores):.4f}")
    print(f"  Std:  {np.std(scores, ddof=1):.4f}")
    print(f"  Min:  {np.min(scores):.4f}")
    print(f"  Max:  {np.max(scores):.4f}")
    print(f"\nZ-scores:")
    print(f"  Mean: {np.mean(z_scores):.4f}")
    print(f"  Std:  {np.std(z_scores, ddof=1):.4f}")
    print(f"  Min:  {np.min(z_scores):.4f}")
    print(f"  Max:  {np.max(z_scores):.4f}")
    print("="*60)
    
    # Group-wise summary if groups provided
    if groups:
        try:
            import pandas as _pd
            df = _pd.DataFrame({
                'group': groups,
                'z': z_scores,
                'score': scores,
            })
            group_stats = df.groupby('group').agg(
                count=('z', 'size'),
                mean_z=('z', 'mean'),
                std_z=('z', 'std'),
                mean_score=('score', 'mean'),
                std_score=('score', 'std'),
            ).reset_index()
            print("\nPer-directory (library) z-score means:")
            for _, row in group_stats.iterrows():
                print(f"  {row['group']}: n={int(row['count'])}, mean_z={row['mean_z']:.4f}, std_z={0.0 if np.isnan(row['std_z']) else row['std_z']:.4f}")
            # Average of per-group mean z-scores (treating each library equally)
            avg_of_group_means = group_stats['mean_z'].mean()
            print(f"\nAverage of per-library mean z-scores: {avg_of_group_means:.4f}")
        except Exception as _:
            pass
    
    # Show top 5 and bottom 5 sequences
    sorted_indices = np.argsort(z_scores)
    print(f"\nTop 5 sequences (highest z-scores):")
    for i in range(min(5, len(scores))):
        idx = sorted_indices[-(i+1)]
        suffix = f" [{groups[idx]}]" if groups else ""
        print(f"  {i+1}. {file_names[idx]}{suffix}: z-score = {z_scores[idx]:.3f}, score = {scores[idx]:.3f}")
    
    print(f"\nBottom 5 sequences (lowest z-scores):")
    for i in range(min(5, len(scores))):
        idx = sorted_indices[i]
        suffix = f" [{groups[idx]}]" if groups else ""
        print(f"  {i+1}. {file_names[idx]}{suffix}: z-score = {z_scores[idx]:.3f}, score = {scores[idx]:.3f}")


def save_results(scores: np.ndarray, z_scores: np.ndarray, file_names: List[str], 
                sequences: List[str], output_file: str = None, groups: List[str] = None):
    """Save results to CSV file."""
    if output_file is None:
        output_file = "score_zscore_results.csv"
    
    df_dict = {
        'filename': file_names,
        'sequence': sequences,
        'original_score': scores,
        'z_score': z_scores
    }
    if groups:
        df_dict['group'] = groups
    df = pd.DataFrame(df_dict)
    
    # Sort by z-score (highest first)
    df = df.sort_values('z_score', ascending=False)
    
    df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")

    # If groups present, also write a small group summary next to the CSV
    if groups:
        summary_path = os.path.splitext(output_file)[0] + "_group_summary.csv"
        group_stats = df.groupby('group').agg(
            count=('z_score', 'size'),
            mean_z=('z_score', 'mean'),
            std_z=('z_score', 'std'),
            mean_score=('original_score', 'mean'),
            std_score=('original_score', 'std'),
        ).reset_index()
        group_stats['avg_of_group_means'] = group_stats['mean_z'].mean()
        group_stats.to_csv(summary_path, index=False)
        print(f"Group summary saved to: {summary_path}")


def main():
    parser = argparse.ArgumentParser(description="Calculate z-scores for ProteinMPNN scores")
    parser.add_argument("--input_dir", type=str, 
                       help="Directory containing .npz files")
    parser.add_argument("--dirs", nargs="+",
                       help="Multiple directories (e.g., 5 libraries for the same paratope)")
    parser.add_argument("--labels", nargs="+",
                       help="Optional labels for each provided --dirs path (must match length)")
    parser.add_argument("--files", nargs="+", 
                       help="Specific .npz files to process")
    parser.add_argument("--output", type=str, default="score_zscore_results.csv",
                       help="Output CSV file (default: score_zscore_results.csv)")
    
    args = parser.parse_args()
    
    if not args.input_dir and not args.files and not args.dirs:
        print("Error: Must specify one of --dirs, --input_dir, or --files")
        return 1
    
    # Determine input paths
    input_paths = []
    labels = None
    if args.dirs:
        input_paths.extend(args.dirs)
        labels = args.labels
    if args.input_dir:
        input_paths.append(args.input_dir)
    if args.files:
        input_paths.extend(args.files)
    
    try:
        # Load scores
        scores, file_names, sequences, groups = load_scores_from_files(input_paths, labels=labels)
        
        if len(scores) == 0:
            print("No valid score data found!")
            return 1
        
        # Calculate z-scores
        z_scores = calculate_z_scores(scores)
        
        # Print summary
        print_summary(scores, z_scores, file_names, groups)
        
        # Save results
        save_results(scores, z_scores, file_names, sequences, args.output, groups)
        
        print("\nZ-score normalization completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())

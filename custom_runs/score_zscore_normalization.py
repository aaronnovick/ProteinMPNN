#!/usr/bin/env python3
"""
Script to normalize ProteinMPNN negative log likelihood scores by calculating z-scores.

This script processes .npz files containing ProteinMPNN scores and calculates z-scores
for normalization. Z-scores are calculated as: (score - mean) / std

Usage:
    python score_zscore_normalization.py --input_dir <path_to_npz_files> --output_dir <output_path>
"""

import numpy as np
import os
import glob
import argparse
from typing import List, Dict, Tuple
import pandas as pd


def load_scores_from_npz_files(input_dir: str) -> Tuple[np.ndarray, List[str], List[str]]:
    """
    Load all scores from .npz files in a directory.
    
    Args:
        input_dir: Directory containing .npz files
        
    Returns:
        Tuple of (scores_array, file_paths, sequence_strings)
    """
    npz_files = glob.glob(os.path.join(input_dir, "**", "*.npz"), recursive=True)
    
    if not npz_files:
        raise ValueError(f"No .npz files found in {input_dir}")
    
    print(f"Found {len(npz_files)} .npz files")
    
    all_scores = []
    file_paths = []
    sequence_strings = []
    
    for npz_file in npz_files:
        try:
            data = np.load(npz_file)
            
            # Extract scores (use 'score' if available, otherwise 'global_score')
            if 'score' in data:
                scores = data['score']
            elif 'global_score' in data:
                scores = data['global_score']
            else:
                print(f"Warning: No score data found in {npz_file}")
                continue
            
            # Handle both scalar and array scores
            if scores.ndim == 0:  # Scalar
                all_scores.append(scores.item())
            else:  # Array - take mean for global score
                all_scores.append(np.mean(scores))
            
            file_paths.append(npz_file)
            
            # Extract sequence string if available
            if 'seq_str' in data:
                seq_str = data['seq_str']
                if isinstance(seq_str, np.ndarray) and seq_str.ndim == 0:
                    sequence_strings.append(seq_str.item())
                else:
                    sequence_strings.append(str(seq_str))
            else:
                sequence_strings.append("")
                
        except Exception as e:
            print(f"Error loading {npz_file}: {e}")
            continue
    
    return np.array(all_scores), file_paths, sequence_strings


def calculate_z_scores(scores: np.ndarray) -> np.ndarray:
    """
    Calculate z-scores for the given scores.
    
    Args:
        scores: Array of scores
        
    Returns:
        Array of z-scores
    """
    mean_score = np.mean(scores)
    std_score = np.std(scores, ddof=1)  # Use sample standard deviation
    
    if std_score == 0:
        print("Warning: Standard deviation is 0, all z-scores will be 0")
        return np.zeros_like(scores)
    
    z_scores = (scores - mean_score) / std_score
    return z_scores


def save_results(scores: np.ndarray, z_scores: np.ndarray, file_paths: List[str], 
                sequence_strings: List[str], output_dir: str) -> None:
    """
    Save the results to CSV and summary files.
    
    Args:
        scores: Original scores
        z_scores: Calculated z-scores
        file_paths: List of file paths
        sequence_strings: List of sequence strings
        output_dir: Output directory
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # Create DataFrame
    df = pd.DataFrame({
        'file_path': file_paths,
        'sequence': sequence_strings,
        'original_score': scores,
        'z_score': z_scores
    })
    
    # Add filename column for easier identification
    df['filename'] = df['file_path'].apply(lambda x: os.path.basename(x))
    
    # Save detailed results
    output_csv = os.path.join(output_dir, "score_zscore_results.csv")
    df.to_csv(output_csv, index=False)
    print(f"Detailed results saved to: {output_csv}")
    
    # Save summary statistics
    summary_stats = {
        'total_sequences': len(scores),
        'mean_original_score': np.mean(scores),
        'std_original_score': np.std(scores, ddof=1),
        'min_original_score': np.min(scores),
        'max_original_score': np.max(scores),
        'mean_z_score': np.mean(z_scores),
        'std_z_score': np.std(z_scores, ddof=1),
        'min_z_score': np.min(z_scores),
        'max_z_score': np.max(z_scores)
    }
    
    summary_df = pd.DataFrame([summary_stats])
    summary_csv = os.path.join(output_dir, "score_zscore_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"Summary statistics saved to: {summary_csv}")
    
    # Print summary to console
    print("\n" + "="*50)
    print("Z-SCORE NORMALIZATION SUMMARY")
    print("="*50)
    print(f"Total sequences processed: {summary_stats['total_sequences']}")
    print(f"\nOriginal scores:")
    print(f"  Mean: {summary_stats['mean_original_score']:.4f}")
    print(f"  Std:  {summary_stats['std_original_score']:.4f}")
    print(f"  Range: [{summary_stats['min_original_score']:.4f}, {summary_stats['max_original_score']:.4f}]")
    print(f"\nZ-scores:")
    print(f"  Mean: {summary_stats['mean_z_score']:.4f}")
    print(f"  Std:  {summary_stats['std_z_score']:.4f}")
    print(f"  Range: [{summary_stats['min_z_score']:.4f}, {summary_stats['max_z_score']:.4f}]")
    print("="*50)


def main():
    parser = argparse.ArgumentParser(description="Normalize ProteinMPNN scores using z-scores")
    parser.add_argument("--input_dir", required=True, 
                       help="Directory containing .npz files with scores")
    parser.add_argument("--output_dir", required=True,
                       help="Output directory for results")
    parser.add_argument("--pattern", default="**/*.npz",
                       help="File pattern to match (default: **/*.npz)")
    
    args = parser.parse_args()
    
    try:
        # Load scores from all .npz files
        print(f"Loading scores from: {args.input_dir}")
        scores, file_paths, sequence_strings = load_scores_from_npz_files(args.input_dir)
        
        if len(scores) == 0:
            print("No valid score data found!")
            return
        
        # Calculate z-scores
        print("Calculating z-scores...")
        z_scores = calculate_z_scores(scores)
        
        # Save results
        print(f"Saving results to: {args.output_dir}")
        save_results(scores, z_scores, file_paths, sequence_strings, args.output_dir)
        
        print("\nZ-score normalization completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())

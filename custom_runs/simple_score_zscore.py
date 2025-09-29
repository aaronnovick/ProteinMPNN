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
from typing import List, Union
import pandas as pd


def load_scores_from_files(input_paths: List[str]) -> tuple:
    """
    Load scores from .npz files.
    
    Args:
        input_paths: List of file paths or directory paths
        
    Returns:
        Tuple of (scores, file_names, sequences)
    """
    all_scores = []
    file_names = []
    sequences = []
    
    # Expand directories to individual files
    file_paths = []
    for path in input_paths:
        if os.path.isdir(path):
            file_paths.extend(glob.glob(os.path.join(path, "**", "*.npz"), recursive=True))
        else:
            file_paths.append(path)
    
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
    
    return np.array(all_scores), file_names, sequences


def calculate_z_scores(scores: np.ndarray) -> np.ndarray:
    """Calculate z-scores for the given scores."""
    mean_score = np.mean(scores)
    std_score = np.std(scores, ddof=1)  # Sample standard deviation
    
    if std_score == 0:
        print("Warning: Standard deviation is 0, all z-scores will be 0")
        return np.zeros_like(scores)
    
    z_scores = (scores - mean_score) / std_score
    return z_scores


def print_summary(scores: np.ndarray, z_scores: np.ndarray, file_names: List[str]):
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
    
    # Show top 5 and bottom 5 sequences
    sorted_indices = np.argsort(z_scores)
    print(f"\nTop 5 sequences (highest z-scores):")
    for i in range(min(5, len(scores))):
        idx = sorted_indices[-(i+1)]
        print(f"  {i+1}. {file_names[idx]}: z-score = {z_scores[idx]:.3f}, score = {scores[idx]:.3f}")
    
    print(f"\nBottom 5 sequences (lowest z-scores):")
    for i in range(min(5, len(scores))):
        idx = sorted_indices[i]
        print(f"  {i+1}. {file_names[idx]}: z-score = {z_scores[idx]:.3f}, score = {scores[idx]:.3f}")


def save_results(scores: np.ndarray, z_scores: np.ndarray, file_names: List[str], 
                sequences: List[str], output_file: str = None):
    """Save results to CSV file."""
    if output_file is None:
        output_file = "score_zscore_results.csv"
    
    df = pd.DataFrame({
        'filename': file_names,
        'sequence': sequences,
        'original_score': scores,
        'z_score': z_scores
    })
    
    # Sort by z-score (highest first)
    df = df.sort_values('z_score', ascending=False)
    
    df.to_csv(output_file, index=False)
    print(f"\nResults saved to: {output_file}")


def main():
    parser = argparse.ArgumentParser(description="Calculate z-scores for ProteinMPNN scores")
    parser.add_argument("--input_dir", type=str, 
                       help="Directory containing .npz files")
    parser.add_argument("--files", nargs="+", 
                       help="Specific .npz files to process")
    parser.add_argument("--output", type=str, default="score_zscore_results.csv",
                       help="Output CSV file (default: score_zscore_results.csv)")
    
    args = parser.parse_args()
    
    if not args.input_dir and not args.files:
        print("Error: Must specify either --input_dir or --files")
        return 1
    
    # Determine input paths
    input_paths = []
    if args.input_dir:
        input_paths.append(args.input_dir)
    if args.files:
        input_paths.extend(args.files)
    
    try:
        # Load scores
        scores, file_names, sequences = load_scores_from_files(input_paths)
        
        if len(scores) == 0:
            print("No valid score data found!")
            return 1
        
        # Calculate z-scores
        z_scores = calculate_z_scores(scores)
        
        # Print summary
        print_summary(scores, z_scores, file_names)
        
        # Save results
        save_results(scores, z_scores, file_names, sequences, args.output)
        
        print("\nZ-score normalization completed successfully!")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())

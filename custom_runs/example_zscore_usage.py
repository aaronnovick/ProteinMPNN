#!/usr/bin/env python3
"""
Example script showing how to use the z-score normalization for ProteinMPNN scores.

This script demonstrates different ways to calculate z-scores for your designed variants.
"""

import os
import sys

# Add the custom_runs directory to the path so we can import our modules
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from simple_score_zscore import load_scores_from_files, calculate_z_scores, print_summary, save_results


def example_1_specific_directory():
    """Example 1: Process all .npz files in a specific directory"""
    print("="*60)
    print("EXAMPLE 1: Processing all .npz files in a directory")
    print("="*60)
    
    # Path to your variant scores directory
    input_dir = "../outputs/my_variants/5UOI/sample_variant_scores/15score/9/score_only/"
    
    if not os.path.exists(input_dir):
        print(f"Directory not found: {input_dir}")
        print("Please update the path to point to your actual .npz files directory")
        return
    
    # Load scores
    scores, file_names, sequences = load_scores_from_files([input_dir])
    
    if len(scores) == 0:
        print("No scores found in the directory")
        return
    
    # Calculate z-scores
    z_scores = calculate_z_scores(scores)
    
    # Print summary
    print_summary(scores, z_scores, file_names)
    
    # Save results
    save_results(scores, z_scores, file_names, sequences, "example1_zscore_results.csv")


def example_2_specific_files():
    """Example 2: Process specific .npz files"""
    print("\n" + "="*60)
    print("EXAMPLE 2: Processing specific .npz files")
    print("="*60)
    
    # List specific files you want to analyze
    specific_files = [
        "../outputs/my_variants/5UOI/sample_variant_scores/15score/9/score_only/5UOI_fasta_100.npz",
        "../outputs/my_variants/5UOI/sample_variant_scores/15score/9/score_only/5UOI_fasta_99.npz",
        "../outputs/my_variants/5UOI/sample_variant_scores/15score/9/score_only/5UOI_fasta_98.npz",
        # Add more files as needed
    ]
    
    # Filter to only existing files
    existing_files = [f for f in specific_files if os.path.exists(f)]
    
    if not existing_files:
        print("None of the specified files exist. Please update the file paths.")
        return
    
    # Load scores
    scores, file_names, sequences = load_scores_from_files(existing_files)
    
    if len(scores) == 0:
        print("No valid scores found in the files")
        return
    
    # Calculate z-scores
    z_scores = calculate_z_scores(scores)
    
    # Print summary
    print_summary(scores, z_scores, file_names)
    
    # Save results
    save_results(scores, z_scores, file_names, sequences, "example2_zscore_results.csv")


def example_3_compare_different_samples():
    """Example 3: Compare z-scores across different samples"""
    print("\n" + "="*60)
    print("EXAMPLE 3: Comparing z-scores across different samples")
    print("="*60)
    
    # Define different sample directories
    sample_dirs = [
        "../outputs/my_variants/5UOI/sample_variant_scores/15score/9/score_only/",
        # Add more sample directories as needed
    ]
    
    all_results = {}
    
    for i, sample_dir in enumerate(sample_dirs):
        if not os.path.exists(sample_dir):
            print(f"Sample directory not found: {sample_dir}")
            continue
        
        print(f"\nProcessing sample {i+1}: {sample_dir}")
        
        # Load scores for this sample
        scores, file_names, sequences = load_scores_from_files([sample_dir])
        
        if len(scores) == 0:
            print(f"No scores found in {sample_dir}")
            continue
        
        # Calculate z-scores
        z_scores = calculate_z_scores(scores)
        
        # Store results
        all_results[f"sample_{i+1}"] = {
            'scores': scores,
            'z_scores': z_scores,
            'file_names': file_names,
            'sequences': sequences
        }
        
        print(f"Sample {i+1} summary:")
        print(f"  Mean score: {scores.mean():.4f} ± {scores.std():.4f}")
        print(f"  Mean z-score: {z_scores.mean():.4f} ± {z_scores.std():.4f}")
    
    # Compare samples
    if len(all_results) > 1:
        print(f"\nComparison across {len(all_results)} samples:")
        for sample_name, data in all_results.items():
            print(f"  {sample_name}: {len(data['scores'])} sequences, "
                  f"mean z-score = {data['z_scores'].mean():.3f}")


def main():
    """Run all examples"""
    print("ProteinMPNN Score Z-Score Normalization Examples")
    print("=" * 60)
    
    # Run examples
    example_1_specific_directory()
    example_2_specific_files()
    example_3_compare_different_samples()
    
    print("\n" + "="*60)
    print("All examples completed!")
    print("Check the generated CSV files for detailed results.")
    print("="*60)


if __name__ == "__main__":
    main()

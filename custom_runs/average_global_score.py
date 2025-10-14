#!/usr/bin/env python3
"""
Average Global Score Calculator

This script calculates the average of the mean global_score from all .npz files 
within each library subdirectory. It searches through subdirectories to find 
individual libraries (e.g., 169, 210, 231) and computes the average for each library.

Note: Files ending with "5UOI_pdb.npz" are automatically excluded from the calculation.

Usage:
    python3 average_global_score.py <directory_path> [--top-n N]
    
Examples:
    python3 average_global_score.py ../McConnell_variants/aaa/outputs/helix2/scores
    python3 average_global_score.py ../McConnell_variants/aaa/outputs/helix2/scores --top-n 10
    python3 average_global_score.py ../McConnell_variants/aaa/outputs/helix2/scores --top-n 25
"""

import numpy as np
import argparse
import os
import sys
from pathlib import Path
import glob

def calculate_average_global_score_per_library(directory_path, top_n=None):
    """
    Calculate the average of mean global_score from all .npz files in each library subdirectory.
    Optionally, only average the top N lowest scores.
    
    Args:
        directory_path (str): Path to the directory containing library subdirectories
        top_n (int, optional): If provided, only average the top N lowest scores. If None, average all scores.
        
    Returns:
        dict: Dictionary with library names as keys and results as values
    """
    # Convert to Path object for easier handling
    dir_path = Path(directory_path)
    
    if not dir_path.exists():
        raise FileNotFoundError(f"Directory '{directory_path}' does not exist")
    
    if not dir_path.is_dir():
        raise ValueError(f"'{directory_path}' is not a directory")
    
    # Find all library subdirectories (directories that contain score_only subdirectories)
    library_dirs = []
    for item in dir_path.iterdir():
        if item.is_dir():
            score_only_dir = item / "score_only"
            if score_only_dir.exists() and score_only_dir.is_dir():
                library_dirs.append(item)
    
    if not library_dirs:
        raise ValueError(f"No library directories found in '{directory_path}'")
    
    print(f"Found {len(library_dirs)} library directories")
    
    library_results = {}
    
    # Process each library directory
    for library_dir in sorted(library_dirs):
        library_name = library_dir.name
        print(f"\nProcessing library: {library_name}")
        
        # Find all .npz files in the score_only subdirectory
        score_only_dir = library_dir / "score_only"
        all_npz_files = list(score_only_dir.glob("*.npz"))
        
        # Filter out 5UOI_pdb.npz files specifically
        npz_files = [f for f in all_npz_files if not f.name.endswith("5UOI_pdb.npz")]
        excluded_files = [f for f in all_npz_files if f.name.endswith("5UOI_pdb.npz")]
        
        if not npz_files:
            print(f"  Warning: No .npz files found in {score_only_dir} (after excluding 5UOI_pdb.npz files)")
            library_results[library_name] = {
                'average_score': None,
                'successful_files': 0,
                'individual_means': [],
                'failed_files': [],
                'excluded_files': excluded_files,
                'total_files': len(all_npz_files)
            }
            continue
        
        print(f"  Found {len(all_npz_files)} total .npz files")
        print(f"  Excluded {len(excluded_files)} 5UOI_pdb.npz files")
        print(f"  Processing {len(npz_files)} files")
        
        individual_means = []
        successful_files = 0
        failed_files = []
        
        # Process each .npz file in this library
        for npz_file in npz_files:
            try:
                # Load the .npz file
                data = np.load(npz_file)
                
                # Check if global_score exists
                if 'global_score' not in data:
                    print(f"  Warning: 'global_score' not found in {npz_file.name}")
                    failed_files.append(str(npz_file))
                    data.close()
                    continue
                
                # Get the global_score array
                global_score = data['global_score']
                
                # Calculate the mean of the global_score array
                mean_score = np.mean(global_score)
                individual_means.append(mean_score)
                successful_files += 1
                
                # Close the file
                data.close()
                
            except Exception as e:
                print(f"  Error processing {npz_file.name}: {e}")
                failed_files.append(str(npz_file))
                continue
        
        if not individual_means:
            print(f"  Warning: No valid .npz files with global_score found in {library_name}")
            average_score = None
        else:
            # Sort the scores to get the lowest values first
            sorted_means = sorted(individual_means)
            
            # Select top N scores if specified, otherwise use all scores
            if top_n is not None and top_n > 0:
                # Use top N lowest scores (first N in sorted list)
                selected_scores = sorted_means[:min(top_n, len(sorted_means))]
                average_score = np.mean(selected_scores)
                print(f"  Library {library_name} average: {average_score:.6f} (from top {len(selected_scores)} of {successful_files} files)")
            else:
                # Calculate the average of all individual means for this library
                average_score = np.mean(individual_means)
                print(f"  Library {library_name} average: {average_score:.6f} (from {successful_files} files)")
        
        library_results[library_name] = {
            'average_score': average_score,
            'successful_files': successful_files,
            'individual_means': individual_means,
            'failed_files': failed_files,
            'excluded_files': excluded_files,
            'total_files': len(all_npz_files),
            'top_n_used': top_n if top_n is not None else None
        }
    
    return library_results

def main():
    """Main function to handle command line arguments and execute the calculation."""
    parser = argparse.ArgumentParser(
        description="Calculate average of mean global_score from .npz files in a directory",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 average_global_score.py McConnell_variants/aaa/outputs/helix2/scores
  python3 average_global_score.py /path/to/npz/directory
  python3 average_global_score.py McConnell_variants/aaa/outputs/helix2/scores --top-n 10
  python3 average_global_score.py /path/to/npz/directory --top-n 25 --verbose
        """
    )
    
    parser.add_argument(
        'directory',
        help='Directory containing library subdirectories with score_only folders'
    )
    
    parser.add_argument(
        '--top-n',
        type=int,
        default=None,
        help='Only average the top N lowest scores (e.g., --top-n 10 for top 10 scores). If not specified, averages all scores.'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Show detailed output including individual file statistics'
    )
    
    parser.add_argument(
        '--output', '-o',
        help='Save results to a file (optional)'
    )
    
    args = parser.parse_args()
    
    try:
        # Calculate the average global score for each library
        library_results = calculate_average_global_score_per_library(args.directory, args.top_n)
        
        # Print results
        print(f"\n{'='*80}")
        print("LIBRARY RESULTS SUMMARY")
        print(f"{'='*80}")
        print(f"Directory: {args.directory}")
        print(f"Total libraries processed: {len(library_results)}")
        if args.top_n:
            print(f"Top N filtering: Using top {args.top_n} lowest scores per library")
        else:
            print("Top N filtering: Using all scores (no filtering)")
        
        # Calculate overall statistics
        valid_libraries = {name: data for name, data in library_results.items() if data['average_score'] is not None}
        if valid_libraries:
            library_averages = [data['average_score'] for data in valid_libraries.values()]
            overall_average = np.mean(library_averages)
            print(f"Overall average across libraries: {overall_average:.6f}")
            print(f"Standard deviation across libraries: {np.std(library_averages):.6f}")
            print(f"Min library average: {np.min(library_averages):.6f}")
            print(f"Max library average: {np.max(library_averages):.6f}")
        
        print(f"\n{'='*80}")
        print("LIBRARY DETAILS")
        print(f"{'='*80}")
        if args.top_n:
            print(f"{'Library':<10} {'Average Score':<15} {'Top N Used':<10} {'Total Files':<12} {'Excluded':<10} {'Failed':<8}")
            print(f"{'-'*10} {'-'*15} {'-'*10} {'-'*12} {'-'*10} {'-'*8}")
        else:
            print(f"{'Library':<10} {'Average Score':<15} {'Files':<8} {'Excluded':<10} {'Failed':<8}")
            print(f"{'-'*10} {'-'*15} {'-'*8} {'-'*10} {'-'*8}")
        
        for library_name, data in sorted(library_results.items()):
            avg_score = f"{data['average_score']:.6f}" if data['average_score'] is not None else "N/A"
            files_processed = data['successful_files']
            excluded_count = len(data['excluded_files'])
            failed_count = len(data['failed_files'])
            
            if args.top_n:
                top_n_used = min(args.top_n, files_processed) if files_processed > 0 else 0
                print(f"{library_name:<10} {avg_score:<15} {top_n_used:<10} {files_processed:<12} {excluded_count:<10} {failed_count:<8}")
            else:
                print(f"{library_name:<10} {avg_score:<15} {files_processed:<8} {excluded_count:<10} {failed_count:<8}")
        
        if args.verbose:
            print(f"\n{'='*80}")
            print("DETAILED LIBRARY STATISTICS")
            print(f"{'='*80}")
            
            for library_name, data in sorted(library_results.items()):
                if data['average_score'] is not None:
                    print(f"\nLibrary: {library_name}")
                    print(f"{'-'*40}")
                    print(f"Average score: {data['average_score']:.6f}")
                    print(f"Files processed: {data['successful_files']}")
                    if args.top_n and data['successful_files'] > 0:
                        top_n_used = min(args.top_n, data['successful_files'])
                        print(f"Top N used: {top_n_used} (out of {data['successful_files']} total)")
                    print(f"Files excluded: {len(data['excluded_files'])}")
                    print(f"Files failed: {len(data['failed_files'])}")
                    
                    if data['individual_means']:
                        print(f"Individual file statistics (sorted by score):")
                        sorted_means = sorted(data['individual_means'])
                        for i, mean_score in enumerate(sorted_means):
                            marker = " *" if args.top_n and i < min(args.top_n, len(sorted_means)) else ""
                            print(f"  File {i+1}: {mean_score:.6f}{marker}")
        
        # Save to file if requested
        if args.output:
            with open(args.output, 'w') as f:
                f.write(f"Directory: {args.directory}\n")
                f.write(f"Total libraries processed: {len(library_results)}\n")
                
                if valid_libraries:
                    f.write(f"Overall average across libraries: {overall_average:.6f}\n")
                    f.write(f"Standard deviation across libraries: {np.std(library_averages):.6f}\n")
                    f.write(f"Min library average: {np.min(library_averages):.6f}\n")
                    f.write(f"Max library average: {np.max(library_averages):.6f}\n")
                
                f.write(f"\nLibrary details:\n")
                for library_name, data in sorted(library_results.items()):
                    avg_score = f"{data['average_score']:.6f}" if data['average_score'] is not None else "N/A"
                    if args.top_n and data['successful_files'] > 0:
                        top_n_used = min(args.top_n, data['successful_files'])
                        f.write(f"{library_name}: {avg_score} (from top {top_n_used} of {data['successful_files']} files)\n")
                    else:
                        f.write(f"{library_name}: {avg_score} (from {data['successful_files']} files)\n")
                
                if args.verbose:
                    f.write(f"\nDetailed statistics:\n")
                    for library_name, data in sorted(library_results.items()):
                        if data['average_score'] is not None:
                            f.write(f"\n{library_name}:\n")
                            f.write(f"  Average: {data['average_score']:.6f}\n")
                            f.write(f"  Files: {data['successful_files']}\n")
                            f.write(f"  Excluded: {len(data['excluded_files'])}\n")
                            f.write(f"  Failed: {len(data['failed_files'])}\n")
                            
                            if data['individual_means']:
                                f.write(f"  Individual scores (sorted by score):\n")
                                sorted_means = sorted(data['individual_means'])
                                for i, mean_score in enumerate(sorted_means):
                                    marker = " *" if args.top_n and i < min(args.top_n, len(sorted_means)) else ""
                                    f.write(f"    File {i+1}: {mean_score:.6f}{marker}\n")
            
            print(f"\nResults saved to: {args.output}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

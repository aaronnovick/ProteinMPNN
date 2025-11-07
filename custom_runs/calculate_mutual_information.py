#!/usr/bin/env python3
"""
Calculate Mutual Information from Excel Dataset

This script calculates mutual information between pairs of continuous variables
in an Excel file. For continuous data, it uses binning to discretize the values
before calculating mutual information.

Usage:
    python3 calculate_mutual_information.py <excel_file> [options]
"""

import pandas as pd
import numpy as np
from sklearn.feature_selection import mutual_info_regression
from sklearn.metrics import mutual_info_score
import argparse
import sys
from itertools import combinations


def discretize_continuous(data, n_bins=10, strategy='uniform'):
    """
    Discretize continuous data into bins.
    
    Parameters:
    -----------
    data : array-like
        Continuous data to discretize
    n_bins : int
        Number of bins to use
    strategy : str
        'uniform' for uniform binning, 'quantile' for quantile-based binning
    
    Returns:
    --------
    discretized : array
        Discretized data with bin indices
    """
    data = np.array(data)
    if strategy == 'uniform':
        bins = np.linspace(data.min(), data.max(), n_bins + 1)
    elif strategy == 'quantile':
        bins = np.quantile(data, np.linspace(0, 1, n_bins + 1))
        bins[0] = data.min()  # Ensure first bin includes minimum
        bins[-1] = data.max()  # Ensure last bin includes maximum
    else:
        raise ValueError(f"Unknown strategy: {strategy}")
    
    discretized = np.digitize(data, bins) - 1
    # Handle edge case where value equals max
    discretized = np.clip(discretized, 0, n_bins - 1)
    return discretized


def calculate_mutual_info_pairwise(df, n_bins=10, strategy='uniform', min_samples=10):
    """
    Calculate mutual information for all pairs of columns in a DataFrame.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with continuous variables
    n_bins : int
        Number of bins for discretization
    strategy : str
        Binning strategy ('uniform' or 'quantile')
    min_samples : int
        Minimum number of samples required (will skip pairs with fewer)
    
    Returns:
    --------
    mi_matrix : pandas.DataFrame
        Symmetric matrix of mutual information values
    """
    columns = df.columns.tolist()
    n_cols = len(columns)
    mi_matrix = np.zeros((n_cols, n_cols))
    
    # Calculate MI for all pairs
    for i, col1 in enumerate(columns):
        for j, col2 in enumerate(columns):
            if i == j:
                # Self-MI is entropy (can be calculated, but typically not useful)
                mi_matrix[i, j] = np.nan
            else:
                # Get data where both columns have valid values
                pair_data = df[[col1, col2]].dropna()
                if len(pair_data) < min_samples:
                    mi_matrix[i, j] = np.nan
                else:
                    # Discretize the data for this pair
                    x = discretize_continuous(pair_data[col1].values, n_bins=n_bins, strategy=strategy)
                    y = discretize_continuous(pair_data[col2].values, n_bins=n_bins, strategy=strategy)
                    mi = mutual_info_score(x, y)
                    mi_matrix[i, j] = mi
    
    # Create DataFrame with column names
    mi_df = pd.DataFrame(mi_matrix, index=columns, columns=columns)
    return mi_df


def calculate_mutual_info_regression(df, target_col, n_bins=10, strategy='uniform'):
    """
    Calculate mutual information between a target column and all other columns.
    Uses regression-based MI which is designed for continuous target variables.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with continuous variables
    target_col : str
        Name of the target column
    n_bins : int
        Number of bins for discretization of features
    strategy : str
        Binning strategy
    
    Returns:
    --------
    mi_scores : pandas.Series
        Mutual information scores for each feature with the target
    """
    if target_col not in df.columns:
        raise ValueError(f"Target column '{target_col}' not found in DataFrame")
    
    # Prepare data
    feature_cols = [col for col in df.columns if col != target_col]
    data = df[[target_col] + feature_cols].dropna()
    
    if len(data) == 0:
        raise ValueError("No valid data after removing NaN values")
    
    y = data[target_col].values
    X = data[feature_cols].values
    
    # Discretize features
    X_discretized = np.zeros_like(X, dtype=int)
    for i, col in enumerate(feature_cols):
        X_discretized[:, i] = discretize_continuous(X[:, i], n_bins=n_bins, strategy=strategy)
    
    # Calculate MI using regression method
    mi_scores = mutual_info_regression(X_discretized, y, discrete_features=True, random_state=42)
    
    return pd.Series(mi_scores, index=feature_cols, name=f'MI_with_{target_col}')


def main():
    parser = argparse.ArgumentParser(
        description='Calculate mutual information from Excel dataset',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Calculate pairwise MI for all columns
  python3 calculate_mutual_information.py data.xlsx
  
  # Use quantile binning with 20 bins
  python3 calculate_mutual_information.py data.xlsx --n_bins 20 --strategy quantile
  
  # Calculate MI with respect to a specific target column
  python3 calculate_mutual_information.py data.xlsx --target_column "score"
  
  # Save results to CSV
  python3 calculate_mutual_information.py data.xlsx --output mi_results.csv
        """
    )
    
    parser.add_argument('excel_file', type=str, help='Path to Excel file')
    parser.add_argument('--sheet_name', type=str, default=0,
                        help='Sheet name or index (default: first sheet)')
    parser.add_argument('--n_bins', type=int, default=10,
                        help='Number of bins for discretization (default: 10)')
    parser.add_argument('--strategy', type=str, default='uniform',
                        choices=['uniform', 'quantile'],
                        help='Binning strategy: uniform or quantile (default: uniform)')
    parser.add_argument('--target_column', type=str, default=None,
                        help='If specified, calculate MI with respect to this column only')
    parser.add_argument('--output', type=str, default=None,
                        help='Output CSV file path (default: print to stdout)')
    parser.add_argument('--min_samples', type=int, default=10,
                        help='Minimum number of samples required (default: 10)')
    
    args = parser.parse_args()
    
    # Read Excel file
    try:
        print(f"Reading Excel file: {args.excel_file}")
        df = pd.read_excel(args.excel_file, sheet_name=args.sheet_name)
        print(f"Loaded {len(df)} rows and {len(df.columns)} columns")
        print(f"Columns: {', '.join(df.columns.tolist())}")
    except Exception as e:
        print(f"Error reading Excel file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Remove non-numeric columns
    numeric_df = df.select_dtypes(include=[np.number])
    if len(numeric_df.columns) < len(df.columns):
        non_numeric = set(df.columns) - set(numeric_df.columns)
        print(f"Warning: Removed non-numeric columns: {', '.join(non_numeric)}")
    
    if len(numeric_df.columns) < 2:
        print("Error: Need at least 2 numeric columns to calculate mutual information", file=sys.stderr)
        sys.exit(1)
    
    # Calculate mutual information
    if args.target_column:
        if args.target_column not in numeric_df.columns:
            print(f"Error: Target column '{args.target_column}' not found or not numeric", file=sys.stderr)
            sys.exit(1)
        
        print(f"\nCalculating MI with respect to target column: {args.target_column}")
        mi_results = calculate_mutual_info_regression(
            numeric_df, 
            args.target_column, 
            n_bins=args.n_bins, 
            strategy=args.strategy
        )
        mi_results = mi_results.sort_values(ascending=False)
        print("\nMutual Information Scores:")
        print(mi_results.to_string())
        
        if args.output:
            mi_results.to_csv(args.output)
            print(f"\nResults saved to: {args.output}")
    else:
        print(f"\nCalculating pairwise mutual information...")
        print(f"Using {args.n_bins} bins with {args.strategy} strategy")
        mi_matrix = calculate_mutual_info_pairwise(
            numeric_df, 
            n_bins=args.n_bins, 
            strategy=args.strategy,
            min_samples=args.min_samples
        )
        
        print("\nMutual Information Matrix:")
        print(mi_matrix.to_string())
        
        # Print summary statistics
        print("\nSummary Statistics:")
        upper_triangle = mi_matrix.where(
            np.triu(np.ones(mi_matrix.shape), k=1).astype(bool)
        )
        mi_values = upper_triangle.values.flatten()
        mi_values = mi_values[~np.isnan(mi_values)]
        
        if len(mi_values) > 0:
            print(f"  Mean MI: {np.mean(mi_values):.4f}")
            print(f"  Median MI: {np.median(mi_values):.4f}")
            print(f"  Min MI: {np.min(mi_values):.4f}")
            print(f"  Max MI: {np.max(mi_values):.4f}")
            
            # Find top pairs
            print("\nTop 10 Variable Pairs by Mutual Information:")
            pairs = []
            for i in range(len(mi_matrix.index)):
                for j in range(i + 1, len(mi_matrix.columns)):
                    val = mi_matrix.iloc[i, j]
                    if not np.isnan(val):
                        pairs.append((mi_matrix.index[i], mi_matrix.columns[j], val))
            
            pairs.sort(key=lambda x: x[2], reverse=True)
            for idx, (col1, col2, mi_val) in enumerate(pairs[:10], 1):
                print(f"  {idx}. {col1} <-> {col2}: {mi_val:.4f}")
        
        if args.output:
            mi_matrix.to_csv(args.output)
            print(f"\nResults saved to: {args.output}")


if __name__ == '__main__':
    main()


#!/bin/bash
# Example usage of the overlayed histogram analysis script

echo "=== Overlayed Histogram Analysis Examples ==="

# Example 1: Use default 5UOI sample directories
echo "Example 1: Analyzing all 5UOI samples (default)"
python3 overlayed_histogram_analysis.py \
    --output_dir ../outputs/my_variants/5UOI/overlayed_histogram_analysis \
    --title "5UOI All Variants Score Comparison"

echo ""

# Example 2: Specify specific sample directories
echo "Example 2: Analyzing specific sample directories"
python3 overlayed_histogram_analysis.py \
    --sample_dirs \
        ../outputs/my_variants/5UOI/sample_variant_scores/1/score_only \
        ../outputs/my_variants/5UOI/sample_variant_scores/2/score_only \
        ../outputs/my_variants/5UOI/sample_variant_scores/3/score_only \
    --output_dir ../outputs/my_variants/5UOI/overlayed_histogram_analysis_1_3 \
    --title "5UOI Variants 1-3 Score Comparison"

echo ""

# Example 3: Use NPZ file patterns
echo "Example 3: Using NPZ file patterns"
python3 overlayed_histogram_analysis.py \
    --npz_patterns \
        "../outputs/my_variants/5UOI/sample_variant_scores/1/score_only/*.npz" \
        "../outputs/my_variants/5UOI/sample_variant_scores/2/score_only/*.npz" \
    --output_dir ../outputs/my_variants/5UOI/overlayed_histogram_analysis_patterns \
    --title "5UOI Variants 1-2 Score Comparison (Pattern-based)"

echo ""
echo "=== Analysis Complete ==="
echo "Check the output directories for the generated overlayed histograms!" 
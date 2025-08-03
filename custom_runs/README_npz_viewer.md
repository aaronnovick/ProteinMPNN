# NPZ File Viewer

A Python script to view and analyze the contents of .npz files.

## Features

- Displays detailed information about arrays in NPZ files
- Shows file size, array shapes, data types, and element counts
- Provides statistical information for numeric arrays (min, max, mean, std, median)
- Shows sample values from arrays
- Supports multiple files and glob patterns
- Handles both numeric and non-numeric data

## Usage

### Basic Usage

```bash
# View a single NPZ file
python npz_viewer.py file.npz

# View multiple NPZ files
python npz_viewer.py file1.npz file2.npz file3.npz

# Use glob patterns to view multiple files
python npz_viewer.py *.npz
python npz_viewer.py outputs/**/*.npz
```

### Examples

```bash
# View a specific score file
python npz_viewer.py outputs/my_variants/5UOI/sample_variant_scores/1/score_only/5UOI_fasta_1.npz

# View all score files in a directory
python npz_viewer.py outputs/my_variants/5UOI/sample_variant_scores/1/score_only/*.npz

# View files from multiple directories
python npz_viewer.py outputs/example_*_outputs/score_only/*.npz
```

## Output Format

The script provides:

1. **File Information**: File path and size
2. **Array Summary**: Table showing all arrays with their shapes, data types, and sizes
3. **Detailed Analysis**: For each array:
   - Shape and data type
   - Memory usage
   - Statistical information (for numeric arrays)
   - Sample values (first/last few elements)

## Example Output

```
============================================================
NPZ File: outputs/my_variants/5UOI/sample_variant_scores/1/score_only/5UOI_fasta_1.npz
============================================================
File size: 1,556 bytes (1.52 KB)

Arrays in file:
Array Name           Shape           Data Type    Size      
-------------------- --------------- ------------ ----------
score                (5,)            float32      5
global_score         (5,)            float32      5
S                    (43,)           int64        43
seq_str              ()              <U43         1

Total elements: 54

============================================================
DETAILED ARRAY INFORMATION
============================================================

Array: 'score'
----------------------------------------
Shape: (5,)
Data type: float32
Size: 5 elements
Memory usage: 20 bytes (0.02 KB)
Min: 1.564774
Max: 1.705446
Mean: 1.627793
Std: 0.046428
Median: 1.616835
All values: [1.6079348 1.6168348 1.5647739 1.6439734 1.7054458]
```

## Requirements

- Python 3.6+
- NumPy

## Installation

No installation required. Just run the script directly:

```bash
python npz_viewer.py your_file.npz
```

## Notes

- The script automatically closes NPZ files after reading them
- For large arrays, only the first and last 5 values are shown
- Non-numeric arrays show unique values instead of statistics
- The script handles errors gracefully and provides informative error messages 
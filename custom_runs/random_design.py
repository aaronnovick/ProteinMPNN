import random
import os
import re

amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def read_fasta_sequences(fasta_file_path):
    """
    Purpose: 
        Read sequences from a FASTA file.
    Parameters:
        fasta_file_path: path to the FASTA file (str)
    Returns:
        list of tuples (header, sequence)
    """
    sequences = []
    with open(fasta_file_path, "r") as f:
        lines = f.readlines()
    
    current_header = None
    current_sequence = ""
    
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            # Save previous sequence if exists
            if current_header and current_sequence:
                sequences.append((current_header, current_sequence))
            # Start new sequence
            current_header = line
            current_sequence = ""
        else:
            current_sequence += line
    
    # Don't forget the last sequence
    if current_header and current_sequence:
        sequences.append((current_header, current_sequence))
    
    return sequences

def extract_sample_number(header):
    """
    Purpose:
        Extract sample number from sequence header.
    Parameters:
        header: sequence header (str)
    Returns:
        sample number (int) or None if not found
    """
    # Look for "sample=X" pattern
    match = re.search(r'sample=(\d+)', header)
    if match:
        return int(match.group(1))
    
    # If no sample number found, return None
    return None

def random_design(sequence, n_designs, positions):
    """
    Purpose: 
        Design a random sequence of length n_designs at the specified positions.
    Parameters:
        sequence: sequence to design (str)
        n_designs: the number of designs to generate (int)
        positions: the positions to design (list of ints)
    Returns:
        a fasta file of the sequences (str)
    """
    fasta_str = ""
    for i in range(n_designs):
        design = list(sequence)
        for pos in positions:
            design[pos] = random.choice(amino_acids)
        fasta_str += f">{i}\n{''.join(design)}\n"
    return fasta_str

def process_multiple_sequences(fasta_file_path, n_designs, positions_dict, base_output_dir):
    """
    Process multiple sequences from a FASTA file and generate random designs for each.
    Parameters:
        fasta_file_path: path to the input FASTA file (str)
        n_designs: number of designs to generate per sequence (int)
        positions_dict: dict mapping sample numbers to positions to randomize, or a single list for all samples
        base_output_dir: base directory for output files (str)
    """
    # Read sequences from FASTA file
    sequences = read_fasta_sequences(fasta_file_path)
    
    print(f"Found {len(sequences)} sequences in {fasta_file_path}")
    
    # Process each sequence
    for i, (header, sequence) in enumerate(sequences):
        print(f"Processing sequence {i+1}/{len(sequences)}: {header[:50]}...")
        
        # Extract sample number from header
        sample_num = extract_sample_number(header)
        
        # Determine which positions to use for this sample
        if isinstance(positions_dict, dict):
            # If positions_dict is a dictionary, look up positions by sample number
            if sample_num is not None and sample_num in positions_dict:
                positions = positions_dict[sample_num]
                print(f"  Using sample-specific positions for sample {sample_num}: {positions}")
            else:
                # Use default positions if sample not found in dict
                default_positions = positions_dict.get('default', [])
                positions = default_positions
                print(f"  Sample {sample_num} not found in positions_dict, using default positions: {positions}")
        else:
            # If positions_dict is a list, use the same positions for all samples
            positions = positions_dict
            print(f"  Using same positions for all samples: {positions}")
        
        if sample_num is not None:
            # Create sample-specific directory
            sample_dir = os.path.join(base_output_dir, str(sample_num))
            os.makedirs(sample_dir, exist_ok=True)
            
            # Output file path
            output_path = os.path.join(sample_dir, "random_designs.fa")
            
            print(f"  Sample number: {sample_num}")
            print(f"  Output directory: {sample_dir}")
        else:
            # For sequences without sample number (like the first one), use a default directory
            default_dir = os.path.join(base_output_dir, "original")
            os.makedirs(default_dir, exist_ok=True)
            output_path = os.path.join(default_dir, "random_designs.fa")
            
            print(f"  No sample number found, using default directory: {default_dir}")
        
        # Generate random designs for this sequence
        designs = random_design(sequence, n_designs, positions)
        
        # Write designs to file
        with open(output_path, "w") as f:
            # Write designs directly without comment lines for ProteinMPNN compatibility
            f.write(designs)
        
        print(f"  Saved {n_designs} designs to {output_path}")

if __name__ == "__main__":
    # Configuration
    fasta_file_path = "../outputs/my_designs/sequences/seqs/5UOI.fa"
    n_designs = 100
    
    # Define position lists once and reference them by name
    position_lists = {
        '1': [18, 19, 20, 22, 23, 26, 27],
        '2': [18, 19, 20, 23, 26, 27],
        '3': [18, 19, 23, 26, 27],
        '4': [18, 19, 23, 27],
        '5': [18, 19, 22, 23, 26, 27],
        '6': [19, 23, 26, 27],
        '7': [19, 20, 23, 26, 27],
        '8': [18, 19, 26, 27],
        '9': [18, 19, 20, 23, 27],
        '10': [18, 19, 23, 26],
        '11': [18, 19, 20, 23, 26],
        '12': [18, 19, 20, 26, 27]
    }
    
    # Map samples to position lists (can reference the same list multiple times)
    positions_dict = {
        1: position_lists['1'],
        2: position_lists['2'],
        3: position_lists['2'],
        4: position_lists['2'],
        5: position_lists['2'],
        6: position_lists['3'],
        7: position_lists['2'],
        8: position_lists['4'],
        9: position_lists['5'],
        10: position_lists['5'],
        11: position_lists['3'],
        12: position_lists['3'],
        13: position_lists['3'],
        14: position_lists['2'],
        15: position_lists['2'],
        16: position_lists['6'],
        17: position_lists['2'],
        18: position_lists['5'],
        19: position_lists['8'],
        20: position_lists['7'],
        21: position_lists['2'],
        22: position_lists['8'],
        23: position_lists['9'],
        24: position_lists['1'],
        25: position_lists['2'],
        26: position_lists['2'],
        27: position_lists['3'],
        28: position_lists['3'],
        29: position_lists['5'],
        30: position_lists['10'],
        31: position_lists['2'],
        32: position_lists['2'],
        33: position_lists['11'],
        34: position_lists['9'],
        35: position_lists['2'],
        36: position_lists['3'],
        37: position_lists['2'],
        38: position_lists['3'],
        39: position_lists['3'],
        40: position_lists['2'],
        41: position_lists['8'],
        42: position_lists['2'],
        43: position_lists['3'],
        44: position_lists['2'],
        45: position_lists['2'],
        46: position_lists['12'],
        47: position_lists['11'],
        48: position_lists['3'],
        49: position_lists['2'],
        50: position_lists['3'],
        # Add more samples as needed
    }
    
    base_output_dir = "../outputs/my_variants/5UOI/sample_random_variants"
    
    # Process all sequences
    process_multiple_sequences(fasta_file_path, n_designs, positions_dict, base_output_dir)
    
    print(f"\nCompleted! All random designs saved to their respective sample directories in {base_output_dir}")
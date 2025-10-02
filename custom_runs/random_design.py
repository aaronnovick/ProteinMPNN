import random
import os
import re
from typing import List, Tuple

AMINO_ACIDS = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
]

def read_fasta_sequences(fasta_file_path: str) -> List[Tuple[str, str]]:
    """
    Read sequences from a FASTA file.

    Args:
        fasta_file_path: Path to the FASTA file.

    Returns:
        List of (header, sequence) tuples.
    """
    sequences = []
    with open(fasta_file_path, "r") as f:
        lines = f.readlines()

    current_header, current_sequence = None, ""
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if current_header and current_sequence:
                sequences.append((current_header, current_sequence))
            current_header, current_sequence = line, ""
        else:
            current_sequence += line

    if current_header and current_sequence:
        sequences.append((current_header, current_sequence))

    return sequences

def extract_sample_number(header: str) -> int | None:
    """
    Extract sample number from a FASTA header.

    Args:
        header: FASTA header string.

    Returns:
        Sample number if found, else None.
    """
    # First try to find sample=X pattern (for 5UOI format)
    match = re.search(r'sample=(\d+)', header)
    if match:
        return int(match.group(1))
    
    # If no sample=X pattern, try to extract a number from the header
    # This handles formats like "HHH_rd4_0169\" by extracting the last number
    numbers = re.findall(r'\d+', header)
    if numbers:
        return int(numbers[-1])  # Use the last number found
    
    return None

def random_design(sequence: str, n_designs: int, positions: List[int], header_prefix: str = "", 
                 position_constraints: dict = None) -> str:
    """
    Generate random designs for a sequence.

    Args:
        sequence: Input sequence.
        n_designs: Number of designs to generate.
        positions: List of positions (0-based) to randomize.
        header_prefix: Prefix to add in FASTA headers.
        position_constraints: Dict mapping position (0-based) to allowed amino acids.
                             If None, all positions use all 20 amino acids.

    Returns:
        FASTA-formatted string of designs.
    """
    fasta_str = ""
    for i in range(n_designs):
        design = list(sequence)
        for pos in set(positions):  # deduplicate positions
            if 0 <= pos < len(sequence):
                # Use position-specific constraints if provided, otherwise use all amino acids
                if position_constraints and pos in position_constraints:
                    allowed_aas = position_constraints[pos]
                    if not allowed_aas:  # Empty list means no constraints (use all)
                        allowed_aas = AMINO_ACIDS
                else:
                    allowed_aas = AMINO_ACIDS
                
                design[pos] = random.choice(allowed_aas)
        fasta_str += f">{header_prefix}design{i}\n{''.join(design)}\n"
    return fasta_str

def process_multiple_sequences(fasta_file_path: str, n_designs: int, positions: List[int], base_output_dir: str, 
                              position_constraints: dict = None, skip_native: bool = True) -> None:
    """
    Process multiple sequences and generate random designs.

    Args:
        fasta_file_path: Path to the input FASTA file.
        n_designs: Number of designs per sequence.
        positions: List of positions (0-based) to randomize for all samples.
        base_output_dir: Base directory for output.
        position_constraints: Dict mapping position (0-based) to allowed amino acids.
        skip_native: If True, skip the first sequence (assumed to be native/wild type).
    """
    sequences = read_fasta_sequences(fasta_file_path)
    print(f"Found {len(sequences)} sequences in {fasta_file_path}")
    
    # Skip the first sequence (native/wild type) if requested
    if skip_native and len(sequences) > 1:
        sequences = sequences[1:]
        print(f"Skipping native sequence. Processing {len(sequences)} design sequences.")
    elif skip_native and len(sequences) == 1:
        print("Warning: Only one sequence found and skip_native=True. No sequences to process.")
        return

    for i, (header, sequence) in enumerate(sequences):
        print(f"Processing sequence {i+1}/{len(sequences)}: {header[:50]}...")

        sample_num = extract_sample_number(header)

        if sample_num is not None:
            sample_dir = os.path.join(base_output_dir, str(sample_num))
            os.makedirs(sample_dir, exist_ok=True)
            output_path = os.path.join(sample_dir, "random_designs.fa")
            header_prefix = f"sample{sample_num}_"
        else:
            default_dir = os.path.join(base_output_dir, "original")
            os.makedirs(default_dir, exist_ok=True)
            output_path = os.path.join(default_dir, "random_designs.fa")
            header_prefix = "original_"

        designs = random_design(sequence, n_designs, positions, header_prefix, position_constraints)

        with open(output_path, "w") as f:
            f.write(designs)

        print(f"  Saved {n_designs} designs to {output_path}")

if __name__ == "__main__":
    fasta_file_path = "../McConnell_variants/aaa/inputs/parents.fa"
    n_designs = 100

    # Define PDB positions
    positions_pdb = [18, 19, 22, 26, 33, 36, 37, 40]

    # Warn if duplicates exist
    if len(positions_pdb) != len(set(positions_pdb)):
        dupes = [p for p in positions_pdb if positions_pdb.count(p) > 1]
        print(f"⚠️ Warning: duplicate positions found in PDB list: {sorted(set(dupes))}")
    
    # Convert to 0-based for Python
    positions = [p - 1 for p in positions_pdb]

    # OPTIONAL: Define position-specific amino acid constraints
    # Format: {position_0_based: [list_of_allowed_amino_acids]}
    # Example: Only allow hydrophobic amino acids at certain positions
    position_constraints = {
        # Examples below:
        # 16: ["A", "V", "L", "I", "M", "F", "W", "Y"],  # Position 17 (0-based 16): only hydrophobic
        # 17: ["K", "R", "H"],  # Position 18 (0-based 17): only basic
        # 18: ["D", "E"],  # Position 19 (0-based 18): only acidic
        # Leave other positions unconstrained (will use all 20 amino acids)
    }
    
    # If you want to use constraints, uncomment the lines above and modify as needed
    # If position_constraints is empty or None, all positions will use all 20 amino acids

    # Skip native sequence (first sequence) when generating variants
    # Set skip_native=False if you want to include the native sequence
    base_output_dir = "../McConnell_variants/aaa/outputs/helix23/variants"
    process_multiple_sequences(fasta_file_path, n_designs, positions, base_output_dir, position_constraints, skip_native=True)
    print(f"\nCompleted! All random designs saved to {base_output_dir}")

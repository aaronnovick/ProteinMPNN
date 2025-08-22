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
    match = re.search(r'sample=(\d+)', header)
    return int(match.group(1)) if match else None

def random_design(sequence: str, n_designs: int, positions: List[int], header_prefix: str = "") -> str:
    """
    Generate random designs for a sequence.

    Args:
        sequence: Input sequence.
        n_designs: Number of designs to generate.
        positions: List of positions (0-based) to randomize.
        header_prefix: Prefix to add in FASTA headers.

    Returns:
        FASTA-formatted string of designs.
    """
    fasta_str = ""
    for i in range(n_designs):
        design = list(sequence)
        for pos in set(positions):  # deduplicate positions
            if 0 <= pos < len(sequence):
                design[pos] = random.choice(AMINO_ACIDS)
        fasta_str += f">{header_prefix}design{i}\n{''.join(design)}\n"
    return fasta_str

def process_multiple_sequences(fasta_file_path: str, n_designs: int, positions: List[int], base_output_dir: str) -> None:
    """
    Process multiple sequences and generate random designs.

    Args:
        fasta_file_path: Path to the input FASTA file.
        n_designs: Number of designs per sequence.
        positions: List of positions (0-based) to randomize for all samples.
        base_output_dir: Base directory for output.
    """
    sequences = read_fasta_sequences(fasta_file_path)
    print(f"Found {len(sequences)} sequences in {fasta_file_path}")

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

        designs = random_design(sequence, n_designs, positions, header_prefix)

        with open(output_path, "w") as f:
            f.write(designs)

        print(f"  Saved {n_designs} designs to {output_path}")

if __name__ == "__main__":
    fasta_file_path = "../outputs/my_designs/sequences/seqs/5UOI.fa"
    n_designs = 100

    # Define PDB positions
    positions_pdb = [17, 18, 19, 22, 27, 29]

    # Warn if duplicates exist
    if len(positions_pdb) != len(set(positions_pdb)):
        dupes = [p for p in positions_pdb if positions_pdb.count(p) > 1]
        print(f"⚠️ Warning: duplicate positions found in PDB list: {sorted(set(dupes))}")
    
    # Convert to 0-based for Python
    positions = [p - 1 for p in positions_pdb]

    base_output_dir = "../outputs/my_variants/5UOI/sample_random_variants"
    process_multiple_sequences(fasta_file_path, n_designs, positions, base_output_dir)
    print(f"\nCompleted! All random designs saved to {base_output_dir}")

import random
amino_acids = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

new_file_path = "../outputs/my_variants/5UOI/sample_random_variants/5/random_designs.fa"

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

if __name__ == "__main__":
    sequence = "KSYEEIAKKLLEKYDVEEEVALRAVKEAGGDLEKAEKLVREPL"
    n_designs = 100
    positions = [18, 19, 23, 26, 27]
    designs = random_design(sequence, n_designs, positions)
    with open(new_file_path, "w") as f:
        f.write(designs)
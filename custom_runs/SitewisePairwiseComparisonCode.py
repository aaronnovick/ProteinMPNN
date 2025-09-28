import pandas as pd


def calculate_pairwise_identity(sequences):
    """
    Calculates the pairwise percent identity for a list of sequences.

    Args:
        sequences (list): A list of strings, where each string is a protein sequence.

    Returns:
        tuple: A tuple containing a pandas DataFrame with the results and a formatted
               string of the results.
    """
    if not all(len(s) == len(sequences[0]) for s in sequences):
        raise ValueError("All sequences must be of the same length.")

    num_seqs = len(sequences)
    identity_matrix = pd.DataFrame(index=range(num_seqs), columns=range(num_seqs))

    for i in range(num_seqs):
        for j in range(i, num_seqs):
            seq1 = sequences[i]
            seq2 = sequences[j]
            length = len(seq1)
            
            # Count identical residues
            matches = sum(c1 == c2 for c1, c2 in zip(seq1, seq2))
            
            # Calculate percent identity
            percent_identity = (matches / length) * 100
            
            # Populate the matrix symmetrically
            identity_matrix.loc[i, j] = percent_identity
            identity_matrix.loc[j, i] = percent_identity
    
    # Format the DataFrame for better readability
    seq_labels = [f"Seq {i+1}" for i in range(num_seqs)]
    identity_matrix.columns = seq_labels
    identity_matrix.index = seq_labels
    identity_matrix = identity_matrix.astype(float).round(2)
    
    # Prepare the formatted table string
    table_string = identity_matrix.to_string()
    
    return identity_matrix, table_string


def read_fasta_sequences(fasta_path):
    """Read sequences from a FASTA file and return a list of sequences.

    Lines beginning with '>' are treated as headers; subsequent non-header lines
    are concatenated as the sequence. Empty lines are ignored.
    """
    sequences = []
    current_seq_parts = []

    with open(fasta_path, 'r') as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq_parts:
                    sequences.append(''.join(current_seq_parts))
                    current_seq_parts = []
            else:
                current_seq_parts.append(line)
        # flush last
        if current_seq_parts:
            sequences.append(''.join(current_seq_parts))

    return sequences


# === INPUT SOURCE ===
# Read sequences from the generated FASTA file instead of manual input.
FASTA_PATH = '/Users/aaronnovick/Desktop/Research_Project/ProteinMPNN/outputs/my_designs/sequences/seqs/5UOI.fa'
all_sequences = read_fasta_sequences(FASTA_PATH)

# Exclude the first sequence (wild type/native)
sequences = all_sequences[1:] if len(all_sequences) > 1 else []

# === RUN THE CODE ===
try:
    if not sequences:
        raise ValueError(f"No design sequences found after excluding wild type in FASTA: {FASTA_PATH}")

    results_df, results_table = calculate_pairwise_identity(sequences)
    
    # 1. Print the results in an easy-to-read table
    print("Pairwise Percent Identity (%)")
    print(results_table)
    print("\n")
    
    # 2. Save the results to an Excel spreadsheet
    excel_filename = '5UOI_sequence_identity_results.xlsx'
    results_df.to_excel(excel_filename, index=True)
    print(f"Results successfully saved to {excel_filename}")

except ValueError as e:
    print(f"Error: {e}")
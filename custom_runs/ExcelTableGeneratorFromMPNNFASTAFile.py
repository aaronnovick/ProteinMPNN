import pandas as pd
import re
import os

def parse_proteinmpnn_fa_to_excel(file_path, output_excel_name):
    """
    Parses a ProteinMPNN .fa file, extracts sequence data and metadata,
    and saves it to an Excel spreadsheet. Skips the first WT sequence.

    Args:
        file_path (str): Path to the .fa file.
        output_excel_name (str): Name of the output Excel file.
    """
    data = []
    
    # Check if the file exists before attempting to open it
    if not os.path.exists(file_path):
        print(f"Error: File not found at '{file_path}'")
        return

    with open(file_path, 'r') as f:
        # Use a state variable to skip the first WT sequence
        skip_wt = True
        
        # Read the file line by line
        lines = f.readlines()
        
        # Iterate through the lines two at a time (header and sequence)
        for i in range(0, len(lines), 2):
            header = lines[i].strip()
            # Ensure the line starts with a '>' and the next line is not a header
            if not header.startswith('>') or i + 1 >= len(lines) or lines[i+1].startswith('>'):
                continue

            # Skip the first (WT) sequence
            if skip_wt:
                skip_wt = False
                continue
            
            sequence = lines[i+1].strip()

            # Use regular expressions to extract metadata
            t_match = re.search(r'T=(\d+\.?\d*)', header)
            sample_match = re.search(r'sample=(\d+)', header)
            score_match = re.search(r'score=(\d+\.?\d*)', header)
            global_score_match = re.search(r'global_score=(\d+\.?\d*)', header)
            seq_recovery_match = re.search(r'seq_recovery=(\d+\.?\d*)', header)

            # Extract values, defaulting to None if not found
            T = float(t_match.group(1)) if t_match else None
            sample = int(sample_match.group(1)) if sample_match else None
            score = float(score_match.group(1)) if score_match else None
            global_score = float(global_score_match.group(1)) if global_score_match else None
            seq_recovery = float(seq_recovery_match.group(1)) if seq_recovery_match else None
            
            # Append the data to the list
            data.append({
                'T': T,
                'Sample': sample,
                'Score': score,
                'Global_Score': global_score,
                'Seq_Recovery': seq_recovery,
                'Sequence': sequence
            })

    # Create a Pandas DataFrame
    df = pd.DataFrame(data)

    # Save the DataFrame to an Excel file
    df.to_excel(output_excel_name, index=False)
    
    print(f"Successfully created '{output_excel_name}' with {len(df)} entries.")

# === EXAMPLE USAGE ===

# Set the base directory for the outputs
base_path = 'C:\\Users\\claud\\ProteinMPNN_Folder\\New_ProteinMPNN\\ProteinMPNN\\output\\seqs'

# Set the customizable part of the filename
distance_cutoff = 'cys_restriction'  # Change this to '8A', '4A', etc.

# Construct the full file paths
fa_file_name = f'1rbp_trimmed_{distance_cutoff}.fa'
fa_file_path = os.path.join(base_path, fa_file_name)
excel_file_name = f'1rbp_trimmed_{distance_cutoff}.xlsx'

parse_proteinmpnn_fa_to_excel(fa_file_path, excel_file_name)
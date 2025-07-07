#!/bin/bash

# Set paths (adjust as needed)
input_pdb_dir="../inputs/my_structures"
output_dir="../outputs/my_designs"
parsed_jsonl="$output_dir/parsed_pdbs.jsonl"
seq_output_dir="$output_dir/sequences"

# Make output dirs if they don't exist
mkdir -p "$output_dir"
mkdir -p "$seq_output_dir"

# Step 1: Parse the structure(s)
python3 ../helper_scripts/parse_multiple_chains.py \
    --input_path="$input_pdb_dir" \
    --output_path="$parsed_jsonl"

# Step 2: Run ProteinMPNN
python3 ../protein_mpnn_run.py \
    --jsonl_path="$parsed_jsonl" \
    --out_folder="$seq_output_dir" \
    --num_seq_per_target 5 \
    --sampling_temp "0.1" \
    --seed 42 \
    --batch_size 1
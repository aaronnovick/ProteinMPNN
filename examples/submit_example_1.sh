#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH --gres=gpu:rtx2080:1
#SBATCH -c 2
#SBATCH --output=example_1.out

source ~/.bash_profile #my purpose is to stop the annoying conda init vs conda source error that is purely cosmetic
echo $SHELL #this is just to print that we're using bash instead of zsh for running commands, can comment out

folder_with_pdbs="../inputs/PDB_monomers/pdbs/"

output_dir="../outputs/example_1_outputs"
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"

python3 ../helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python3 ../protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --out_folder $output_dir \
        --num_seq_per_target 2 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1

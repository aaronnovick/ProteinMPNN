#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH --gres=gpu:rtx2080:1
#SBATCH -c 3
#SBATCH --output=example_3.out

#source activate mlfold
source ~/.bash_profile
echo $SHELL

path_to_PDB="../inputs/PDB_complexes/pdbs/3HTN.pdb"

output_dir="../outputs/example_3_outputs"
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

chains_to_design="A B"

python3 ../protein_mpnn_run.py \
        --pdb_path $path_to_PDB \
        --pdb_path_chains "$chains_to_design" \
        --out_folder $output_dir \
        --num_seq_per_target 2 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1

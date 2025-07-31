#!/bin/bash
#SBATCH -p gpu
#SBATCH --mem=32g
#SBATCH --gres=gpu:rtx2080:1
#SBATCH -c 3
#SBATCH --output=example_3_from_fasta.out

#source activate mlfold
source ~/.bash_profile
echo $SHELL
conda activate mlfold

path_to_PDB="../inputs/my_structures/5UOI/5UOI.pdb"
path_to_fasta="../outputs/my_variants/5UOI/sample_random_variants/1/random_designs.fa"

output_dir="../outputs/my_variants/5UOI/sample_variant_scores/1"
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

chains_to_design="A"

python3 ../protein_mpnn_run.py \
        --path_to_fasta $path_to_fasta \
        --pdb_path $path_to_PDB \
        --pdb_path_chains "$chains_to_design" \
        --out_folder $output_dir \
        --num_seq_per_target 5 \
        --sampling_temp "0.1" \
        --score_only 1 \
        --seed 42 \
        --batch_size 1

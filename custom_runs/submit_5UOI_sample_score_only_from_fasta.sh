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
base_variants_dir="../outputs/my_variants/5UOI/sample_random_variants"
base_output_dir="../outputs/my_variants/5UOI/sample_variant_scores"
chains_to_design="A"

# Loop through all sample directories
for sample_dir in "$base_variants_dir"/*/; do
    if [ -d "$sample_dir" ]; then
        # Extract sample number from directory name
        sample_name=$(basename "$sample_dir")
        
        # Skip if it's not a number (like "original" directory)
        if [[ "$sample_name" =~ ^[0-9]+$ ]]; then
            echo "Processing sample: $sample_name"
            
            # Define paths for this sample (remove trailing slash to avoid double slashes)
            sample_dir_clean="${sample_dir%/}"
            path_to_fasta="$sample_dir_clean/random_designs.fa"
            output_dir="$base_output_dir/$sample_name"
            
            # Create output directory if it doesn't exist
            if [ ! -d "$output_dir" ]; then
                mkdir -p "$output_dir"
            fi
            
            # Check if the fasta file exists
            if [ -f "$path_to_fasta" ]; then
                echo "  Processing: $path_to_fasta"
                echo "  Output: $output_dir"
                
                # Check if fasta file is not empty
                if [ -s "$path_to_fasta" ]; then
                    # Create a clean fasta file without comment lines for ProteinMPNN
                    clean_fasta="$output_dir/clean_random_designs.fa"
                    grep -v '^#' "$path_to_fasta" > "$clean_fasta"
                    
                    python3 ../protein_mpnn_run.py \
                        --path_to_fasta "$clean_fasta" \
                        --pdb_path "$path_to_PDB" \
                        --pdb_path_chains "$chains_to_design" \
                        --out_folder "$output_dir" \
                        --num_seq_per_target 10 \
                        --sampling_temp "0.1" \
                        --score_only 1 \
                        --seed 42 \
                        --batch_size 1
                    
                    echo "  Completed sample $sample_name"
                else
                    echo "  Warning: Fasta file is empty for sample $sample_name: $path_to_fasta"
                fi
            else
                echo "  Warning: Fasta file not found for sample $sample_name: $path_to_fasta"
            fi
            
            echo "----------------------------------------"
        else
            echo "Skipping non-numeric directory: $sample_name"
        fi
    fi
done

echo "All samples processed!"

#!/bin/bash
#SBATCH --job-name=csgs
#SBATCH --partition=schultzlab
#SBATCH --output=/work/hs325/csgs2026/src/logs/csgs.out
#SBATCH --error=/work/hs325/csgs2026/src/logs/csgs.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hs325@duke.edu

source /hpc/group/schultzlab/hs325/miniconda3/etc/profile.d/conda.sh
conda activate gsAI

# run scripts here TBA
python3 01_hpt.py \
    -o "output directory from 01" \
    -f "data file path" \
    -g "generation_choice" 

python3 02_crossval.py \
    -i "output directory from 01" \
    -o "output directory path" \
    -f "data file path" \
    -g "generation_choice" 


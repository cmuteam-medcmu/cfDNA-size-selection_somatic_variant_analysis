#!/bin/bash

#SBATCH --job-name=TLOD
#SBATCH --account=o250002
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --mem=80G
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# =========================
# Example usage:
# =========================
# TNM_VCF_LIST=/path/to/vcf_TNM_file_list.txt \
# TO_VCF_LIST=/path/to/vcf_TO_file_list.txt \
# OUTPUT_DIR=/path/to/output_directory \
# bash run_sweep_TLOD.sh

module load python/3.12.7

python sweep_TLOD.py
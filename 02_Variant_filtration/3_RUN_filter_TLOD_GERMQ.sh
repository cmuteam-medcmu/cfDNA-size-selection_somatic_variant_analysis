#!/bin/bash

#SBATCH --job-name=filterVAR
#SBATCH --account=o250002
#SBATCH --partition=short
#SBATCH --time=4:00:00
#SBATCH --mem=80G
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# =========================
# Example usage:
# =========================
# VCF_LIST=/path/to/vcf_TO_file_list.txt \
# OUTPUT_DIR=/path/to/output_directory \
# bash run_filter_TLOD_GERMQ.sh

module load python/3.12.7

python filter_TLOD_GERMQ.py
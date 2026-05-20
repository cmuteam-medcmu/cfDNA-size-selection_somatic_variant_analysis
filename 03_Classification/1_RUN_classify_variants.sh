#!/bin/bash

#SBATCH --job-name=classify_variants
#SBATCH --account=<project-account>
#SBATCH --partition=short
#SBATCH --time=04:00:00
#SBATCH --mem=100G
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

module load python/3.12.7

# =========================================================
# Example usage
# =========================================================
#
# export VCF_DIR="/path/to/annotated_vcfs"
# export OUT_ROOT="/path/to/output_dir"
# export COSMIC_TSV="/path/to/COSMIC.tsv"
#
# python classify_variants.py
#
# =========================================================

export VCF_DIR="/path/to/annotated_vcfs"

export OUT_ROOT="/path/to/output_dir"

export COSMIC_TSV="/path/to/CancerMutationCensus_AllData_v99_GRCh38.tsv"

python classify_variants.py
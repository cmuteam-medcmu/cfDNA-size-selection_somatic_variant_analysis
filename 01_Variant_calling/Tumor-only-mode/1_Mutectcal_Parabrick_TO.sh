#!/bin/bash

#SBATCH --job-name=parabrick_TO
#SBATCH --account=o250002
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=47
#SBATCH --mem=188G
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# =========================
# Example usage:
# =========================
# TUMOR_ID=Sample_T \
# TUMOR_BAM=/data/tumor.bam \
# REF=/ref/GRCh38.fa \
# PON=/ref/pon.vcf.gz \
# GERMLINE=/ref/gnomad.vcf.gz \
# INTERVAL=/ref/interval.list \
# OUT_VCF=output.vcf \
# bash tumor_only.sh

module load apptainer

TUMOR_ID=${TUMOR_ID}

REF=${REF}
TUMOR_BAM=${TUMOR_BAM}
PON=${PON}
GERMLINE=${GERMLINE}
INTERVAL=${INTERVAL}
OUT_VCF=${OUT_VCF}

apptainer exec --bind /common \
    --bind /project/o250002_cfVAR/ \
    --nv \
    /common/sif/clara-parabricks/4.4.0.sif pbrun \
        prepon \
            --in-pon-file $PON

apptainer exec --bind /common \
    --bind /project/o250002_cfVAR/ \
    --nv \
    /common/sif/clara-parabricks/4.4.0.sif pbrun \
        mutectcaller \
            --ref $REF \
            --tumor-name ${TUMOR_ID} \
            --in-tumor-bam $TUMOR_BAM \
            --pon $PON \
            --mutect-germline-resource $GERMLINE \
            --interval-file $INTERVAL \
            --out-vcf $OUT_VCF

apptainer exec --bind /common \
    --bind /project/o250002_cfVAR/ \
    --nv \
    /common/sif/clara-parabricks/4.4.0.sif pbrun \
        postpon \
            --in-vcf $OUT_VCF \
            --in-pon-file $PON \
            --out-vcf ${OUT_VCF%.vcf}_withpon.vcf
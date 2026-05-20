#!/bin/bash

#SBATCH --job-name=variantcall_1
#SBATCH --account=o250002
#SBATCH --partition=compute
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=40
#SBATCH --mem=100G
#SBATCH --output=output_%j.log
#SBATCH --error=error_%j.log

# =========================
# Example usage:
# =========================
# TUMOR=Sample_T \
# TUMOR_ID=Sample_T \
# NORMAL_ID=Sample_N \
# BAM_TUMOR=/data/tumor.bam \
# BAM_NORMAL=/data/normal.bam \
# REF=/ref/GRCh38.fa \
# GERMLINE=/ref/gnomad.vcf.gz \
# V_GERM=/ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
# INTERVAL=/ref/wgs_calling_regions.interval_list \
# OUT_PATH=/path/to/output_directory \
# bash tumor_normal_filtermutectcalls.sh

## Variables
ref=${REF}
germline_resource=${GERMLINE}
V_germ_resource=${V_GERM}
interval=${INTERVAL}

out_path=${OUT_PATH}

tumor=${TUMOR}
tumor_ID=${TUMOR_ID}
normal_ID=${NORMAL_ID}

bam_tumor=${BAM_TUMOR}
bam_normal=${BAM_NORMAL}

memory=100
cpu=40

module load python
module load gatk/4.6.1.0
module load java-jdk/21.0.3+9

## Estimate contamination in tumor sample
gatk --java-options "-Xms8G -Xmx${memory}G -XX:ParallelGCThreads=${cpu}" GetPileupSummaries \
    -I ${bam_tumor} \
    -V ${V_germ_resource} \
    -L ${interval} \
    -O ${out_path}/${tumor}_getpileupsummaries-Tumor.table 2>&1 >> ${out_path}/${tumor}_getpileup-report-Tumor.txt

## Estimate contamination in matched normal sample
gatk --java-options "-Xms8G -Xmx${memory}G -XX:ParallelGCThreads=${cpu}" GetPileupSummaries \
    -I ${bam_normal} \
    -V ${V_germ_resource} \
    -L ${interval} \
    -O ${out_path}/${tumor}_getpileupsummaries-Normal.table 2>&1 >> ${out_path}/${tumor}_getpileup-report-Normal.txt

## Calculate contamination
gatk --java-options "-Xms8G -Xmx${memory}G -XX:ParallelGCThreads=${cpu}" CalculateContamination \
    -I ${out_path}/${tumor}_getpileupsummaries-Tumor.table \
    --matched ${out_path}/${tumor}_getpileupsummaries-Normal.table \
    -O ${out_path}/${tumor}_calculatecontamination.table 2>&1 >> ${out_path}/${tumor}_calcontaminate-report.txt

## Filter Mutect2 calls
gatk --java-options "-Xms8G -Xmx${memory}G -XX:ParallelGCThreads=${cpu}" FilterMutectCalls \
    -R ${ref} \
    -V ${out_path}/${tumor}_TNM_withpon.vcf \
    --contamination-table ${out_path}/${tumor}_calculatecontamination.table \
    --stats ${out_path}/${tumor}_TNM.vcf.stats \
    -O ${out_path}/${tumor}_TNM_somatic.vcf 2>&1 >> ${out_path}/${tumor}_TNM_filter-report.txt

## Select PASS variants only
head -n 1 ${out_path}/${tumor}_TNM_somatic.vcf > ${out_path}/${tumor}_TNM_somatic_PASS.vcf

grep -v '##' ${out_path}/${tumor}_TNM_somatic.vcf | grep '#' | \
    sed "s/${tumor_ID}/tumor/g" | \
    sed "s/${normal_ID}/normal/g" >> ${out_path}/${tumor}_TNM_somatic_PASS.vcf

grep -v '#' ${out_path}/${tumor}_TNM_somatic.vcf | \
    grep "PASS" >> ${out_path}/${tumor}_TNM_somatic_PASS.vcf
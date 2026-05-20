#!/bin/bash

#SBATCH --job-name=variantcall_TO
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
# BAM=/data/tumor.bam \
# REF=/ref/GRCh38.fa \
# GERMLINE=/ref/gnomad.vcf.gz \
# V_GERM=/ref/1000G.vcf.gz \
# INTERVAL=/ref/interval.list \
# OUT_PATH=/output/dir \
# bash tumor_only_filter.sh

## Variables
ref=${REF}
germline_resource=${GERMLINE}
V_germ_resource=${V_GERM}
interval=${INTERVAL}
out_path=${OUT_PATH}
tumor=${TUMOR}
bam=${BAM}
memory=100
cpu=40

module load gatk/4.6.1.0
module load java-jdk/21.0.3+9

## Estimate Contamination
gatk --java-options "-Xms8G -Xmx${memory}G -XX:ParallelGCThreads=${cpu}" GetPileupSummaries \
    -I ${bam} \
    -V ${V_germ_resource} \
    -L ${interval} \
    -O ${out_path}/${tumor}_getpileupsummaries.table 2>&1 >> ${out_path}/${tumor}_getpileup-report.txt 

gatk --java-options "-Xms8G -Xmx${memory}G -XX:ParallelGCThreads=${cpu}" CalculateContamination \
    -I ${out_path}/${tumor}_getpileupsummaries.table \
    -O ${out_path}/${tumor}_calculatecontamination.table 2>&1 >> ${out_path}/${tumor}_calcontaminate-report.txt 

## Filter Mutect2 Calls
gatk --java-options "-Xms8G -Xmx${memory}G -XX:ParallelGCThreads=${cpu}" FilterMutectCalls \
    -R ${ref} \
    -V ${out_path}/${tumor}_withpon.vcf \
    --contamination-table ${out_path}/${tumor}_calculatecontamination.table \
    --stats ${out_path}/${tumor}.vcf.stats \
    -O ${out_path}/${tumor}_somatic.vcf 2>&1 >> ${out_path}/${tumor}_filter-report.txt 

## Select PASS variants
head -n 1 ${out_path}/${tumor}_somatic.vcf > ${out_path}/${tumor}_somatic_PASS.vcf
grep -v '##' ${out_path}/${tumor}_somatic.vcf | grep '#' >> ${out_path}/${tumor}_somatic_PASS.vcf
grep -v '#' ${out_path}/${tumor}_somatic.vcf | grep "PASS" >> ${out_path}/${tumor}_somatic_PASS.vcf
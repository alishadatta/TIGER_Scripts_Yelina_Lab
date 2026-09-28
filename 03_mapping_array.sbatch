#!/bin/bash
#SBATCH --job-name=map500
#SBATCH --array=1-500
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=mapping_logs/%A_%a.out
#SBATCH --error=mapping_logs/%A_%a.err

source ~/.bashrc

# BWA 0.7.17 from the GitHub pipeline environment
conda activate tiger_mapping

# Put working samtools 1.23 from bcftools_env on PATH
export PATH="/users/ltwb0962/.conda/envs/bcftools_env/bin:$PATH"

set -eo pipefail

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" fastq/samples_500.txt)

echo "Processing sample: ${SAMPLE}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "BWA: $(bwa 2>&1 | grep Version)"
echo "Samtools: $(samtools --version | head -1)"

bwa mem \
    -t 10 \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA\tLB:1\tPU:ILLUMINA" \
    TAIR10.fa \
    "fastq/results/trimmed/${SAMPLE}.R1.trimmed.fastq.gz" \
    "fastq/results/trimmed/${SAMPLE}.R2.trimmed.fastq.gz" \
    | \
samtools sort \
    -@ 10 \
    -o "results/mapped/${SAMPLE}.sorted.bam" \
    -T "results/mapped/${SAMPLE}.tmp"


#!/bin/bash
#SBATCH --job-name=fastqc500
#SBATCH --array=1-500
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --output=fastqc_logs/%A_%a.out
#SBATCH --error=fastqc_logs/%A_%a.err

source ~/.bashrc
conda activate qc_env

set -eo pipefail

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples_500.txt)

echo "Processing sample: ${SAMPLE}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"

fastqc \
    --threads 2 \
    --outdir results/trimmed_fastqc \
    "results/trimmed/${SAMPLE}_val_1.fq.gz" \
    "results/trimmed/${SAMPLE}_val_2.fq.gz"

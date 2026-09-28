#!/bin/bash
#SBATCH --job-name=trim500
#SBATCH --array=1-500
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=2:00:00
#SBATCH --output=trim_logs/%A_%a.out
#SBATCH --error=trim_logs/%A_%a.err

source ~/.bashrc
conda activate trim_reads

set -eo pipefail

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples_500.txt)

echo "Processing sample: ${SAMPLE}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"

trim_galore \
    --cores 4 \
    --gzip \
    --paired \
    --output_dir results/trimmed \
    --basename "${SAMPLE}" \
    "${SAMPLE}_1.fastq.gz" \
    "${SAMPLE}_2.fastq.gz"

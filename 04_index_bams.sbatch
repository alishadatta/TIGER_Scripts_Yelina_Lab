#!/bin/bash
#SBATCH --job-name=index500
#SBATCH --output=logs/index500_%A_%a.log
#SBATCH --error=logs/index500_%A_%a.err
#SBATCH --array=1-500
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --chdir=/scratch/ltwb0962/Tiger_Pipeline

set -euo pipefail

BASE=/scratch/ltwb0962/Tiger_Pipeline
SAMTOOLS=/scratch/leedsomics_tools/samtools-1.23/samtools
BAM_DIR=$BASE/results/mapped
SAMPLE_LIST=$BASE/input/sample_list_500.txt

BAM_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
BAM="$BAM_DIR/$BAM_NAME"

echo "Indexing: $BAM"

"$SAMTOOLS" index "$BAM"

echo "Finished: $BAM"

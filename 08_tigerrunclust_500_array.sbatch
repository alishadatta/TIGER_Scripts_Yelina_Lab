#!/bin/bash
#SBATCH --job-name=TIGER500
#SBATCH --array=1-500
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=/scratch/ltwb0962/Tiger_Pipeline/logs/tiger_run_%A_%a.log
#SBATCH --error=/scratch/ltwb0962/Tiger_Pipeline/logs/tiger_run_%A_%a.err

cd /scratch/ltwb0962/Tiger_Pipeline

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" input/sample_list_500.txt)
SAMPLE=${SAMPLE%.sorted.bam}

echo "=========================================="
echo "SLURM task: $SLURM_ARRAY_TASK_ID"
echo "Sample: $SAMPLE"
echo "Start: $(date)"
echo "=========================================="

bash TIGER_Module/tiger-pipeline/workflow/TIGER_Scripts-for-distribution/tigerrunclust_AD.sh \
    /scratch/ltwb0962/Tiger_Pipeline/TIGER_Module/tiger-pipeline/workflow/TIGER_Scripts-for-distribution \
    /scratch/ltwb0962/Tiger_Pipeline/input/tiger500 \
    "$SAMPLE"

STATUS=$?

if [[ $STATUS -eq 0 ]]; then
    echo "Finished TIGER successfully: $SAMPLE"

elif [[ -s "input/tiger500/$SAMPLE/rough_COs_${SAMPLE}.breaks.txt" ]] && \
     [[ -s "input/tiger500/$SAMPLE/rough_COs_${SAMPLE}.refined.breaks.txt" ]] && \
     [[ -s "input/tiger500/$SAMPLE/${SAMPLE}_refined.corrected.breaks" ]]; then

    echo "TIGER analysis completed for $SAMPLE, but plotting failed."
    echo "Core TIGER output files are present."
    echo "Finished sample: $SAMPLE"

else
    echo "FAILED sample: $SAMPLE (exit code $STATUS)"
    exit $STATUS
fi

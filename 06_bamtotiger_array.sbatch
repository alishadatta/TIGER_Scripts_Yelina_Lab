#!/bin/bash
#SBATCH --job-name=tiger500
#SBATCH --output=logs/tiger500_%A_%a.log
#SBATCH --error=logs/tiger500_%A_%a.err
#SBATCH --array=1-500
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --chdir=/scratch/ltwb0962/Tiger_Pipeline

set -euo pipefail

BASE=/scratch/ltwb0962/Tiger_Pipeline

SAMTOOLS=$HOME/.conda/envs/bcftools_env/bin/samtools
BCFTOOLS=$HOME/.conda/envs/bcftools_env/bin/bcftools

REF=$BASE/TAIR10.fa
SAMPLE_LIST=$BASE/input/sample_list_500.txt
BAM_DIR=$BASE/results/mapped

COMP_MARKERS=$BASE/markers/complete_marker.txt
CORR_MARKERS=$BASE/markers/corrected_marker_provisional.txt
COMP_REGIONS=$BASE/input/tiger500/complete_regions.txt

OUT_DIR=$BASE/input/tiger500

BAM_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST")
BAM=$BAM_DIR/$BAM_NAME
SAMPLE=${BAM_NAME%.sorted.bam}

echo "=========================================="
echo "Sample: $SAMPLE"
echo "BAM:    $BAM"
echo "=========================================="

if [[ ! -s "$BAM" ]]; then
    echo "ERROR: BAM not found: $BAM"
    exit 1
fi

VCF=$OUT_DIR/${SAMPLE}.vcf.gz
CALLS=$OUT_DIR/${SAMPLE}.calls.txt
COMP_OUT=$OUT_DIR/${SAMPLE}.complete.txt
CORR_OUT=$OUT_DIR/${SAMPLE}.corrected.txt

echo "Running samtools mpileup..."

"$SAMTOOLS" mpileup \
    -E \
    --skip-indels \
    --output-tags DP,DP4 \
    -f "$REF" \
    -R "$COMP_REGIONS" \
    "$BAM" \
    -Ou | \
"$BCFTOOLS" call \
    --multiallelic-caller \
    --skip-variants indels \
    -Oz \
    -o "$VCF"

"$BCFTOOLS" index -f "$VCF"

echo "Extracting allele counts..."

"$BCFTOOLS" query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%DP4\n' \
    "$VCF" | \
awk '
BEGIN { OFS="\t" }
{
    split($5,d,",")
    ref_reads=d[1]+d[2]
    alt_reads=d[3]+d[4]

    if (substr($1,1,3) == "Chr")
        chr=substr($1,4)
    else
        chr=$1

    print chr,$2,$3,ref_reads,$4,alt_reads
}
' > "$CALLS"

echo "Creating complete TIGER input..."

perl "$BASE/TIGER_Module/tiger-pipeline/workflow/TIGER_Scripts-for-distribution/get_subset_2.pl" \
    -c "$CALLS" \
    -p "$COMP_MARKERS" \
    -o "$COMP_OUT"

echo "Creating corrected TIGER input..."

perl "$BASE/TIGER_Module/tiger-pipeline/workflow/TIGER_Scripts-for-distribution/get_subset_2.pl" \
    -c "$CALLS" \
    -p "$CORR_MARKERS" \
    -o "$CORR_OUT"

echo "=========================================="
echo "Finished: $SAMPLE"
echo "=========================================="

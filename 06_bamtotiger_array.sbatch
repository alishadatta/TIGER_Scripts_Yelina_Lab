#!/bin/bash
#SBATCH --job-name=tiger500
#SBATCH --output=logs/tiger500_%A_%a.log
#SBATCH --error=logs/tiger500_%A_%a.err
#SBATCH --array=1-500
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --chdir=/scratch/ltwb0962/Tiger_Pipeline

set -euo pipefail

BASE=/scratch/ltwb0962/Tiger_Pipeline

BCFTOOLS=$HOME/.conda/envs/bcftools_env/bin/bcftools

REF=$BASE/TAIR10.fa
SAMPLE_LIST=$BASE/input/sample_list_500.txt
BAM_DIR=$BASE/results/mapped

COMP_ALLELES=$BASE/markers/complete_marker_alleles.txt
CORR_ALLELES=$BASE/markers/corrected_marker_provisional_alleles.txt

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

if [[ ! -s "$COMP_ALLELES" ]]; then
    echo "ERROR: Complete allele file not found"
    exit 1
fi

if [[ ! -s "$CORR_ALLELES" ]]; then
    echo "ERROR: Corrected allele file not found"
    exit 1
fi

TMP_COUNTS=$(mktemp)
TMP_COMPLETE=$(mktemp)
TMP_CORRECTED=$(mktemp)

trap 'rm -f "$TMP_COUNTS" "$TMP_COMPLETE" "$TMP_CORRECTED"' EXIT

echo "Running bcftools mpileup + call..."

"$BCFTOOLS" mpileup \
    -f "$REF" \
    -R "$COMP_REGIONS" \
    "$BAM" \
    -Ou | \
"$BCFTOOLS" call \
    -m \
    --skip-variants indels \
    -Ou | \
"$BCFTOOLS" query \
    -f '%CHROM\t%POS\t%DP4\n' | \
awk '
BEGIN { OFS="\t" }
{
    split($3,d,",")
    ref_reads=d[1]+d[2]
    alt_reads=d[3]+d[4]

    print substr($1,4),$2,ref_reads,alt_reads
}
' > "$TMP_COUNTS"

echo "Building complete TIGER input..."

awk '
BEGIN { OFS="\t" }

NR==FNR {
    allele[$1"\t"$2]=$3"\t"$4
    next
}

{
    key=$1"\t"$2
    count[key]=$3"\t"$4
}

END {
    for (key in allele) {
        split(key,k,"\t")
        split(allele[key],a,"\t")

        if (key in count) {
            split(count[key],c,"\t")
            print k[1],k[2],a[1],c[1],a[2],c[2]
        }
        else {
            print k[1],k[2],a[1],0,a[2],0
        }
    }
}
' "$COMP_ALLELES" "$TMP_COUNTS" |
sort -k1,1n -k2,2n > "$TMP_COMPLETE"

echo "Building corrected TIGER input..."

awk '
BEGIN { OFS="\t" }

NR==FNR {
    allele[$1"\t"$2]=$3"\t"$4
    next
}

{
    key=$1"\t"$2
    count[key]=$3"\t"$4
}

END {
    for (key in allele) {
        split(key,k,"\t")
        split(allele[key],a,"\t")

        if (key in count) {
            split(count[key],c,"\t")
            print k[1],k[2],a[1],c[1],a[2],c[2]
        }
        else {
            print k[1],k[2],a[1],0,a[2],0
        }
    }
}
' "$CORR_ALLELES" "$TMP_COUNTS" |
sort -k1,1n -k2,2n > "$TMP_CORRECTED"

mv "$TMP_COMPLETE" "$OUT_DIR/input_complete_${SAMPLE}.txt"
mv "$TMP_CORRECTED" "$OUT_DIR/input_corrected_${SAMPLE}.txt"

echo "Checking output..."

COMP_N=$(wc -l < "$OUT_DIR/input_complete_${SAMPLE}.txt")
CORR_N=$(wc -l < "$OUT_DIR/input_corrected_${SAMPLE}.txt")

echo "Complete markers:  $COMP_N"
echo "Corrected markers: $CORR_N"

if [[ "$COMP_N" -ne 484712 ]]; then
    echo "ERROR: Unexpected complete marker count!"
    exit 1
fi

if [[ "$CORR_N" -ne 200998 ]]; then
    echo "ERROR: Unexpected corrected marker count!"
    exit 1
fi

echo "=========================================="
echo "Finished: $SAMPLE"
echo "=========================================="

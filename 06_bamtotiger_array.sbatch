#!/bin/bash
#SBATCH --job-name=tiger_beth
#SBATCH --output=logs/tiger_beth_%A_%a.log
#SBATCH --error=logs/tiger_beth_%A_%a.err
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

BETH_MASK=$BASE/markers/beth.mask.txt
BETH_REGIONS=$BASE/markers/beth.mask.regions.txt

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

if [[ ! -s "$BETH_MASK" ]]; then
    echo "ERROR: beth.mask.txt not found"
    exit 1
fi

if [[ ! -s "$BETH_REGIONS" ]]; then
    echo "ERROR: beth.mask.regions.txt not found"
    exit 1
fi

TMP_COUNTS=$(mktemp)

trap 'rm -f "$TMP_COUNTS"' EXIT

echo "Running bcftools mpileup + call..."

"$BCFTOOLS" mpileup \
    -f "$REF" \
    -R "$BETH_REGIONS" \
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

echo "Building Beth TIGER input..."

awk '
BEGIN { OFS="\t" }

NR==FNR {
    count[$1"\t"$2]=$3"\t"$4
    next
}

{
    key=$1"\t"$2

    if (key in count) {
        split(count[key],c,"\t")
        print $1,$2,$3,c[1],$5,c[2]
    } else {
        print $1,$2,$3,0,$5,0
    }
}
' "$TMP_COUNTS" "$BETH_MASK" \
| sort -k1,1n -k2,2n \
> "$OUT_DIR/input_beth_${SAMPLE}.txt"

echo "Checking output..."

N=$(wc -l < "$OUT_DIR/input_beth_${SAMPLE}.txt")

echo "Beth markers: $N"

if [[ "$N" -ne 132892 ]]; then
    echo "ERROR: Unexpected marker count!"
    exit 1
fi

echo "=========================================="
echo "Finished: $SAMPLE"
echo "=========================================="

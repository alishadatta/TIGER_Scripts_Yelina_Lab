#!/bin/bash
# tigerrunclust.sh
# Runs TIGER pipeline for one sample

set -euo pipefail

# --------------------------------------------------
# Arguments
# --------------------------------------------------

TIGER=$1
DATA=$2
SAMPLE=$3

# --------------------------------------------------
# Paths
# --------------------------------------------------

SAMPLE_DIR=$DATA/$SAMPLE

mkdir -p "$SAMPLE_DIR"

echo "=========================================="
echo "Running TIGER for sample: $SAMPLE"
echo "TIGER directory: $TIGER"
echo "Data directory:  $DATA"
echo "Output directory: $SAMPLE_DIR"
echo "=========================================="

# --------------------------------------------------
# Java memory
# --------------------------------------------------

export _JAVA_OPTIONS=-Xmx1024M

# --------------------------------------------------
# TIGER input files
# --------------------------------------------------

CORRECTED=$DATA/input_corrected_${SAMPLE}.txt
COMPLETE=$DATA/input_complete_${SAMPLE}.txt

if [[ ! -s "$CORRECTED" ]]; then
    echo "ERROR: corrected input not found:"
    echo "$CORRECTED"
    exit 1
fi

if [[ ! -s "$COMPLETE" ]]; then
    echo "ERROR: complete input not found:"
    echo "$COMPLETE"
    exit 1
fi

# --------------------------------------------------
# Base caller
# --------------------------------------------------

echo "[$SAMPLE] Running base caller..."

cd "$TIGER"

java -jar base_caller.jar \
    -r "$CORRECTED" \
    -o "$SAMPLE_DIR/allele_count_base_call_${SAMPLE}.txt" \
    -n bi

# --------------------------------------------------
# Allele frequency estimator
# --------------------------------------------------

echo "[$SAMPLE] Running allele frequency estimator..."

java -jar allele_freq_estimator.jar \
    -r "$CORRECTED" \
    -o "$SAMPLE_DIR/frequencies_for_bmm_${SAMPLE}.txt" \
    -n bi \
    -w 1000

# --------------------------------------------------
# Beta mixture model
# --------------------------------------------------

echo "[$SAMPLE] Running beta mixture model..."

R --slave --vanilla --args \
    "$SAMPLE_DIR/frequencies_for_bmm_${SAMPLE}.txt" \
    "$SAMPLE_DIR/bmm.intersections_${SAMPLE}.txt" \
    < "$TIGER/beta_mixture_model.R"

# --------------------------------------------------
# Prepare HMM probability input
# --------------------------------------------------

echo "[$SAMPLE] Preparing HMM probability input..."

perl "$TIGER/prep_prob.pl" \
    -s "$SAMPLE" \
    -m "$CORRECTED" \
    -b "$SAMPLE_DIR/allele_count_base_call_${SAMPLE}.txt" \
    -c "$TIGER/TAIR10_chrSize.txt" \
    -o "$SAMPLE_DIR/file_for_probabilities_${SAMPLE}.txt"

# --------------------------------------------------
# Calculate HMM probabilities
# --------------------------------------------------

echo "[$SAMPLE] Calculating HMM probabilities..."

perl "$TIGER/hmm_prob.pl" \
    -s "$SAMPLE_DIR/frequencies_for_bmm_${SAMPLE}.txt" \
    -p "$SAMPLE_DIR/file_for_probabilities_${SAMPLE}.txt" \
    -o "$SAMPLE_DIR/INDEX${SAMPLE}" \
    -a "$SAMPLE_DIR/bmm.intersections_${SAMPLE}.txt" \
    -c "$TIGER/TAIR10_chrSize.txt"

# --------------------------------------------------
# Run HMM
# --------------------------------------------------

echo "[$SAMPLE] Running HMM..."

java -jar "$TIGER/hmm_play.jar" \
    -r "$SAMPLE_DIR/allele_count_base_call_${SAMPLE}.txt" \
    -o "$SAMPLE_DIR/hmm.out_${SAMPLE}.txt" \
    -t bi \
    -z "$SAMPLE_DIR/INDEX${SAMPLE}_hmm_model"

# --------------------------------------------------
# Rough crossover positions
# --------------------------------------------------

echo "[$SAMPLE] Finding rough crossover positions..."

perl "$TIGER/prepare_break.pl" \
    -s "$SAMPLE" \
    -m "$CORRECTED" \
    -b "$SAMPLE_DIR/hmm.out_${SAMPLE}.txt" \
    -c "$TIGER/TAIR10_chrSize.txt" \
    -o "$SAMPLE_DIR/rough_COs_${SAMPLE}.txt"

# --------------------------------------------------
# Refine crossover positions using complete markers
# --------------------------------------------------

echo "[$SAMPLE] Refining crossover positions..."

perl "$TIGER/refine_recombination_break.pl" \
    "$COMPLETE" \
    "$SAMPLE_DIR/rough_COs_${SAMPLE}.breaks.txt"

# --------------------------------------------------
# Smooth crossover positions
# --------------------------------------------------

echo "[$SAMPLE] Smoothing crossover positions..."

perl "$TIGER/breaks_smoother.pl" \
    -b "$SAMPLE_DIR/rough_COs_${SAMPLE}.refined.breaks.txt" \
    -o "$SAMPLE_DIR/${SAMPLE}_refined.corrected.breaks"

# --------------------------------------------------
# Plot / visualise
# --------------------------------------------------

echo "[$SAMPLE] Generating plot..."

cd "$SAMPLE_DIR"

R --slave --vanilla --args \
    "$SAMPLE" \
    "$SAMPLE_DIR/sample_${SAMPLE}.visual_out.pdf" \
    "$SAMPLE_DIR/rough_COs_${SAMPLE}.breaks.txt" \
    "$SAMPLE_DIR/rough_COs_${SAMPLE}.refined.breaks.txt" \
    "$SAMPLE_DIR/${SAMPLE}_refined.corrected.breaks" \
    "$SAMPLE_DIR/frequencies_for_bmm_${SAMPLE}.txt" \
    "$SAMPLE_DIR/INDEX${SAMPLE}_sliding_window.breaks.txt" \
    < "$TIGER/plot_genotyping.R"

echo "=========================================="
echo "Finished TIGER for sample: $SAMPLE"
echo "=========================================="

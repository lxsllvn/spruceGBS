#!/bin/bash

# This script implements a parameter sweep pipeline for testing ANGSD settings
# across a ~100 Mbp filtered region using experimental samples.

# ------------------------
# REQUIRED INPUTS
# ------------------------
# 1. REF: Path to reference genome
# 2. DOMAIN: domain name (e.g., southern, siberia)
# 3. BAMLIST: list of BAM files for the entire domain (e.g., southern_bams.txt)
# 4. REGION: list of scaffolds for ~1 Mbp filtered region
# 5. SITES: sites to analyze, tab-delimited 1-indexed format for -sites
# 6. OUTDIR: base output directory
# 7. Per-population BAM lists should be stored in "${DOMAIN}_populations/"

# ------------------------
# OUTPUTS
# ------------------------
# 1. counts: domain-level read count matrix (*.counts.gz)
# 2. *.gl.mafs.gz: per-site minor allele frequencies
# 3. *.gl.hwe.gz: F and HWE-based statistics (used to derive Hobs/Hexp)
# 4. *.gl.beagle.gz: genotype likelihoods in BEAGLE format
# 5. thetas.persite.txt: per-site, logscale π and other theta estimators
# 6. thetas.summary.pestPG: per-scaffold summary of π and other theta estimators

# ------------------------
# USER-DEFINED PARAMETERS
# ------------------------

BAQ_LIST=(0 1 2)
C_LIST=(0 50)
MINQ_LIST=(20 30)
MINMAPQ_LIST=(20 30 40)

DEFAULT_FILTERS="-remove_bads 1 -uniqueOnly 1 -only_proper_pairs 1 -skipTriallelic 1 -setMinDepthInd 3"

# Load required modules 
ml GCC/10.2.0
ml angsd/0.935

# ------------------------
# SCRIPT ARGUMENTS
# ------------------------

REF="$1"
DOMAIN="$2"
BAMLIST="$3"
REGION="$4"
SITES="$5"
OUTDIR="$6"

mkdir -p "$OUTDIR"

for baq in "${BAQ_LIST[@]}"; do
  for C in "${C_LIST[@]}"; do
    for minQ in "${MINQ_LIST[@]}"; do
      for minMapQ in "${MINMAPQ_LIST[@]}"; do

        PARAM_ID="baq${baq}_C${C}_q${minQ}_mq${minMapQ}"
        PREFIX_BASE="${OUTDIR}/${DOMAIN}_${PARAM_ID}"
        mkdir -p "$PREFIX_BASE"

        echo "Running ANGSD parameter set: $PARAM_ID"

        # Step 1: Calculate domain-wide read count matrix and genotype likelihoods
        angsd -bam "$BAMLIST" \
          -ref "$REF" \
          -anc "$REF" \
          -rf "$REGION" \
          -sites "$SITES" \
          -GL 1 \
          -doMajorMinor 1 \
          -doMaf 1 \
          -doCounts 1 \
          -dumpCounts 2 \
          -doGlf 2 \
          -baq $baq \
          -C $C \
          -minQ $minQ \
          -minMapQ $minMapQ \
          $DEFAULT_FILTERS \
          -out "$PREFIX_BASE/domain_sfs"

        # Step 2: Run per-population analyses
        for POP_BAMLIST in ${DOMAIN}_populations/*.txt; do
          POP=$(basename "$POP_BAMLIST" .txt)
          POP_PREFIX="$PREFIX_BASE/$POP"
          mkdir -p "$POP_PREFIX"

          angsd -bam "$POP_BAMLIST" \
            -ref "$REF" \
            -anc "$REF" \
            -rf "$REGION" \
            -sites "$SITES" \
            -GL 1 \
            -doSaf 1 \
            -doMajorMinor 1 \
            -doMaf 1 \
            -doHWE 1 \
            -doCounts 1 \
            -doGlf 2 \
            -baq $baq \
            -C $C \
            -minQ $minQ \
            -minMapQ $minMapQ \
            $DEFAULT_FILTERS \
            -out "$POP_PREFIX/gl"
            
          # Estimate site frequency spectrum (SFS) from "sample allele frequencies" (SAF)
          # note: SAFs are the per-site likelihood surfaces over allele counts, not frequencies.
          realSFS "$POP_PREFIX/gl.saf.idx" -fold 1 > "$POP_PREFIX/gl.sfs" 
        
          # Calculate per-site theta using the SFS as a prior 
          realSFS saf2theta "$POP_PREFIX/gl.saf.idx" \
            -sfs "$POP_PREFIX/gl.sfs" \
            -fold 1 \
            -outname "$POP_PREFIX/saf2theta"
        
          # Estimate Tajimas D and other statistics estimates per window
          # Produces "theta.summary.pestPG"
          thetaStat do_stat "$POP_PREFIX/saf2theta.thetas.idx" \
            -outnames "$POP_PREFIX/thetas.summary"
        
          # Extract the logscale per site thetas estimates
          thetaStat print "$POP_PREFIX/saf2theta.thetas.idx" \
            > "$POP_PREFIX/thetas.persite.txt"
        done
      done
    done
  done
done

echo "Parameter sweep complete for $DOMAIN."

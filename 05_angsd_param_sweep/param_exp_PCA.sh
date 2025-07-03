#!/bin/bash

ml GCC/13.2.0
ml SciPy-bundle/2023.11

INPUT="$1"
OUTPUT="$2"

pcangsd \
  -b "${INPUT}.beagle.gz" \
  -o "${OUTPUT}.Pcangsd" \
  --sites-save --snp-weights --pcadapt --selection --iter 500 --maf-iter 1000

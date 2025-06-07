#!/bin/bash

ml GCC/13.2.0
ml SciPy-bundle/2023.11

INPUT="$1"
OUTPUT="$2"

pcangsd \
  -b "${INPUT}" \
  -o "${OUTPUT}" \
  --sites_save --snp_weights --pcadapt --selection --iter 500 --maf_iter 1000

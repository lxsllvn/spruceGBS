#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J pcangsd
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --constraint=skylake
#SBATCH -t 0-01:00:00

ml GCC/13.2.0
ml SciPy-bundle/2023.11

INPUT="$1"
OUTPUT="$2"
/home/l/lxsllvn/Public/python-modules/bin/bin/pcangsd \
  -b "${INPUT}.beagle.gz" \
  -o "${OUTPUT}.Pcangsd" \
  --sites_save --snp_weights --pcadapt --selection --iter 500 --maf_iter 1000

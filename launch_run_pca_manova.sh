#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J R_pca_summary
#SBATCH --output=r_pca.out
#SBATCH --error=r_pca.err
#SBATCH --constraint=skylake
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 0-01:00:00

ml GCC/10.2.0  CUDA/11.1.1  OpenMPI/4.0.5
ml R/4.0.4

Rscript  ${SCRIPTS}/run_pca_manova.R \
  siberia \
  siberia_results \
  parameter_exp_siberia_bamlist \
  sequenced_samples_metadata.csv

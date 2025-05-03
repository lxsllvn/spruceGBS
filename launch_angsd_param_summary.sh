#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J R_param_summary
#SBATCH --output=r_param_summary.out
#SBATCH --error=r_param_summary.err
#SBATCH --constraint=skylake
#SBATCH -n 1
#SBATCH -c 3
#SBATCH -t 0-01:00:00

ml GCC/10.2.0  CUDA/11.1.1  OpenMPI/4.0.5
ml R/4.0.4

# Rscript ${SCRIPTS}/angsd_param_exp_summary.R results_folder meta_data.csv domain_bamlist out_name
 
#Rscript ${SCRIPTS}/angsd_param_exp_summary.R southern_results sequenced_samples_metadata.csv parameter_exp_southern_bamlist southern



Rscript ${SCRIPTS}/summary_with_beagle.R \
  siberia_results \
  sequenced_samples_metadata.csv \
  parameter_exp_siberia_bamlist \
  siberia

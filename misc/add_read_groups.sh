#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J addreplaceRG
#SBATCH -n 1
#SBATCH -c 6
#SBATCH -t 0-02:00:00

# Load required modules
ml GCC/13.2.0               
ml SAMtools/1.19.2        

# Get the sample name from the first command-line argument
SAMPLE="$1"

samtools addreplacerg -r "ID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" -o "${SAMPLE}.intersect.sorted.RG.bam" "${SAMPLE}.intersect.sorted.bam"

samtools index "${SAMPLE}.intersect.sorted.RG.bam"

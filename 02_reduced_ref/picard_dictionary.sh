#!/bin/bash
#SBATCH -A hpc2n2024-161
#SBATCH -J seq_dict
#SBATCH --output=seq_dict.out
#SBATCH --error=seq_dict.err
#SBATCH -n 1
#SBATCH -c 6
#SBATCH -t 1-00:00:00

# Load required modules
ml Java/17.0.6

# Path to picard.jar
PICARDPATH="/pfs/proj/nobackup/fs/projnb10/wanglab-2020/Alex/alignments"

# Path to reference genome
REF="${SPRUCE_PROJECT}/ref/picea_newref.fa"

# Path for output dictionary
OUTPUT="${SPRUCE_PROJECT}/ref/picea_newref.dict"

# Create sequence dictionary for indel realignment with GATK
java -jar "${PICARDPATH}/picard.jar" \
    CreateSequenceDictionary \
    R="$REF" \
    O="$OUTPUT"
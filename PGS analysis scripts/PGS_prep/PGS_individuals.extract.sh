#!/bin/sh
#SBATCH --job-name=GRM_individuals
#SBATCH --time=00:05:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10000


######################################################################
# Create a relatedness file using GCTA.
# https://yanglab.westlake.edu.cn/software/gcta/#GREMLanalysis

# Load GCTA tool in Linux bash
PATH=/groups/umcg-lifelines/tmp02/projects/ov21_0398/Software/gcta-1.94.3-linux-kernel-3-x86_64/:$PATH
# Activate it with this command: gcta64

#---------------------------------------
# Defining directories and paths and load modules
#---------------------------------------
WorkingDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_prep"
GRMDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/GRMs/GRMs_ugli2_gsa_cyto/fam.and.unrel"

#---------------------------------------
# Use GCTA to create a relatedness file
#---------------------------------------
gcta64 \
--grm $GRMDir/GRM.All \
--grm-cutoff 0.125 \
--make-grm \
--out $WorkingDir/output/relatedness.0.125_individuals







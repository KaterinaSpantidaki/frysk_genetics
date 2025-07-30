#!/bin/sh
#SBATCH --job-name=GCTA_snph2_per_pheno
#SBATCH --time=14:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=150000

######################################################################
# Calculate SNP heritability using GCTA.
# https://yanglab.westlake.edu.cn/software/gcta/#GREMLanalysis

# Load GCTA tool in Linux bash
PATH=/groups/umcg-lifelines/tmp02/projects/ov21_0398/Software/gcta-1.94.3-linux-kernel-3-x86_64/:$PATH
# Activate it with this command: gcta64

#---------------------------------------
# Defining directories and paths and load modules
#---------------------------------------

WorkingDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk"

# Binary heritability
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# for ugli 1 and 2

gcta64 \
--mgrm $WorkingDir/GRMs/GRMs_ugli2_gsa_cyto/bK/GRM.bK.02.mgrm.txt \
--pheno $WorkingDir/heritability/input_files/pheno_ugli_cyto_foreign_language_learning_cogn_filtered.txt \
--covar $WorkingDir/heritability/input_files/covar_binary_ugli1and2_cyto_cogn_filtered.txt \
--qcovar $WorkingDir/heritability/input_files/covar_ugli_cyto_quant_ugli1and2_accents_immitation_learning_cogn_filtered.txt \
--reml-no-constrain \
--out $WorkingDir/heritability/output_files/SNPh2.bK02.UGLI1and2.CYTO.foreign_language_learning_cogn_filtered \
--thread-num 10




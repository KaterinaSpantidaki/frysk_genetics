#!/bin/sh
#SBATCH --job-name=GREML_gencor
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=200000

#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 21 May 2024.
# title: "Calculate genetic correlation using GCTA".
# date: 10 April 2025.
# data: Lifelines 
#---
# Script to calculate genetic correlation using GCTA for language learning phenotypic data: n_fluent_languages, foreign_accents_understanding, foreign_speech_imitation, foreign_language_learning.
# https://yanglab.westlake.edu.cn/software/gcta/#BivariateGREMLanalysis
#---

#---------------------------------------
# Defining directories and paths and load modules.
#---------------------------------------

# Load GCTA tool in Linux bash
PATH=/groups/umcg-lifelines/tmp02/projects/ov21_0398/Software/gcta-1.94.3-linux-kernel-3-x86_64/:$PATH
# Activate it with this command: gcta64

# Set working directory.
WorkingDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/"


#---------------------------------------
# Calculate genetic correlations
#---------------------------------------
gcta64 \
 --reml-bivar \
 --mgrm $WorkingDir/GRMs/GRMs_ugli2_gsa_cyto/bK/GRM.bK.02.mgrm.txt \
 --pheno $WorkingDir/gen_cor/input_files/pheno_ugli_cyto_bivar_cogn_filtered_UI.txt \
 --covar $WorkingDir/gen_cor/input_files/covar_binary_ugli1and2_cyto_bivar_cogn_filtered_UI.txt \
 --qcovar $WorkingDir/gen_cor/input_files/covar_quant_ugli1and2_cyto_bivar_cogn_filtered_UI.txt \
 --thread-num 10 \
 --out $WorkingDir/gen_cor/output_files/GREML_bivar.bk02_cogn_filtered_UI
 


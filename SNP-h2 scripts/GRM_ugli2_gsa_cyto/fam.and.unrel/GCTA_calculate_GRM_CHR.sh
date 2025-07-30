#!/bin/sh
#SBATCH --job-name=__CHR___GCTA_GRM
#SBATCH --time=3:30:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=200000

#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 23 April 2024.
# title: "GRM construction from plink files using GCTA - GRM per chromosome".
# date: 14 April 2025.
# data: Lifelines 
#---
# Script for GRM generation per chromosome (to be submitted).
# Run this script to run all scripts and create all GRMs: /groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/GCTA/GRM_ugli2_gsa_cyto/fam.and.unrel/GCTA_calculate_GRM_make_and_run_all_chr.sh
#---

#---------------------------------------
# Defining directories and paths and load modules.
#---------------------------------------

# Set working directory.
WorkingDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/GRMs"

# Load GCTA tool in Linux bash.
PATH=/groups/umcg-lifelines/tmp02/projects/ov21_0398/Software/gcta-1.94.3-linux-kernel-3-x86_64/:$PATH
# Activate it with this command: gcta64

#---------------------------------------
# Calculate GRM for each chromosome for total dataset (UGLI2, GSA and CYTO combined).
#---------------------------------------

gcta64 --bfile $WorkingDir/GRMs_ugli2_gsa_cyto/All_intermediate_files/__CHR__.All.filtered \
       --chr __CHR__ \
       --make-grm \
       --out $WorkingDir/GRMs_ugli2_gsa_cyto/All_intermediate_files/GRM.__CHR__.All.filtered \
       --thread-num 10

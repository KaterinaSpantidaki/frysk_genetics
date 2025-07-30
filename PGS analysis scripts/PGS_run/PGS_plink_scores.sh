#!/bin/sh
#SBATCH --job-name=PGS_plink
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=100000

#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 11 November 2024.
# title: "Obtain individual-level polygenic scores by concatenating output files from all chromosomes with plink".
# date: 15 November 2025.
# data: Lifelines 
#---

#---------------------------------------
# Defining directories and paths and load modules
#---------------------------------------
module load PLINK/2.0-alpha4.5-20230813
module load BCFtools/1.17-GCCcore-11.3.0
module load PLINK/1.9-beta6-20190617

WorkingDir=/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_run
GenoDir=/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_prep/output_com

#---------------------------------------
# Run PLINK --scores to obtain individual PGS scores.
#---------------------------------------
predictor=$1
outcome=$2

plink --bfile $GenoDir/"${outcome}"_PRScs_geno \
--score $WorkingDir/output/"${predictor}"_"${outcome}"_pst_eff_a1_b0.5_phiauto_all.txt 2 4 6 \
--out $WorkingDir/output_pgs_plink/"${predictor}"_"${outcome}"



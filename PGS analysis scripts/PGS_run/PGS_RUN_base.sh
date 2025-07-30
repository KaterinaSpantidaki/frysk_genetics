#!/bin/sh
#SBATCH --job-name=calculate_pgs
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=100000

#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 11 November 2024.
# title: "Calculate PGS".
# date: 20 May 2025.
# data: Lifelines 
#---

#---------------------------------------
# Defining directories and paths and load modules
#---------------------------------------
module load Anaconda3/2024.02-1
module load PLINK/1.9-beta6-20190617   

PRScsDir=/groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/PGS/PGS_run/PRScs
RefDir=/groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/PGS/PGS_run/lang_soc
GenoDir=/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_prep/output_com
SumstatDir=/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_prep/sumstats
OutDir=/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_run/output

cd $PRScsDir

#---------------------------------------
# Run PRScs
#---------------------------------------
predictor=$1
outcome=$2
n_gwas=$3

python PRScs.py \
--ref_dir=$RefDir/ldblk_1kg_eur \
--bim_prefix=$GenoDir/"${outcome}_PRScs_geno" \
--sst_file=$SumstatDir/"${predictor}"_PRScs_filtered_flipped.txt \
--n_gwas=$n_gwas \
--out_dir=$OutDir/"${predictor}"_"${outcome}"

wait

#---------------------------------------
# PRScs outputs one file per chromosome. Concatenate all these files into a single file.
#---------------------------------------
cat $OutDir/"${predictor}"_"${outcome}"_pst_eff_a1_b0.5_phiauto_chr*.txt > $OutDir/"${predictor}"_"${outcome}"_pst_eff_a1_b0.5_phiauto_all.txt






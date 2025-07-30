#!/bin/sh
#SBATCH --job-name=combine_chrom
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=100000

#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 11 November 2024.
# title: "Prep of genotype input files for PGS - combine genotype chromosome files".
# date: 14 May 2025.
# data: Lifelines 
#---

#---------------------------------------
# Defining directories and paths and load modules
#---------------------------------------
WorkingDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_prep"

module load PLINK/2.0-alpha4.5-20230813
module load BCFtools/1.17-GCCcore-11.3.0
module load PLINK/1.9-beta6-20190617


# Create a list of plink files for all chromosomes
ls $WorkingDir/output_PLINK/All.chr*.filtered.bed > $WorkingDir/output_com/All.chr.filtered_plink.files_list.txt
sed -i 's|.bed||g' $WorkingDir/output_com/All.chr.filtered_plink.files_list.txt
chmod 700 $WorkingDir/output_com/All.chr.filtered_plink.files_list.txt

# Merge files with plink
plink --merge-list $WorkingDir/output_com/All.chr.filtered_plink.files_list.txt \
  --make-bed \
  --out $WorkingDir/output_com/All.chr.filtered



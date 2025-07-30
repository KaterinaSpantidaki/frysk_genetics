#!/bin/sh
#SBATCH --job-name=plink_per.trait
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=100000

#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 11 November 2024.
# title: "Prep of genotype input files for PGS - obtain plink files per trait".
# date: 20 May 2025.
# data: Lifelines 
#---

#---------------------------------------
# Defining directories and paths and load modules
#---------------------------------------
module load PLINK/2.0-alpha4.5-20230813
module load BCFtools/1.17-GCCcore-11.3.0
module load PLINK/1.9-beta6-20190617

WorkingDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_prep/output_com"

#---------------------------------------
# Generate plink file per trait
#---------------------------------------
trait=$1

plink \
--bfile $WorkingDir/All.chr.filtered \
--keep $WorkingDir/"${trait}_snp_list.txt" \
--make-bed \
--out $WorkingDir/"${trait}_PRScs_geno"



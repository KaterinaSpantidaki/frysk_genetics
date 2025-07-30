#!/bin/sh
#SBATCH --job-name=chr17_PGS_prep_geno_input_filter.PLINK
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=100000

#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 2 November 2024.
# title: "Prep of genotype input files for PGS - filter with PLINK".
# date: 16 May 2025.
# data: Lifelines 
#---
# QC steps followed from following link: https://choishingwan.github.io/PRS-Tutorial/target/


#---------------------------------------
# Defining directories and paths and load modules
#---------------------------------------
module load PLINK/2.0-alpha4.5-20230813
module load BCFtools/1.17-GCCcore-11.3.0
module load PLINK/1.9-beta6-20190617

WorkingDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_prep"
GRMDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/GRMs/GRMs_ugli2_gsa_cyto/All_intermediate_files"


#---------------------------------------
# Use PLINK to remove related individuals
#---------------------------------------
plink \
--bfile $GRMDir/17.All.filtered \
--keep $WorkingDir/output/relatedness.0.125_individuals.grm.id \
--make-bed \
--out $WorkingDir/output_PLINK/All_ind.exc_chr17

#---------------------------------------
# Filter SNPs with PLINK
#---------------------------------------
#Filter for Maf, HWE, SNPs with high missingness, individuals with high missingness, 
plink \
--bfile $WorkingDir/output_PLINK/All_ind.exc_chr17 \
--maf 0.01 \
--hwe 0.000001 \
--geno 0.01 \
--mind 0.01 \
--make-bed \
--out $WorkingDir/output_PLINK/All_ind.exc_PLINK_chr17

# Check and remove ambiguous SNPs
awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || ($5 == "C" && $6 == "G") || ($5 == "G" && $6 == "C") {print $2}' $WorkingDir/output_PLINK/All_ind.exc_PLINK_chr17.bim > $WorkingDir/output_PLINK/ambiguous_snps_chr17.txt

plink \
--bfile $WorkingDir/output_PLINK/All_ind.exc_PLINK_chr17 \
--exclude $WorkingDir/output_PLINK/ambiguous_snps_chr17.txt \
--make-bed \
--out $WorkingDir/output_PLINK/All.chr17.filtered




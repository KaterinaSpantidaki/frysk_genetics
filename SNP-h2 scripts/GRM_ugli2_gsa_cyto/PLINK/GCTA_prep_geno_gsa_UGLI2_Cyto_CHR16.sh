#!/bin/sh
#SBATCH --job-name=chr16_GCTA_prep_geno_input
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=40000

#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 9 April 2023.
# title: "Prep of genotype input files for GCTA".
# date: 10 April 2025.
# data: Lifelines 
#---
# Script to prepare the genotype plink input fils for the GRM construction GCTA - ALL CHROMOSOMES, of language learning phenotypic data.
#---


#---------------------------------------
# Overview of steps
#---------------------------------------

#1. Filter GSA vcf, UGLI2 vcf and CYTO vcf for imputation quality and MAF.
#2. Convert both to PLINK best guess data
#3. Merge, keeping only overlapping SNPs
#4. Filter merged file on HWE and MAF


#---------------------------------------
# Defining directories and paths and load modules.
#---------------------------------------

# Set UGLI2, GSA and CYTO vcf files directory. 
VCFInputDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/VCFs"  

# Set working directory.
WorkingDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/GRMs"

# Load necessary modules. 
module load PLINK/2.0-alpha4.5-20230813
module load BCFtools/1.17-GCCcore-11.3.0
module load PLINK/1.9-beta6-20190617

#---------------------------------------
# Filter UGLI2, GSA and CYTO vcfs for imputation quality and MAF.
#---------------------------------------

bcftools view -i 'INFO>.8 & MAF>.01' -f .,PASS -Oz $VCFInputDir/16.UGLI2_r2.vcf.gz > $WorkingDir/16.UGLI2_r2.filtered.vcf.gz

bcftools view -i 'INFO>.8 & MAF>.01' -f .,PASS -Oz $VCFInputDir/16.GSAr2.vcf.gz > $WorkingDir/16.GSAr2.filtered.vcf.gz

bcftools view -i 'INFO>.8 & MAF>.01' -f .,PASS -Oz $VCFInputDir/16.pbwt_reference_impute_renamed.vcf.gz > $WorkingDir/16.pbwt_reference_impute_renamed.filtered.vcf.gz

#---------------------------------------
# Convert to PLINK bfiles
#---------------------------------------

plink --vcf $WorkingDir/16.UGLI2_r2.filtered.vcf.gz --make-bed --out $WorkingDir/16.UGLI2_r2
# 577490 variants loaded from .bim file.
# 28149 people (0 males, 0 females, 28149 ambiguous) 

plink --vcf $WorkingDir/16.GSAr2.filtered.vcf.gz --make-bed --out $WorkingDir/16.GSAr2
# 571135 variants loaded from .bim file.
# 36339 people (0 males, 0 females, 36339 ambiguous) loaded from .fam.

plink --vcf $WorkingDir/16.pbwt_reference_impute_renamed.filtered.vcf.gz --make-bed --out $WorkingDir/16.pbwt_reference_impute_renamed
# 519925 variants loaded from .bim file.
# 15422 people (0 males, 0 females, 15422 ambiguous) loaded from .fam.

#---------------------------------------
# Merge with UGLI2, GSA and CYTO data.
#---------------------------------------

# Dry run bmerge for UGLI2 and GSA to identify SNPs PLINK will fail on.
plink \
--bfile $WorkingDir/16.UGLI2_r2 \
--bmerge $WorkingDir/16.GSAr2 \
--merge-mode 6 \
--out $WorkingDir/16.UGLI2_GSA.merge.failures	

# Select variants with multiple positions
# Add them to the list of missnps for UGLI2 and GSA.
fgrep \'rs $WorkingDir/16.UGLI2_GSA.merge.failures.log |\
awk '{print $7}' |\
sed "s/'//g" | sed "s/\.//g" >> $WorkingDir/16.UGLI2_GSA.merge.failures.missnp   

# Remove variants with multiple positions and missnps from UGLI2 and GSA.
plink \
--bfile $WorkingDir/16.GSAr2 \
--exclude $WorkingDir/16.UGLI2_GSA.merge.failures.missnp  \
--make-bed \
--out $WorkingDir/16.GSAr2.for_merge
# 552659/571135 variants pass filters and QC.

plink \
--bfile $WorkingDir/16.UGLI2_r2 \
--exclude $WorkingDir/16.UGLI2_GSA.merge.failures.missnp  \
--make-bed \
--out $WorkingDir/16.UGLI2_r2.for_merge
# 558935/577490 variants pass filters and QC.

# Merge cleaned UGLI2 and GSA.
plink \
--bfile $WorkingDir/16.GSAr2.for_merge \
--bmerge $WorkingDir/16.UGLI2_r2.for_merge \
--make-bed \
--out $WorkingDir/16.GSAr2_UGLI2_r2.cleaned
# 565521 variants and 64488 people pass filters and QC.

# Dry run bmerge for CYTO and cleaned merged UGLI2 and GSA to identify SNPs PLINK will fail on.
plink \
--bfile $WorkingDir/16.GSAr2_UGLI2_r2.cleaned \
--bmerge $WorkingDir/16.pbwt_reference_impute_renamed \
--merge-mode 6 \
--out $WorkingDir/16.All.merge.failures	

# Select variants with multiple positions.
# Add them to the list of missnps for CYTO and cleaned merged UGLI2 and GSA.
fgrep \'rs $WorkingDir/16.All.merge.failures.log |\
awk '{print $7}' |\
sed "s/'//g" | sed "s/\.//g" >> $WorkingDir/16.All.merge.failures.missnp   

# Remove variants with multiple positions and missnps from CYTO.
plink \
--bfile $WorkingDir/16.pbwt_reference_impute_renamed \
--exclude $WorkingDir/16.All.merge.failures.missnp  \
--make-bed \
--out $WorkingDir/16.pbwt_reference_impute_renamed.for_merge
# 504183/519925 variants pass filters and QC.

# Identify overlapping variants between UGLI2, GSA and CYTO.
grep -Ff $WorkingDir/16.UGLI2_r2.for_merge.bim $WorkingDir/16.GSAr2.for_merge.bim > $WorkingDir/16.UGLI2_r2_GSAr2.overlapping.variants.txt
# 545.134

grep -Ff $WorkingDir/16.UGLI2_r2_GSAr2.overlapping.variants.txt $WorkingDir/16.pbwt_reference_impute_renamed.for_merge.bim > $WorkingDir/16.All.overlapping.variants.txt
# 497.073

# Merge UGLI2 and GSA keeping only overlapping SNPs.
plink \
--bfile $WorkingDir/16.GSAr2.for_merge \
--bmerge $WorkingDir/16.UGLI2_r2.for_merge \
--extract $WorkingDir/16.All.overlapping.variants.txt \
--make-bed \
--out $WorkingDir//16.GSAr2_UGLI2_r2

# Merge UGLI2, GSA and CYTO keeping only overlapping SNPs.
plink \
--bfile $WorkingDir/16.GSAr2_UGLI2_r2 \
--bmerge $WorkingDir/16.pbwt_reference_impute_renamed.for_merge \
--extract $WorkingDir/16.All.overlapping.variants.txt \
--make-bed \
--out $WorkingDir//16.All

#---------------------------------------
# Filter HWE and MAF.
#---------------------------------------
# Filter rare variants and HWE.
plink \
--bfile $WorkingDir/16.All \
--maf 0.01 \
--hwe 0.0001 \
--make-bed \
--out $WorkingDir/16.All.filtered
# 494.184/497.073

#---------------------------------------
# Remove non-required output files.
#---------------------------------------
# remove files not needed anymore
#rm $WorkingDir/PLINK/All.chr16.for_merge.*
#rm $WorkingDir/PLINK/GSA.chr16.for_merge.*
#rm $WorkingDir/PLINK/UGLI2.chr16.for_merge.*
#rm $WorkingDir/PLINK/All3.chr16.merge.*
#rm $WorkingDir/PLINK/UGLI1and2.chr16.merge.*
#rm $WorkingDir/PLINK/chr16.overlapping.variants.txt
#rm $WorkingDir/PLINK/UGLI2.chr16.bed
#rm $WorkingDir/PLINK/UGLI2.chr16.bim
#rm $WorkingDir/PLINK/UGLI2.chr16.fam
#rm $WorkingDir/PLINK/UGLI2.chr16.log
#rm $WorkingDir/PLINK/UGLI1and2.chr16.bed
#rm $WorkingDir/PLINK/UGLI1and2.chr16.bim
#rm $WorkingDir/PLINK/UGLI1and2.chr16.fam
#rm $WorkingDir/PLINK/UGLI1and2.chr16.log
#rm $WorkingDir/PLINK/UGLI1and2.chr16.nosex



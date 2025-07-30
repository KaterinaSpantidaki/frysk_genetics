#!/bin/sh
#SBATCH --job-name=mergeGCTA_GRM
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=200000

#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 23 April 2024.
# title: "Merge GRMs per chromosome, and create GRM based on IBS kernels".
# date: 16 April 2025.
# data: Lifelines 
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
# Merge the GRM output of each individual chromosome to one file.
#---------------------------------------
# store list of grms in textfile.
ls -d $WorkingDir/GRMs_ugli2_gsa_cyto/All_intermediate_files/GRM.*.All.filtered.grm.N.bin | sed 's/.grm.N.bin//' > $WorkingDir/GRMs_ugli2_gsa_cyto/fam.and.unrel/List_GRM_All.txt

# merge grms.
gcta64 --mgrm $WorkingDir/GRMs_ugli2_gsa_cyto/fam.and.unrel/List_GRM_All.txt --make-grm --out $WorkingDir/GRMs_ugli2_gsa_cyto/fam.and.unrel/GRM.All


#---------------------------------------
# Create an additional GRM from the GRM above to allow GREML in family data.
#---------------------------------------

gcta64 --grm $WorkingDir/GRMs_ugli2_gsa_cyto/fam.and.unrel/GRM.All --make-bK 0.02 --make-grm --out $WorkingDir/GRMs_ugli2_gsa_cyto/bK/GRM.All.bK.02
#gcta64 --grm $WorkingDir/GRM/fam.and.unrel/GRM.All --make-bK 0.05 --make-grm --out $WorkingDir/GRM/bK/GRM.All.bK.05

#---------------------------------------
# Create text file with location of full GRM and bK adjusted GRM.
#---------------------------------------
#Required for GREML heritability analyses
file_locations=(
  "/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/GRMs/GRMs_ugli2_gsa_cyto/fam.and.unrel/GRM.All"
  "/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/GRMs/GRMs_ugli2_gsa_cyto/bK/GRM.All.bK.02"
)

printf "%s\n" "${file_locations[@]}" > $WorkingDir/GRMs_ugli2_gsa_cyto/bK/GRM.bK.02.mgrm.txt



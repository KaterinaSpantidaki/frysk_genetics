---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 24 April 2024.
# title: "GRM construction - create script for all chromosomes and run all".
# date: 14 April 2025.
# data: Lifelines 
---
# Script to run all GRM generation per chromosome scripts and create all GRMs.
# Script to run per chromosome: /groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/GCTA/GRM/fam.and.unrel/GCTA_calculate_GRM_CHR.sh
---

# Set the directory for the scripts. 
ScriptDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/GCTA/GRM_ugli2_gsa_cyto/fam.and.unrel"

# Create scripts per chromosome
for i in {1..22}
do
sed 's/__CHR__/'$i'/g' $ScriptDir/GCTA_calculate_GRM_CHR.sh > $ScriptDir/GCTA_calculate_GRM_chr$i.sh
done

# Set permissions
cd $ScriptDir/
chmod 700 *.sh

# Run
for i in {1..22}
do
sbatch $ScriptDir/GCTA_calculate_GRM_chr$i.sh
done

# Check progress.
module load cluster-utils
cqueue -u umcg-aspantidaki

# Cancel jobs
#scancel -u umcg-aspantidaki


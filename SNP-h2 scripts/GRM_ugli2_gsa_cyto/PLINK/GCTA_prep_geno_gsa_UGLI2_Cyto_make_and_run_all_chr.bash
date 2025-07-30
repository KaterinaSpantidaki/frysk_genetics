---
# author: "Katerina Spantidaki"
# original author(s): Else Eising 10 February 2023, Danielle Admiraal 12 April 2023.
# title: "Make and run the prep script of the genotype input files for GCTA - ALL CHROMOSOMES".
# date: 11 April 2025.
# data: Lifelines 
---
# Script to make and run the scripts of preparing the genotype input files for GCTA - ALL CHROMOSOMES, of language learning phenotypic data.
---

# Set the directory for the scripts. 
ScriptDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/GCTA/SNP_h2/GRM_ugli2_gsa_cyto/PLINK"

# Create scripts per chromosome.
for i in {1..22} X; do
	sed 's/__CHR__/'$i'/g' $ScriptDir/GCTA_prep_geno_gsa_UGLI2_Cyto_CHR.sh > $ScriptDir/GCTA_prep_geno_gsa_UGLI2_Cyto_CHR$i.sh
done


# Set permissions.
cd $ScriptDir/
chmod 700 *.sh

# Run the scripts.
for i in {1..22} X; do
sbatch $ScriptDir/GCTA_prep_geno_gsa_UGLI2_Cyto_CHR$i.sh
done

# Check progress.
module load cluster-utils
cqueue -u umcg-aspantidaki

# Cancel jobs
#scancel -u umcg-aspantidaki




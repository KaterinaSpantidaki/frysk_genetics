#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 2 November 2024.
# title: "Prep of genotype input files for PGS per trait".
# date: 14 May 2025.
# data: Lifelines 
#---

#---------------------------------------
# Overview of steps
#---------------------------------------
#STEPS:
#1. Merging of UGLI1, UGLI2 and CYTO genotype platforms has already been done in these scripts: /groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/GCTA/SNP_h2/GRM_ugli2_gsa_cyto/PLINK/
#2. Filter out related individuals.
#3. Create PLINK files per trait.
#4. Remove files to maintain space.
#---------------------------------------
# STEP2: Filtering of related individuals
#---------------------------------------
# Remove individuals that have a first or second degree relative with the cutoff 0.125 command. 
ScriptDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/PGS/PGS_prep"
cd $ScriptDir
sbatch PGS_individuals.extract.sh
# After pruning the GRM, there are 37547 individuals (42363 individuals removed).

#---------------------------------------
# STEP3: Filtering with PLINK
#---------------------------------------
ScriptDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/PGS/PGS_prep"
cd $ScriptDir

# Create scripts per chromosome.
for i in {1..22}
do
sed -e 's/__CHR__/'"$i"'/g' -e 's/__\.All\./'"$i"'.All./g' $ScriptDir/PGS_prep_geno_data_filter.PLINK_CHR.sh > $ScriptDir/PGS_prep_geno_data_filter.PLINK_chr$i.sh
done

# set permissions
cd $ScriptDir/
chmod 700 *.sh

# run
for i in {1..22}
do
sbatch $ScriptDir/PGS_prep_geno_data_filter.PLINK_chr$i.sh
done

# Compile all chromosomal PLINK files into one.
ScriptDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/PGS/PGS_prep"
cd $ScriptDir
chmod 700 *.sh
sbatch combine_chrom.sh

#---------------------------------------
# STEP4: Create PLINK files per trait
#---------------------------------------
# Generate plink files per trait with only individuals that have scores on that trait and are unrelated (in R). 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load R.
module load R/4.2.2-foss-2022a-bare
R

# Load necessary libraries. 
library(dplyr)
library(data.table)

### Read the phenotype data of each trait and remove the individuals per trait that dont have a pheno score. The columns include: FID, IID , pheno score.
# Foreign accents understanding.
indiv_with_pheno_scores_foreign_accents_understanding <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/pheno_ugli_cyto_foreign_accents_understanding_cogn_filtered.txt", header = FALSE, col.names = c("FID", "IID", "pheno_score"))
indiv_with_pheno_scores_foreign_accents_understanding_no_nas <- indiv_with_pheno_scores_foreign_accents_understanding %>% filter(!is.na(pheno_score))
# 18084/19600, 1516 people were removed. 

# Foreign language learning.
indiv_with_pheno_scores_foreign_language_learning <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/pheno_ugli_cyto_foreign_language_learning_cogn_filtered.txt", header = FALSE, col.names = c("FID", "IID", "pheno_score"))
indiv_with_pheno_scores_foreign_language_learning_no_nas <- indiv_with_pheno_scores_foreign_language_learning %>% filter(!is.na(pheno_score))
# 14479/19600, 5121 people were removed. 

# Foreign speech imitation.
indiv_with_pheno_scores_foreign_speech_imitation <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/pheno_ugli_cyto_foreign_speech_imitation_cogn_filtered.txt", header = FALSE, col.names = c("FID", "IID", "pheno_score"))
indiv_with_pheno_scores_foreign_speech_imitation_no_nas <- indiv_with_pheno_scores_foreign_speech_imitation %>% filter(!is.na(pheno_score))
# 14467/19600, 5133 people were removed. 

# Number of fluent languages.
indiv_with_pheno_scores_n_fluent_languages <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/pheno_ugli_cyto_n_fluent_languages_cogn_filtered.txt", header = FALSE, col.names = c("FID", "IID", "pheno_score"))
indiv_with_pheno_scores_n_fluent_languages_no_nas <- indiv_with_pheno_scores_n_fluent_languages %>% filter(!is.na(pheno_score))
# 4788/19600, 14812 people were removed. 

### Read the unrelated individuals GRM data. The columns include: FID, IID.
unrelated_indiv <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_prep/output/relatedness.0.125_individuals.grm.id", header = FALSE, col.names = c("FID", "IID"))
# 37547.

### Keep only the unrelated individuals per trait.
# Foreign accents understanding.
unrel_accents <- indiv_with_pheno_scores_foreign_accents_understanding_no_nas %>%
  filter(IID %in% unrelated_indiv$IID)
# 8784/18084, 9300 people were removed. 
  
# Foreign language learning.
unrel_lang_learning <- indiv_with_pheno_scores_foreign_language_learning_no_nas %>%
  filter(IID %in% unrelated_indiv$IID)
# 7009/14479, 7470 people were removed.   

# Foreign speech imitation.
unrel_speech_imit <- indiv_with_pheno_scores_foreign_speech_imitation_no_nas %>%
  filter(IID %in% unrelated_indiv$IID)
# 7006/14467, 7461 people were removed. 

# Number of fluent languages.
unrel_fluent_langs <- indiv_with_pheno_scores_n_fluent_languages_no_nas %>%
  filter(IID %in% unrelated_indiv$IID)
# 2282/4788, 2506 people were removed. 

### Save the resulting snp list per trait. The columns include: FID, IID.
trait_dfs <- list(foreign_accents_understanding = unrel_foreign_accents, foreign_language_learning = unrel_lang_learning, foreign_speech_imitation = unrel_speech_imit, n_fluent_languages = unrel_fluent_langs)
for (trait in names(trait_dfs)) {
  data_out <- trait_dfs[[trait]] %>%
    select(FID, IID)
  write.table(
    data_out,
    paste0(trait, "_snp_list.txt"),
    row.names = FALSE,
    quote = FALSE,
    col.names = FALSE,
    sep = " ")}

### Generate Plink files per trait with only those individuals.
cd /groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/PGS/PGS_prep/

traits=(foreign_accents_understanding foreign_language_learning foreign_speech_imitation n_fluent_languages)
for trait in "${traits[@]}"; do
    sbatch plink_per.trait.sh "$trait"
	echo "$trait"
done


#---------------------------------------
# Continue with PGS analyses in this script: "/groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/PGS/PGS_run/PGS_steps.sh"
#---------------------------------------

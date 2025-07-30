---
# author: "Katerina Spantidaki"
# title: "Prepare pheno and covar files for GCTA analysis"
# date: "2025-04-01"
---
# Script to prepare the input files for the GREML analyses of foreign language acquisition phenotypic data.
---

# Project directory. 

# Nibbler
cd /groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/

# ------------------- R related ------------------- #
# Load R in Linux bash.
module load R/4.2.2-foss-2022a-bare
R

# Load necessary libraries. 
library(dplyr)
library(data.table)


### Read in the phenotype data for the four main traits. 
pheno <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/main_pheno_scores_age_gender_cogn_filtered.txt", header=TRUE)

# Read in genotyping linkage files IDs
linker_gsa <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/linkage_files/gsa_linkage_file.csv", header=TRUE, sep=",")
linker_ugli2 <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/linkage_files/affymetrix_linkage_file-v1_ov21_0398.csv", header=TRUE, sep=",")
genotyped_individuals <- unique(c(linker_ugli2$project_pseudo_id, linker_gsa$project_pseudo_id))

# Add column to phenotyped data with genotyped individuals (TRUE/FALSE)
pheno$genotyped <- pheno$project_pseudo_id %in% genotyped_individuals

# Add PCs 
#~~~~~~~~~~~~~~~
pcs_gsa <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/PCs/PCA_eur.UGLI.eigenvec")
colnames(pcs_gsa) <- c("gsa_ID","gsa_ID2","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10","PC11","PC12","PC13","PC14","PC15","PC16","PC17","PC18","PC19","PC20")
pcs_ugli2 <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/PCs/UGLI2_PCs.txt", header=TRUE)


#---------------------------------------------------------------------------
# Combine linker and pcs files of arrays
#---------------------------------------------------------------------------
# Give PC files same column names and define batch; merge files
# Keep only the first 10 PCs.
colnames(pcs_gsa)[1:12] <- colnames(pcs_ugli2)[1:12]
pcs_gsa$batch <- "gsa"
pcs_ugli2$batch <- "ugli2"
pcs <- rbind(pcs_gsa[,c(1:12, 23)], pcs_ugli2[,c(1:12, 23)])

# Repeat for linker file
colnames(linker_gsa)[2] <- "IID"
colnames(linker_ugli2)[2] <- "IID"
linker <- rbind(linker_gsa[,c(1:2)], linker_ugli2[,c(1:2)])


#---------------------------------------------------------------------------
# Combine phenotype and genotype data
#---------------------------------------------------------------------------
# Combine phenotype file, linker file, pc file, keep only those that overlap.
data1 <- merge(pcs, linker, by="IID", all=FALSE)
dim(pcs)
dim(linker)
dim(data1)

data2 <- merge(data1, pheno, by="project_pseudo_id", all=FALSE)
dim(data1) # 64553    14
dim(pheno) # 35762     9
dim(data2) # 15928    22



# Construct IID column that matches the IID naming in the GRM.
# GSA data IID is an UGLI data number for both FID and IID columns.
# UGLI2 data IID is a unique number and FID is 1. Use linker data.
data2$IID_GRM <- paste(data2$IID)
data2$FID_GRM <- 1
data2$FID_GRM[which(data2$batch == "gsa")] <- data2$FID[which(data2$batch == "gsa")]

#---------------------------------------------------------------------------
# Make covariate files 
#---------------------------------------------------------------------------

# Create binary covariate files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
covar_ugli_binary <- data2 %>% select(., c("FID_GRM","IID_GRM","gender","batch"))
write.table(covar_ugli_binary, "/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/covar_binary_ugli1and2_cogn_filtered.txt", row.names = FALSE, col.names = FALSE, quote=FALSE, sep="\t")

# Create quantitative covariate files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# Due to age that is different for every trait, every trait requires it's own quantitative covariate file.				   
covar_ugli_quant_accents_immitation_learning <- data2 %>% select("FID_GRM", "IID_GRM", "age_speech", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
write.table(covar_ugli_quant_accents_immitation_learning, "/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/covar_ugli_quant_ugli1and2_accents_immitation_learning_cogn_filtered.txt", row.names = FALSE, col.names = FALSE, quote=FALSE, sep="\t")

covar_ugli_quant_n_fluent_languages <- data2 %>% select("FID_GRM", "IID_GRM", "age_multilingualism", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
write.table(covar_ugli_quant_n_fluent_languages, "/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/covar_ugli_quant_ugli1and2_n_fluent_languages_cogn_filtered.txt", row.names = FALSE, col.names = FALSE, quote=FALSE, sep="\t")

#---------------------------------------------------------------------------
# Make pheno input file for analyses (FID_GRM, IID_GRM, pheno data) + run argument order file
#---------------------------------------------------------------------------
# Make pheno input file, write output
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
output_herit = "/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/"

trait_columns <- c("n_fluent_languages", "foreign_accents_understanding", 
                   "foreign_speech_imitation", "foreign_language_learning")
trait_subsets <- list()
for(trait in trait_columns) {
  trait_subsets[[trait]] <- data2 %>%
    select("FID_GRM", "IID_GRM", all_of(trait))
  write.table(trait_subsets[[trait]], 
              file = paste0(output_herit, "pheno_ugli_", trait, "_cogn_filtered",".txt"),  
              row.names = FALSE, 
              col.names = FALSE,
			  quote=FALSE,			  
              sep = "\t")
}

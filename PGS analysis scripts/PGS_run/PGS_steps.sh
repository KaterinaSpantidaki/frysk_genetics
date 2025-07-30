#---
# author: "Katerina Spantidaki"
# original author: Danielle Admiraal 11 November 2024.
# title: "Running of PGS analysis".
# date: 20 May 2025.
# data: Lifelines 
#---

#---------------------------------------
# Overview of the steps to follow in this scipt:
#---------------------------------------
#STEPS:
#1. Environment/software set-up
#2. Run PRScs - Adjust BETA's (standardized posterior SNP effect sizes) based on LD structure genotype.
#3. Obtain individual-level polygenic scores by concatenating output files from all chromosomes with plink.
#4. Run regression between individual PGS scores and phenotype scores.

#---------------------------------------
# STEP1: Setup of required environment.
#---------------------------------------
# Obtain 1000 genomes ref panel and then drag it in your nibbler directory.
cd /data/workspaces/lag/workspaces/lg-lifelines/working_data/SCRIPTS/PGS/frysk/PGS_run/lang_soc
cp /data/workspaces/lag/workspaces/lg-lifelines/working_data/SCRIPTS/PGS/wils/PGS/ref_panels/ldblk_1kg_eur.tar.gz .
tar -zxvf ldblk_1kg_eur.tar.gz

# Clone PRScs software into your working directory and then drag it in your nibbler directory.
cd /data/workspaces/lag/workspaces/lg-lifelines/working_data/SCRIPTS/PGS/frysk/PGS_run/
git clone https://github.com/getian107/PRScs.git

#---------------------------------------
# STEP2: Run PRScs - Adjust BETA's (standardized posterior SNP effect sizes) based on LD structure genotype.
#---------------------------------------
cd /groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/PGS/PGS_run/

# Summary statistics.
predictor=(schoolgrades_E4_sumstats)  

# Phenotype data of lifelines.
outcome=(foreign_accents_understanding foreign_language_learning foreign_speech_imitation n_fluent_languages) 

# Sample size of your GWAS, so for you this is the sample size of the F4 factor from Rajagopal data.
n_gwas=30982   			

for trait in "${outcome[@]}"; do
    sbatch PGS_RUN_base.sh "$predictor" "$trait" "$n_gwas"
	echo predicting "$trait" with "$predictor"
done

#---------------------------------------
# STEP3: Obtain individual-level polygenic scores by concatenating output files from all chromosomes with plink.
#---------------------------------------
cd /groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/PGS/PGS_run/

# Summary statistics.
predictor=(schoolgrades_E4_sumstats)  

# Phenotype data of lifelines.
outcome=(foreign_accents_understanding foreign_language_learning foreign_speech_imitation n_fluent_languages) 

for trait in "${outcome[@]}"; do
    sbatch PGS_plink_scores.sh "$predictor" "$trait" 
	echo predicting "$trait" with "$predictor"
done

#---------------------------------------
# STEP4: PGS regression analysis between individual PRScs (per trait) and phenotype scores - base analysis (in R).
#---------------------------------------
# Based on manual Reyna: https://www.notion.so/PRS-Manual-Clapbeat-d644fc46e59b41229d178d5fb03fb73b
# Based on ordinal logistic regression video: https://www.youtube.com/watch?v=rrRrI9gElYA
# Load the latest version of R since the rcompanion package needs R (≥ 4.4.0).
module load R/4.4.0-foss-2022a-bare 
R

# Load necessary libraries. 
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ordinal)  #ordinal regression package
library(rcompanion) #pseudo R square 
library(MASS) #plyr method (for getting data that allows the test of proportional odds)
library(brant)# test of proportional odds

# Set working directory
setwd("/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_regression/")

### Step 1: Make a functin (data_prep) that reads in PRS scores, all pheno traits and control covariates and creates z-scores for all of them, returning a data_PRS df for each trait.
### Each data_PRS file should contain all participants per trait, their zscored PRS score, their zscored age, their zscored PCs, their zscored phenotype and their gender.
data_prep <- function(trait, age_col) {
  PRScsDir <- "/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_run/output_pgs_plink/"
  PRS_scores <- read.table(paste0(PRScsDir, "schoolgrades_E4_sumstats_", trait, ".profile"), header = TRUE)
  PRS_scores$PRS_z <- scale(PRS_scores$SCORE)
  PRS_z <- PRS_scores[, c("IID", "PRS_z")]
  data_in_scores <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/main_pheno_ugli_cyto_cogn_filtered.txt", header = TRUE)
  pheno_data <- dplyr::select(data_in_scores, IID, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10,
                              all_of(trait), all_of(age_col), gender)						
  data_PRS <- merge(pheno_data, PRS_z, by = "IID")
  for (pc in paste0("PC", 1:10)) {
    data_PRS[[paste0("z", pc)]] <- scale(data_PRS[[pc]])
  } 
  data_PRS$z_age <- scale(data_PRS[[age_col]])
  return(data_PRS)
}

# Run that funciton for all traits and save the resulting dfs in a list (data_PRS_list). Make sure to include the right age covariates per trait.
traits <- c("foreign_accents_understanding", "foreign_speech_imitation", "foreign_language_learning", "n_fluent_languages")
data_PRS_list <- list()
for (trait in traits) {
  if (trait == "n_fluent_languages") {
    age_col <- "age_multilingualism"
  } else {
    age_col <- "age_speech"
  }
  data_PRS_list[[trait]] <- data_prep(trait = trait, age_col = age_col)
}

# Save each PGS_trait df individually. 
for (trait in names(data_PRS_list)) {
  write.csv(data_PRS_list[[trait]], file = paste0("/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_run/PGS_dfs/PGS_", trait, ".csv"), row.names = FALSE)}
     
# Plot distribution of PRS zscores to visualise and see any discrepancies.
PRS_z_distribution_plots <- list()
for (trait in names(data_PRS_list)) {
  data_PRS <- data_PRS_list[[trait]]
  p <- ggplot(data_PRS, aes(x = PRS_z)) +
    geom_histogram(binwidth = 0.1) +
    labs(x = "PRS Z-score", title = paste(trait)) +
    theme_minimal()
  PRS_z_distribution_plots[[trait]] <- p
}
pdf("PRS_z_distribution_plots_language_acquisition.pdf", width = 10, height = 8)
grid.arrange(grobs = PRS_z_distribution_plots, ncol = 2, nrow = 2)
dev.off()


### Step 2: Make a functin (data_prep) that runs ordinal logistic regression analysis for continous variables as predictors (PRS_z of each trait) and ordinal traits (language acquisition) as outcomes. 
# see: https://www.bookdown.org/rwnahhas/RMPH/blr-ordinal.html
ordinal_logistic_regression <- function(trait) {
  data_PRS <- data_PRS_list[[trait]]  
  
  # Convert the outcome column to an ordered factor (ordinal)
  data_PRS[[trait]] <- factor(data_PRS[[trait]], ordered = TRUE)
  
  # Define the null model as a control to check how much of the variance in an outcome trait can be explained using known covariates only as predictors 
  # (demographic (gender),population structure (e.g., zPC1 to zPC5), age (z_age)), without any genetic information (PRS_z). 
  model_null <- clm((data_PRS[[trait]])~gender + zPC1 + zPC2 + zPC3 + zPC4 + zPC5 + z_age, data = data_PRS, link = "logit")
  
  # Define the full model to check if adding genetic predisposition (as captured by the PRS_z) significantly contributes to the explained variance in an outcome trait above and beyond the covariates.  
  model_PGS <- clm((data_PRS[[trait]])~PRS_z + gender + zPC1 + zPC2 + zPC3 + zPC4 + zPC5 + z_age, data = data_PRS, link = "logit")

  # Take the coefficients from the test and compute the exponential function. Then extracts the odds ratio for PRS_z.
  OR_exp <- exp(coef(model_PGS))  
  OR_PGS <- OR_exp["PRS_z"]       

  # Calculates the 95% confidence intervals, and compute these into exponential function.
  CI_exp <- exp(confint(model_PGS)) 
  OR_CI_lower <- CI_exp["PRS_z",1]
  OR_CI_upper <- CI_exp["PRS_z",2]
  
  # Compare model fit improvement using an anova; significant means that the models are different.
  # see: https://www.r-bloggers.com/2015/08/evaluating-logistic-regression-models/
  anova_test <- anova(model_null, model_PGS)
  Chi_2_p <- anova_test$`Pr(>Chisq)`[2]
  
  # Calculate McFadden’s R2: gives relative improvement in explanatory power.
  # see https://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/ 
  # McFadden's R².
  McFadden <- nagelkerke(fit  = model_PGS,
           null = model_null)  
  McFadden <- McFadden$Pseudo.R.squared.for.model.vs.null
  McFadden_R2 <- McFadden[1,1]
  
  # Perform the test of proporsional odds(parallel regression) assumption.
  modelt <- polr((data_PRS[[trait]])~PRS_z + gender + zPC1 + zPC2 + zPC3 + zPC4 + zPC5 + z_age, data = data_PRS,  Hess = TRUE)
  brant_output <- capture.output(brant(modelt))
  if (trimws(tail(brant_output, 1)) == "H0: Parallel Regression Assumption holds") {
    brant_result <- "Pass"
  } else {
    brant_result <- "Fail"
  }

  # Get the number of individuals per trait.
  n_ind <- nrow(data_PRS)
  
  # Adjust how the output file will look, including all the avobe measurements. 
  data_out <- data.frame(predictor = as.character(),
                      outcome = as.character(),
                      class_outcome = as.character(),
                      OR_PGS = as.numeric(),
                      OR_CI_lower = as.numeric(),
                      OR_CI_upper = as.numeric(),
                      Chi_2_p = as.numeric(),
                      McFadden_R2 = as.numeric(),
					  brant_test = as.character(),
                      n_ind = as.numeric())

  data_out <- data_out %>% add_row(predictor = paste0("PRS_z_", trait),
                               outcome = trait,
                               class_outcome = "ordinal",
                               OR_PGS = OR_PGS,
                               OR_CI_lower = OR_CI_lower,
                               OR_CI_upper = OR_CI_upper,
                               Chi_2_p = Chi_2_p,
                               McFadden_R2 = McFadden_R2,
							   brant_test = brant_result,
                               n_ind = n_ind)  
  return(data_out)
}

# Run ordinal regressions for each trait and collect results in one file.
PRS_z_ordinal_regression <- data.table()
for (trait in traits) {
  data_out <- ordinal_logistic_regression(trait)
  PRS_z_ordinal_regression <- rbind(PRS_z_ordinal_regression, data_out)
}

# --- This is the step you asked for: Apply the Bonferroni Correction ---
# Apply the Bonferroni correction to the Chi_2_p values and also add a column with weather the Chi_2_p_Bonferroni_adjusted values are significant. 
apply_bonferroni_to_results <- function(PRS_z_ordinal_regression, original_alpha = 0.05) {
  num_tests <- nrow(PRS_z_ordinal_regression)
  # Apply Bonferroni correction to the Chi_2_p column
  # p.adjust multiplies each p-value by the number of tests and caps it at 1.
  PRS_z_ordinal_regression$Chi_2_p_Bonferroni_adjusted <- p.adjust(PRS_z_ordinal_regression$Chi_2_p, method = "bonferroni")
  # Add a logical column to see which traits are significant after adjustment.
  PRS_z_ordinal_regression$Significant_Bonferroni <- PRS_z_ordinal_regression$Chi_2_p_Bonferroni_adjusted < original_alpha
  message(paste0("Applied Bonferroni correction for ", num_tests, " tests."))
  message(paste0("Original alpha level (family-wise error rate): ", original_alpha))
  message(paste0("Individual p-values will be considered significant if they are less than the adjusted threshold of: ", round(original_alpha / num_tests, 5)))
  return(PRS_z_ordinal_regression)
}
PRS_z_ordinal_regression_adjusted <- apply_bonferroni_to_results(PRS_z_ordinal_regression)

## Applied Bonferroni correction for 4 tests.
## Original alpha level (family-wise error rate): 0.05
## Individual p-values will be considered significant if they are less than the adjusted threshold of: 0.0125

# Adjust the formatting of the resulting regression numeric values. 
format_smart <- function(x) {
  if (is.numeric(x)) {
    sapply(x, function(v) {
      if (abs(v) < 0.001 && v != 0) {
        formatC(v, format = "e", digits = 2) # Use scientific notation for very small numbers
      } else {
        format(round(v, 3), nsmall = 2) # Round to 3 decimal places and ensure 2 decimal places display
      }
    })
  } else {
    x
  }
}

PRS_z_ordinal_regression_adjusted[] <- lapply(PRS_z_ordinal_regression_adjusted, format_smart)

# Save the final dataframe.
fwrite(PRS_z_ordinal_regression_adjusted, "PRS_z_ordinal_regression.txt", row.names = FALSE, quote = FALSE, col.names = TRUE, sep = "\t")

# Save the txt file with aligned columns in Linux.
column -t -s $'\t' PRS_z_ordinal_regression.txt > temp.txt && mv temp.txt PRS_z_ordinal_regression.txt



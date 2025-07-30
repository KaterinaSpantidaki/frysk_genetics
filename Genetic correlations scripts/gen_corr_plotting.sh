#---
# author: "Katerina Spantidaki"
# title: "Plotting of Genetic correlations results".
# date: 10 June 2025.
# data: Lifelines 
#---

# Trait abbreviations:
# Number of fluent languages (N)
# Language learning (L)
# Accents understanding (U)
# Speech imitation (I)


# ------------------------------------------------ #
# 1: Retract the rg values and the sample sizes from the output files of the genetic correlation analysis.
# ------------------- In Linux ------------------- #

# Define directory
WorkingDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/gen_cor"

# List of .hsq files withg the genetic correlations results.
files=(
  "GREML_bivar.bk02_cogn_filtered_IN.hsq"
  "GREML_bivar.bk02_cogn_filtered_LN.hsq"
  "GREML_bivar.bk02_cogn_filtered_UL.hsq"
  "GREML_bivar.bk02_cogn_filtered_LI.hsq"
  "GREML_bivar.bk02_cogn_filtered_UI.hsq"
  "GREML_bivar.bk02_cogn_filtered_UN.hsq")


# Loop through each file.
for file in "${files[@]}"; do
  input_file="$WorkingDir/output_files/$file"
  
  # Remove .hsq extension and replace dots with underscores to generate output filename.
  output_file="$WorkingDir/output_files/$(echo ${file%.hsq} | tr '.' '_').txt"

  # Convert .hsq to .txt (this step is effectively just copying the file).
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $0}' "$input_file" > "$output_file"

  # Keep only rows 1, 17, and 20.
  sed -i -n -e '1p' -e '17p' -e '20p' "$output_file"
done


# -------------------------------------------- #
# 2: Make a dataframe with the rg values, the 95% CI and the sample sizes so that you can then plot it.
# ------------------- In R ------------------- #

# Load R (≥ 4.4.0).
module load R/4.4.0-foss-2022a-bare 
R

# Load necessary libraries. 
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggstance) 
library(corrplot)
library(Hmisc)
library(DescTools)


# Set working directory
setwd("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/gen_cor/gen_corr_plots/")

# Read the genetic correlation txt files. 
GREML_bivar_bk02_cogn_filtered_IN <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/gen_cor/output_files/GREML_bivar_bk02_cogn_filtered_IN.txt", sep = "\t", header = TRUE, fill = TRUE)
GREML_bivar_bk02_cogn_filtered_LN <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/gen_cor/output_files/GREML_bivar_bk02_cogn_filtered_LN.txt", sep = "\t", header = TRUE, fill = TRUE)
GREML_bivar_bk02_cogn_filtered_UL <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/gen_cor/output_files/GREML_bivar_bk02_cogn_filtered_UL.txt", sep = "\t", header = TRUE, fill = TRUE)
GREML_bivar_bk02_cogn_filtered_LI <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/gen_cor/output_files/GREML_bivar_bk02_cogn_filtered_LI.txt", sep = "\t", header = TRUE, fill = TRUE)
GREML_bivar_bk02_cogn_filtered_UI <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/gen_cor/output_files/GREML_bivar_bk02_cogn_filtered_UI.txt", sep = "\t", header = TRUE, fill = TRUE)
GREML_bivar_bk02_cogn_filtered_UN <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/gen_cor/output_files/GREML_bivar_bk02_cogn_filtered_UN.txt", sep = "\t", header = TRUE, fill = TRUE)

# Form the language acquisition dataframe. 
GREML_bivar_bk02_language_acquisition <- data.frame(
  Trait_pair = c(
    "Speech imitation - Number of fluent languages",
    "Language learning - Number of fluent languages",
    "Accents understanding - Language learning",
    "Language learning - Speech imitation",
    "Accents understanding - Speech imitation",
    "Accents understanding - Number of fluent languages"),
  rg = c(
    GREML_bivar_bk02_cogn_filtered_IN[1, 2],
    GREML_bivar_bk02_cogn_filtered_LN[1, 2],
    GREML_bivar_bk02_cogn_filtered_UL[1, 2],
    GREML_bivar_bk02_cogn_filtered_LI[1, 2],
    GREML_bivar_bk02_cogn_filtered_UI[1, 2],
    GREML_bivar_bk02_cogn_filtered_UN[1, 2]),
  SE = c(
    GREML_bivar_bk02_cogn_filtered_IN[1, 3],
    GREML_bivar_bk02_cogn_filtered_LN[1, 3],
    GREML_bivar_bk02_cogn_filtered_UL[1, 3],
    GREML_bivar_bk02_cogn_filtered_LI[1, 3],
    GREML_bivar_bk02_cogn_filtered_UI[1, 3],
    GREML_bivar_bk02_cogn_filtered_UN[1, 3]),
  Count = c(
    GREML_bivar_bk02_cogn_filtered_IN[2, 2],
    GREML_bivar_bk02_cogn_filtered_LN[2, 2],
    GREML_bivar_bk02_cogn_filtered_UL[2, 2],
    GREML_bivar_bk02_cogn_filtered_LI[2, 2],
    GREML_bivar_bk02_cogn_filtered_UI[2, 2],
    GREML_bivar_bk02_cogn_filtered_UN[2, 2]))

# Calculate the confidence interval (CI) range based on the standard error (SE) in order to add it as errorbars in the Genetic correlations plot.
GREML_bivar_bk02_language_acquisition$CI_lower <- (GREML_bivar_bk02_language_acquisition$rg) - (1.96 * GREML_bivar_bk02_language_acquisition$SE)
GREML_bivar_bk02_language_acquisition$CI_upper <- (GREML_bivar_bk02_language_acquisition$rg) + (1.96 * GREML_bivar_bk02_language_acquisition$SE)


# -------------------------------------------- #
# 3: Make a forest plot with the genetic correlations between the 4 main traits.
# ------------------- In R ------------------- #

png("language_acquisition_gen_corr_forest_plot.png", width = 1300, height = 600)
ggplot(GREML_bivar_bk02_language_acquisition, aes(x = rg, y = Trait_pair)) +
  geom_pointrangeh(aes(xmin = CI_lower, xmax = CI_upper), color = "#800026", fill = "#800026", size = 1.5, shape = 21, stroke = 0.5, fatten = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 1) +
  scale_x_continuous(breaks = seq(-10, 10, by = 2), limits = c(-10, 10), expand = c(0, 0)) + 
  scale_y_discrete(labels = function(x) {
    idx <- match(x, GREML_bivar_bk02_language_acquisition$Trait_pair)
    parse(text = paste0("atop(bold(", shQuote(x), "), (n==", GREML_bivar_bk02_language_acquisition$Count[idx], "))"))
  }) + 
  labs(title = "Genetic Correlations between Language Acquisition Traits", 
       x = "Genetic Correlation (rg)", 
       y = "Language Acquisition Trait Pair") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title.x = element_text(margin = margin(t = 30), hjust = 0.35), 
    axis.title.y = element_text(margin = margin(r = 30)), 
    legend.title = element_text(size = 14, face = "bold", margin = margin(b = 1)),
    plot.title = element_text(face = "bold", margin = margin(b = 40)), 
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    plot.margin = unit(c(0.5, 1.5, 0.5, 0.5), "cm"))
dev.off()



# -------------------------------------------- #
# 4: Make a heatmap with the genetic correlations between the 4 main traits.
# ------------------- In R ------------------- #

# Split Trait_pair into Trait1 and Trait2.
GREML_bivar_bk02_language_acquisition <- GREML_bivar_bk02_language_acquisition %>%
  separate(Trait_pair, into = c("Trait1", "Trait2"), sep = " - ")

# Rename the 4 main traits with shorter and clearer names for the corr heatmap.   
rename_traits <- c("Accents understanding" = "foreign_accents_understanding", "Speech imitation" = "foreign_speech_imitation", "Language learning" = "foreign_language_learning", "Number of fluent languages" = "n_fluent_languages")
GREML_bivar_bk02_language_acquisition <- GREML_bivar_bk02_language_acquisition %>%
  mutate(Trait1 = recode(Trait1, !!!rename_traits), Trait2 = recode(Trait2, !!!rename_traits))

# Create symmetric matrices.
# Empty matrices.
traits <- c("foreign_speech_imitation", "foreign_language_learning", "foreign_accents_understanding", "n_fluent_languages")
n <- length(traits)
rg_mat <- matrix(NA, n, n, dimnames = list(traits, traits))
pval_mat <- matrix(1, n, n, dimnames = list(traits, traits))  # default: not significant

# Fill in values
for (i in 1:nrow(GREML_bivar_bk02_language_acquisition)) {
  t1 <- GREML_bivar_bk02_language_acquisition$Trait1[i]
  t2 <- GREML_bivar_bk02_language_acquisition$Trait2[i]
  rg_val <- GREML_bivar_bk02_language_acquisition$rg[i]
  
  # Insert into matrix symmetrically.
  rg_mat[t1, t2] <- rg_val
  rg_mat[t2, t1] <- rg_val
  
  # Determine significance: CI excludes 0.
  sig <- (GREML_bivar_bk02_language_acquisition$CI_lower[i] > 0) | (GREML_bivar_bk02_language_acquisition$CI_upper[i] < 0)
  if (sig) {
    pval_mat[t1, t2] <- 0.01
    pval_mat[t2, t1] <- 0.01
  }
}

diag(rg_mat) <- 1
diag(pval_mat) <- 0

# Make the corr heatmap.
png("language_acquisition_gen_corr_heatmap.png", width = 2000, height = 1450, res = 150)
corrplot(rg_mat, 
         is.corr = TRUE, 
         cl.pos = 'b',
         method = 'circle',
         addgrid.col = 'grey', 
         addCoef.col = 'black',
         tl.col = "black",
         diag = FALSE,
         tl.srt = 45,            
         tl.cex = 1,
         p.mat = pval_mat, 
         sig.level = 0.05,
         pch.col = '#636363',
         pch.cex = 6)
dev.off()


# -------------------------------------------- #
# 5: Make a combined heatmap with the genetic correlations and the phenotypic correlations between the 4 main traits.
# ------------------- In R ------------------- #

rownames(rg_mat) <- colnames(rg_mat) <- c(
  "Speech imitation",
  "Language learning",
  "Accents understanding",
  "Fluent languages (n)")
rownames(pval_mat) <- colnames(pval_mat) <- rownames(rg_mat)


png("combined_corrplot_half_upper_lower.png", width = 3000, height = 2000, res = 150)
par(lwd = 3)
corrplot(rg_mat, 
         is.corr = TRUE, 
         method = 'circle',
		 type = "upper",
         addgrid.col = 'goldenrod2', 
         addCoef.col = 'black',
         tl.col = "black",
		 tl.pos = 'd',
		 tl.cex = 1.5,  
		 cl.cex = 1.5,          		 
         number.cex = 2, 
         diag = TRUE,
         p.mat = pval_mat, 
         sig.level = 0.05,
         pch.col = 'white',
         pch.cex = 7,
		 bg = '#FFFDE7')

# This corr matrix was taken from the script: "/groups/umcg-lifelines/tmp02/projects/ov21_0398/Scripts/frysk/pheno_corr/phenotypic_analysis_correlations_3.sh"	 
corrplot(corr_matrix_4_main_traits_plotting, 
		 add = TRUE,
         is.corr = TRUE, 
         cl.pos = 'n',
         method = 'circle',
		 type = 'lower',
         addgrid.col = 'darkgreen', 
         addCoef.col = 'black',
		 tl.pos = 'n',      
         diag = FALSE,           
         p.mat = p_values_4_main_traits_plotting, 
         sig.level = 0.05,
		 number.cex = 2, 
         pch.col = 'white',
         pch.cex = 7,
		 bg = '#E8F5E9')
dev.off()




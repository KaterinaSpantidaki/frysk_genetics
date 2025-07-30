---
# author: "Katerina Spantidaki"
# title: "Phenotypic correlations"
# date: "2025-04-04"
---
# Script for SNP-h2 bK 0.02 plotting of language learning phenotypic data.
---

# ------------------- In Linux ------------------- #

# Define directory
WorkingDir="/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability"

# List of .hsq files (relative to output_files directory)
files=(
  "SNPh2.bK02.UGLI1and2.CYTO.foreign_speech_imitation_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.CYTO.foreign_speech_imitation_non_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.foreign_speech_imitation_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.foreign_speech_imitation_non_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.CYTO.n_fluent_languages_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.CYTO.n_fluent_languages_non_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.n_fluent_languages_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.n_fluent_languages_non_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.CYTO.foreign_accents_understanding_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.CYTO.foreign_accents_understanding_non_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.foreign_accents_understanding_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.foreign_accents_understanding_non_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.CYTO.foreign_language_learning_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.CYTO.foreign_language_learning_non_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.foreign_language_learning_cogn_filtered.hsq"
  "SNPh2.bK02.UGLI1and2.foreign_language_learning_non_cogn_filtered.hsq"
)

# Loop through each file
for file in "${files[@]}"; do
  input_file="$WorkingDir/output_files/$file"
  
  # Remove .hsq extension and replace dots with underscores to generate output filename
  output_file="$WorkingDir/retracted_output_files/$(echo ${file%.hsq} | tr '.' '_').txt"

  # Convert .hsq to .txt
  awk 'BEGIN {FS="\t"; OFS="\t"} {print $0}' "$input_file" > "$output_file"

  # Keep only the usefull lines.
  sed -i '/^$/d;9,14d' "$output_file"
done


# ------------------- In R ------------------- #

# Load R in Linux bash.
module load R/4.2.2-foss-2022a-bare
R

# Load necessary libraries. 
library(dplyr)
library(data.table)
library(ggplot2)

# ---------
### Make a dataframe with the heritabilities for the UGLI1and2 & CYTO of cognitively filtered individuals. 
SNPh2_bK02_UGLI1and2_CYTO_foreign_accents_understanding_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_CYTO_foreign_accents_understanding_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_CYTO_foreign_language_learning_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_CYTO_foreign_language_learning_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_CYTO_foreign_speech_imitation_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_CYTO_foreign_speech_imitation_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_CYTO_n_fluent_languages_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_CYTO_n_fluent_languages_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)

SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered <- data.frame(
  Trait = c("accents_ugli_cyto_cogn_filtered", "learning_ugli_cyto_cogn_filtered", "imitation_ugli_cyto_cogn_filtered", "multi_ugli_cyto_cogn_filtered"),
  SNP_h2 = c((SNPh2_bK02_UGLI1and2_CYTO_foreign_accents_understanding_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_language_learning_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_speech_imitation_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_CYTO_n_fluent_languages_cogn_filtered[5, 2])),
  SE = c((SNPh2_bK02_UGLI1and2_CYTO_foreign_accents_understanding_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_language_learning_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_speech_imitation_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_CYTO_n_fluent_languages_cogn_filtered[5, 3])),
  Count = c((SNPh2_bK02_UGLI1and2_CYTO_foreign_accents_understanding_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_language_learning_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_speech_imitation_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_CYTO_n_fluent_languages_cogn_filtered[7,2])))
  

# Calculate the confidence interval (CI) range based on the standard error (SE) for the UGLI1and2 & CYTO of cognitively filtered individuals, in order to add it as errorbars in the SNP heritability barplot.
SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered$CI_lower <- (SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered$SNP_h2) - (1.96 * SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered$SE)
SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered$CI_upper <- (SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered$SNP_h2) + (1.96 * SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered$SE)

# Save the df.
SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered <- fwrite(SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered, "SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered.csv", quote = TRUE)

# ---------
### Make a dataframe with the heritabilities for the UGLI1and2 & CYTO of non-cognitively filtered individuals. 
SNPh2_bK02_UGLI1and2_CYTO_foreign_accents_understanding_non_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_CYTO_foreign_accents_understanding_non_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_CYTO_foreign_language_learning_non_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_CYTO_foreign_language_learning_non_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_CYTO_foreign_speech_imitation_non_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_CYTO_foreign_speech_imitation_non_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_CYTO_n_fluent_languages_non_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_CYTO_n_fluent_languages_non_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)

SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_non_cogn_filtered <- data.frame(
  Trait = c("accents_ugli_cyto_non_cogn_filtered", "learning_ugli_cyto_non_cogn_filtered", "imitation_ugli_cyto_non_cogn_filtered", "multi_ugli_cyto_non_cogn_filtered"),
  SNP_h2 = c((SNPh2_bK02_UGLI1and2_CYTO_foreign_accents_understanding_non_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_language_learning_non_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_speech_imitation_non_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_CYTO_n_fluent_languages_non_cogn_filtered[5, 2])),
  SE = c((SNPh2_bK02_UGLI1and2_CYTO_foreign_accents_understanding_non_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_language_learning_non_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_speech_imitation_non_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_CYTO_n_fluent_languages_non_cogn_filtered[5, 3])),
  Count = c((SNPh2_bK02_UGLI1and2_CYTO_foreign_accents_understanding_non_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_language_learning_non_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_CYTO_foreign_speech_imitation_non_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_CYTO_n_fluent_languages_non_cogn_filtered[7,2])))
  
# Calculate the confidence interval (CI) range based on the standard error (SE) for the UGLI1and2 & CYTO of non-cognitively filtered individuals, in order to add it as errorbars in the SNP heritability barplot.
SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_non_cogn_filtered$CI_lower <- (SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_non_cogn_filtered$SNP_h2) - (1.96 * SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_non_cogn_filtered$SE)
SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_non_cogn_filtered$CI_upper <- (SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_non_cogn_filtered$SNP_h2) + (1.96 * SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_non_cogn_filtered$SE)

# ---------
### Make a dataframe with the heritabilities for the UGLI1and2 of cognitively filtered individuals. 
SNPh2_bK02_UGLI1and2_foreign_accents_understanding_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_foreign_accents_understanding_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_foreign_language_learning_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_foreign_language_learning_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_foreign_speech_imitation_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_foreign_speech_imitation_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_n_fluent_languages_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_n_fluent_languages_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)

SNPh2_bK02_UGLI1and2_language_acquisition_cogn_filtered <- data.frame(
  Trait = c("accents_ugli_cogn_filtered", "learning_ugli_cogn_filtered", "imitation_ugli_cogn_filtered", "multi_ugli_cogn_filtered"),
  SNP_h2 = c((SNPh2_bK02_UGLI1and2_foreign_accents_understanding_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_foreign_language_learning_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_foreign_speech_imitation_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_n_fluent_languages_cogn_filtered[5, 2])),
  SE = c((SNPh2_bK02_UGLI1and2_foreign_accents_understanding_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_foreign_language_learning_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_foreign_speech_imitation_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_n_fluent_languages_cogn_filtered[5, 3])),
  Count = c((SNPh2_bK02_UGLI1and2_foreign_accents_understanding_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_foreign_language_learning_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_foreign_speech_imitation_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_n_fluent_languages_cogn_filtered[7,2])))
  
# Calculate the confidence interval (CI) range based on the standard error (SE) for the UGLI1and2 of cognitively filtered individuals, in order to add it as errorbars in the SNP heritability barplot.
SNPh2_bK02_UGLI1and2_language_acquisition_cogn_filtered$CI_lower <- (SNPh2_bK02_UGLI1and2_language_acquisition_cogn_filtered$SNP_h2) - (1.96 * SNPh2_bK02_UGLI1and2_language_acquisition_cogn_filtered$SE)
SNPh2_bK02_UGLI1and2_language_acquisition_cogn_filtered$CI_upper <- (SNPh2_bK02_UGLI1and2_language_acquisition_cogn_filtered$SNP_h2) + (1.96 * SNPh2_bK02_UGLI1and2_language_acquisition_cogn_filtered$SE)

# ---------
### Make a dataframe with the heritabilities for the UGLI1and2 of non-cognitively filtered individuals. 
SNPh2_bK02_UGLI1and2_foreign_accents_understanding_non_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_foreign_accents_understanding_non_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_foreign_language_learning_non_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_foreign_language_learning_non_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_foreign_speech_imitation_non_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_foreign_speech_imitation_non_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)
SNPh2_bK02_UGLI1and2_n_fluent_languages_non_cogn_filtered <- read.table("SNPh2_bK02_UGLI1and2_n_fluent_languages_non_cogn_filtered.txt", sep = "\t", header = TRUE, fill = TRUE)

SNPh2_bK02_UGLI1and2_language_acquisition_non_cogn_filtered <- data.frame(
  Trait = c("accents_ugli_non_cogn_filtered", "learning_ugli_non_cogn_filtered", "imitation_ugli_non_cogn_filtered", "multi_ugli_non_cogn_filtered"),
  SNP_h2 = c((SNPh2_bK02_UGLI1and2_foreign_accents_understanding_non_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_foreign_language_learning_non_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_foreign_speech_imitation_non_cogn_filtered[5, 2]), (SNPh2_bK02_UGLI1and2_n_fluent_languages_non_cogn_filtered[5, 2])),
  SE = c((SNPh2_bK02_UGLI1and2_foreign_accents_understanding_non_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_foreign_language_learning_non_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_foreign_speech_imitation_non_cogn_filtered[5, 3]), (SNPh2_bK02_UGLI1and2_n_fluent_languages_non_cogn_filtered[5, 3])),
  Count = c((SNPh2_bK02_UGLI1and2_foreign_accents_understanding_non_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_foreign_language_learning_non_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_foreign_speech_imitation_non_cogn_filtered[7,2]), (SNPh2_bK02_UGLI1and2_n_fluent_languages_non_cogn_filtered[7,2])))
  
# Calculate the confidence interval (CI) range based on the standard error (SE) for the UGLI1and2 of non-cognitively filtered individuals, in order to add it as errorbars in the SNP heritability barplot.
SNPh2_bK02_UGLI1and2_language_acquisition_non_cogn_filtered$CI_lower <- (SNPh2_bK02_UGLI1and2_language_acquisition_non_cogn_filtered$SNP_h2) - (1.96 * SNPh2_bK02_UGLI1and2_language_acquisition_non_cogn_filtered$SE)
SNPh2_bK02_UGLI1and2_language_acquisition_non_cogn_filtered$CI_upper <- (SNPh2_bK02_UGLI1and2_language_acquisition_non_cogn_filtered$SNP_h2) + (1.96 * SNPh2_bK02_UGLI1and2_language_acquisition_non_cogn_filtered$SE)

# ---------
# Merge the four SNP heritability data frames.
SNPh2_bK02_language_acquisition <- bind_rows(
  SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered %>%
    select(Trait, SNP_h2, SE, Count, CI_lower, CI_upper),
  SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_non_cogn_filtered %>%
    select(Trait, SNP_h2, SE, Count, CI_lower, CI_upper),
  SNPh2_bK02_UGLI1and2_language_acquisition_cogn_filtered %>%
    select(Trait, SNP_h2, SE, Count, CI_lower, CI_upper),
  SNPh2_bK02_UGLI1and2_language_acquisition_non_cogn_filtered %>%
    select(Trait, SNP_h2, SE, Count, CI_lower, CI_upper))

# ---------
# Adjust the merged dataframe so it can be plotted.
SNPh2_bK02_language_acquisition_plotting <- SNPh2_bK02_language_acquisition %>%
  mutate(Batch = ifelse(grepl("cyto", Trait), "ugli&cyto", "ugli"),
  Cogn_filter = ifelse(grepl("non_cogn", Trait), "Including cognitively impaired", "Excluding cognitively impaired"),
  Trait = recode(Trait, 
                       "accents_ugli_cyto_cogn_filtered" = "foreign_accents_understanding",
                       "learning_ugli_cyto_cogn_filtered" = "foreign_language_learning",
                       "imitation_ugli_cyto_cogn_filtered" = "foreign_speech_imitation",
                       "multi_ugli_cyto_cogn_filtered" = "n_fluent_languages",
                       "accents_ugli_cyto_non_cogn_filtered" = "foreign_accents_understanding",
                       "learning_ugli_cyto_non_cogn_filtered" = "foreign_language_learning",
                       "imitation_ugli_cyto_non_cogn_filtered" = "foreign_speech_imitation",
                       "multi_ugli_cyto_non_cogn_filtered" = "n_fluent_languages",
                       "accents_ugli_cogn_filtered" = "foreign_accents_understanding",
                       "learning_ugli_cogn_filtered" = "foreign_language_learning",
                       "imitation_ugli_cogn_filtered" = "foreign_speech_imitation",
                       "multi_ugli_cogn_filtered" = "n_fluent_languages",
                       "accents_ugli_non_cogn_filtered" = "foreign_accents_understanding",
                       "learning_ugli_non_cogn_filtered" = "foreign_language_learning",
                       "imitation_ugli_non_cogn_filtered" = "foreign_speech_imitation",
                       "multi_ugli_non_cogn_filtered" = "n_fluent_languages")) %>%
  mutate(Trait = factor(Trait, levels = unique(Trait[order(-SNP_h2)]))) %>%
  arrange(Cogn_filter, desc(SNP_h2))  

# Make the barplot
png("SNPh2_bK02_language_acquisition_barplot.png", width = 1200, height = 800)
ggplot(SNPh2_bK02_language_acquisition_plotting, aes(x = Trait, y = SNP_h2, fill = Batch)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +  
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.2, color = "black", linewidth = 0.3, position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width = 0.8), size = 4) +  
  facet_wrap(~ Cogn_filter, scales = "free_x") +  
  labs(x = "Trait", y = "SNP-h2 (bK02)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 14), 
    axis.title.x = element_text(size = 18, margin = margin(t = 5)), 
    axis.title.y = element_text(size = 18, margin = margin(r = 20)),  
    strip.text = element_text(size = 15, face = "bold"), 
    strip.background = element_blank(),  
    legend.text = element_text(size = 16),  
    legend.title = element_text(size = 16))
dev.off()




#---------------------------------------
# Plot the SNPh2 of the UGLI1and2 & CYTO of cognitively filtered individuals. 
#---------------------------------------

# Read the df.
SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered <- fread("SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered.csv", stringsAsFactors = FALSE) 

# Adjust the df for plotting.
SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered_plotting <- SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered %>%
  mutate(Trait = recode(Trait, 
                       "accents_ugli_cyto_cogn_filtered" = "foreign_accents_understanding",
                       "learning_ugli_cyto_cogn_filtered" = "foreign_language_learning",
                       "imitation_ugli_cyto_cogn_filtered" = "foreign_speech_imitation",
                       "multi_ugli_cyto_cogn_filtered" = "n_fluent_languages")) %>%
  mutate(Trait = factor(Trait, levels = unique(Trait[order(-SNP_h2)]))) %>%
  arrange(desc(SNP_h2))
 
# Make the barplot with errorbars.
png("SNPh2_bK02_language_acquisition_barplot_errorbars.png", width = 500, height = 600)
ggplot(SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered_plotting, aes(x = Trait, y = SNP_h2)) + 
  geom_bar(stat = "identity", width = 0.5, fill = "#800026") +  
  geom_errorbar(aes(ymin = SNP_h2, ymax = CI_upper), width = 0.05, color = "black", linewidth = 0.3) +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width = 0.8), size = 4) +   
  labs(x = "Trait", y = "SNP-h2 (bK02)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 14), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14, margin = margin(r = 20)),  
    strip.text = element_text(size = 15, face = "bold"), 
    strip.background = element_blank(),  
    legend.text = element_text(size = 16),  
    legend.title = element_text(size = 16))
dev.off()

# Make the barplot with significance "*".
SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered_plotting$Significance <- ""
for (i in 1:nrow(SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered_plotting)) {
  lower <- SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered_plotting$CI_lower[i]
  upper <- SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered_plotting$CI_upper[i]
  
  # Check if CI excludes zero. If it does, it is significant.
  if (lower > 0 | upper < 0) {
    SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered_plotting$Significance[i] <- "*"
  }
}

png("SNPh2_bK02_language_acquisition_barplot_significance.png", width = 600, height = 700)
ggplot(SNPh2_bK02_UGLI1and2_CYTO_language_acquisition_cogn_filtered_plotting, aes(x = Trait, y = SNP_h2)) + 
  geom_bar(stat = "identity", width = 0.5, fill = "#800026") +  
  geom_text(aes(label = Significance), vjust = -1.2, size = 6, color = "black") +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(width = 0.8), size = 4.5) +   
  labs(x = "Trait", y = "SNP-h2 (bK02)") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  
    axis.text.y = element_text(size = 14), 
    axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14, margin = margin(r = 20)),  
    strip.text = element_text(size = 15, face = "bold"), 
    strip.background = element_blank(),  
    legend.text = element_text(size = 16),  
    legend.title = element_text(size = 16))
dev.off()

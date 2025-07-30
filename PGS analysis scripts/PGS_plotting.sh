#---
# author: "Katerina Spantidaki"
# title: "Plotting of PGS analysis results".
# date: 06 June 2025.
# data: Lifelines 
#---

#---------------------------------------
# Overview of the different type of plots:  
#---------------------------------------
#1. PGS Scores Plot
#2. Pseudo R² Plot
#3. Odds Ratio (OR) Plot

# Load R (≥ 4.4.0).
module load R/4.4.0-foss-2022a-bare 
R

# Load necessary libraries. 
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(purrr)

# Set working directory
setwd("/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_plots/")


#---------------------------------------
# 1: Plot the PGS scores to visualise how they vary across each trait category. E.g. Are there higher PGS scores associated with a higher rank in the trait?
#---------------------------------------
# Read the dfs of each trait and tis PGS score.
PGS_foreign_language_learning <- read.csv("/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_run/PGS_dfs/PGS_foreign_language_learning.csv", header = TRUE)
PGS_foreign_accents_understanding <- read.csv("/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_run/PGS_dfs/PGS_foreign_accents_understanding.csv", header = TRUE)
PGS_foreign_speech_imitation <- read.csv("/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_run/PGS_dfs/PGS_foreign_speech_imitation.csv", header = TRUE)
PGS_n_fluent_languages <- read.csv("/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_run/PGS_dfs/PGS_n_fluent_languages.csv", header = TRUE)

# Read the df for all traits and their pheno scores. 
PhenoData <- fread("/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/main_pheno_ugli_cyto_cogn_filtered.txt", header = T)

# Combine those dfs into one df that includes the pheno scores for each trait and the PGS scores for each trait. 
PRS_pheno_scores <- PhenoData %>%
  left_join(PGS_foreign_language_learning %>% select(IID, PRS_z) %>% rename(PRS_z_foreign_language_learning = PRS_z), by = "IID") %>%
  left_join(PGS_foreign_accents_understanding %>% select(IID, PRS_z) %>% rename(PRS_z_foreign_accents_understanding = PRS_z), by = "IID") %>%
  left_join(PGS_foreign_speech_imitation %>% select(IID, PRS_z) %>% rename(PRS_z_foreign_speech_imitation = PRS_z), by = "IID") %>%
  left_join(PGS_n_fluent_languages %>% select(IID, PRS_z) %>% rename(PRS_z_n_fluent_languages = PRS_z), by = "IID") %>%
  select(IID, n_fluent_languages, PRS_z_n_fluent_languages, foreign_speech_imitation, PRS_z_foreign_speech_imitation, foreign_language_learning, PRS_z_foreign_language_learning, foreign_accents_understanding, PRS_z_foreign_accents_understanding)

# Define trait columns and corresponding PRS_z columns.
traits <- c("n_fluent_languages", "foreign_speech_imitation",
            "foreign_language_learning", "foreign_accents_understanding")

prs_cols <- c("PRS_z_n_fluent_languages", "PRS_z_foreign_speech_imitation",
              "PRS_z_foreign_language_learning", "PRS_z_foreign_accents_understanding")

# Add decile columns for each trait's PRS_z score.
for (prs_col in prs_cols) {
  decile_col <- paste0(prs_col, "_decile")
  PRS_pheno_scores[[decile_col]] <- ntile(PRS_pheno_scores[[prs_col]], 10)
}

# Compute distribution of trait scores for min and max deciles only.
plot_data <- map2_df(traits, prs_cols, function(trait_col, prs_col) {
  decile_col <- paste0(prs_col, "_decile")
  
  PRS_pheno_scores %>%
    filter(!is.na(.data[[decile_col]]), !is.na(.data[[trait_col]]),
           .data[[decile_col]] %in% c(1, 10)) %>%
    group_by(Trait = trait_col, Decile = .data[[decile_col]], 
             TraitCategory = factor(.data[[trait_col]], 
                                   levels = c("1", "2", "3", "4"),
                                   labels = c("Difficult", "A bit difficult", "Average", "Easy"))) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(Trait, Decile) %>%
    mutate(Percentage = Count / sum(Count) * 100) %>%
    ungroup()
})

plot_data$Decile <- factor(plot_data$Decile, levels = c(1, 10), labels = c("Bottom 10% PGS", "Top 10% PGS"))
plot_data_others <- plot_data %>%
  filter(Trait %in% c("foreign_accents_understanding", "foreign_language_learning"))

score_colors <- c(
  "Difficult" = "#800026",     
  "A bit difficult" = "#B03A3A",
  "Average" = "#4C7D7D",       
  "Easy" = "#005F5F")

trait_labels <- c(
  "foreign_accents_understanding" = "Foreign accents understanding",
  "foreign_language_learning" = "Foreign language learning")

# Plot just the two significant traits in a stacked barplot. 
png("two_sign_traits_score_PGS_stackedbarplot.png", width = 1000, height = 1000) 
ggplot(plot_data_others, aes(x = Decile, y = Percentage, fill = TraitCategory)) +
  geom_col(position = "stack") +
  facet_wrap(~Trait, scales = "free_y", labeller = labeller(Trait = trait_labels)) +
  scale_fill_manual(values = score_colors) +
  labs(title = "Trait scores distribution across the bottom and top 10% PGS", x = NULL, y = "Individuals per trait score (%)", fill = "Trait Score") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 40)),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 18, margin = margin(r = 30)),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
	strip.text = element_text(size = 18, margin = margin(t = 20, b = 1)))
dev.off()


# Plot the distribution of the scores of onlt the two significant traits. 
traits_subset <- c("foreign_language_learning", "foreign_accents_understanding")
trait_long <- PRS_pheno_scores %>%
  select(all_of(traits_subset)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Trait",
    values_to = "Score"
  ) %>%
  filter(!is.na(Score)) %>%
  mutate(
    Score = factor(Score,
                   levels = c("4", "3", "2", "1"),
                   labels = c("Easy", "Average", "A bit difficult", "Difficult")))

plot_data <- trait_long %>%
  group_by(Trait, Score) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Trait) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

png("distribution_pheno_scores.png", width = 1000, height = 800)
ggplot(plot_data, aes(x = Score, y = Count, fill = Trait)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("foreign_language_learning" = "#800026", 
                               "foreign_accents_understanding" = "#005F5F")) +
  labs(
    title = "Distribution of trait scores",
    x = "Score",
    y = "Count",
    fill = "Trait") +
  theme_minimal() +
  theme(
   plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(b = 20)), 
   axis.text.x = element_text(size = 16, angle = 45, hjust = 1),
   axis.text.y = element_text(size = 14),                 
   axis.title.x = element_text(size = 18, margin = margin(t = 15)),  
   axis.title.y = element_text(size = 18, margin = margin(r = 15)),  
   legend.title = element_text(size = 16, face = "bold"),
   legend.text = element_text(size = 14))
dev.off()



#---------------------------------------
# 2: Plot the R2 to visualise the precentage of varation explained by my PGS E4 factor score. 
#---------------------------------------
# Read the PGS regression data generated by the PGS_steps.sh script. 
pgs_regress <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_regression/PRS_z_ordinal_regression.txt", header = TRUE)

# Convert the predictor column into foreign_over_native_language_performance_PGS PGS for all the outcome traits. 
pgs_regress[, "predictor"] <- "foreign_over_native_language_performance_PGS"  

# Convert McFadden's R² into a percentage.
pgs_regress$McFadden_R2_percent <- pgs_regress$McFadden_R2 * 100

# Define significance labels.
pgs_regress$Significance <- ifelse(pgs_regress$Chi_2_p < 0.05, "*", "")

# Make the barplot.
png("R2_foreign_over_native_language_performance_PGS_z_language_acquisition_plot.png", width = 1000, height = 800) 
max_r2 <- max(pgs_regress$McFadden_R2_percent, na.rm = TRUE)
pgs_regress$outcome <- factor(pgs_regress$outcome, levels = pgs_regress$outcome[order(pgs_regress$McFadden_R2_percent)])
ggplot(pgs_regress, aes(x = outcome, y = McFadden_R2_percent, fill = predictor)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("foreign_over_native_language_performance_PGS" = "#800026"), labels = "PGS: Foreign over\nNative language\nperformance") +
  labs(title = "Percentage of variation explained by foreign over native language performance", x = "Language acquisition trait", y = "McFadden's R² (%)", fill = "Predictor") +
  theme_minimal() +
   theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
	axis.title.x = element_text(margin = margin(t = 18)),
	axis.title.y = element_text(margin = margin(r = 18)),
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)) +
  geom_text(aes(label = Significance, y = McFadden_R2_percent + (max_r2 * 0.05)), size = 8) + 
  ylim(0, max_r2 * 1.1)
dev.off()

### Make the same plot but now adjust the p-value threshold by Bonferroni correction. ### 
# Read the PGS regression data generated by the PGS_steps.sh script. 
pgs_regress <- read.table("/groups/umcg-lifelines/tmp02/projects/ov21_0398/PGS/frysk/PGS_regression/PRS_z_ordinal_regression.txt", header = TRUE)

# Convert the predictor column into foreign_over_native_language_performance_PGS PGS for all the outcome traits. 
pgs_regress[, "predictor"] <- "foreign_over_native_language_performance_PGS"  

# Convert McFadden's R² into a percentage.
pgs_regress$McFadden_R2_percent <- pgs_regress$McFadden_R2 * 100

# Define significance labels based on the TRUE/FALSE 'Significant_Bonferroni' column.
pgs_regress$Significance <- ifelse(pgs_regress$Significant_Bonferroni, "*", "")


# Make the barplot.
png("R2_Bonferroni_adj_foreign_over_native_language_performance_PGS_z_language_acquisition_plot.png", width = 1200, height = 800) 
max_r2 <- max(pgs_regress$McFadden_R2_percent, na.rm = TRUE)
pgs_regress$outcome <- factor(pgs_regress$outcome, levels = pgs_regress$outcome[order(pgs_regress$McFadden_R2_percent)])
ggplot(pgs_regress, aes(x = outcome, y = McFadden_R2_percent, fill = predictor)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("foreign_over_native_language_performance_PGS" = "#800026"), labels = "PGS: Foreign over\nNative language\nperformance") +
  labs(title = "Percentage of variation explained by foreign over native language performance (Bonferroni corrected, n=4)", x = "Language acquisition trait", y = "McFadden's R² (%)", fill = "Predictor") +
  theme_minimal() +
   theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
	axis.title.x = element_text(margin = margin(t = 18)),
	axis.title.y = element_text(margin = margin(r = 18)),
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)) +
  geom_text(aes(label = Significance, y = McFadden_R2_percent + (max_r2 * 0.05)), size = 8) + 
  ylim(0, max_r2 * 1.1)
dev.off()

#---------------------------------------
# 3: Plot the OR to visualise the direction and strength of the effect of PGS E4 factor scores on different traits.
#---------------------------------------
# Make the forest plot.
png("OR_foreign_over_native_PGS_z_language_acquisition_plot_pointrange.png", width = 1000, height = 400)
ggplot(pgs_regress, aes(x = OR_PGS, y = outcome, color = predictor)) +
  geom_point(size = 4) +
  geom_segment(aes(x = OR_CI_lower, xend = OR_CI_upper, y = outcome, yend = outcome), linewidth  = 1.5) +
  scale_color_manual(values = c("foreign_over_native_language_performance_PGS" = "#800026"), labels = "PGS: Foreign over\nNative language\nperformance") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(breaks = seq(0.5, 1.5, by = 0.5), limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(title = "Odds Ratios associated with foreign over native language performance", x = "OR", y = "Language acquisition trait", color = "Predictor") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.title = element_text(size = 14, face = "bold", margin = margin(b = 1)),
    plot.title = element_text(face = "bold", hjust = 0.8),
        axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13))
dev.off()

### Make the same plot but now adjust the p-value threshold by Bonferroni correction. ### 
png("OR_foreign_over_native_PGS_z_Bonferroni_adj_language_acquisition_plot.png", width = 1200, height = 400)
ggplot(pgs_regress, aes(x = OR_PGS, y = outcome, color = predictor)) +
  geom_point(size = 4) +
  geom_segment(aes(x = OR_CI_lower, xend = OR_CI_upper, y = outcome, yend = outcome), linewidth  = 1.5) +
  geom_text(aes(label = Significance), nudge_y = 0.125, size = 5, color = "black", fontface = "bold") +       
  scale_color_manual(values = c("foreign_over_native_language_performance_PGS" = "#800026"), labels = "PGS: Foreign over\nNative language\nperformance") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  scale_x_continuous(breaks = seq(0.5, 1.5, by = 0.5), limits = c(0.5, 1.5), expand = c(0, 0)) +
  labs(title = "Odds Ratios associated with foreign over native language performance (Bonferroni corrected, n=4)", x = "OR", y = "Language acquisition trait", color = "Predictor") +
  theme_minimal(base_size = 15) +
  theme(
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)),
    legend.title = element_text(size = 14, face = "bold", margin = margin(b = 1)),
    plot.title = element_text(face = "bold", hjust = 0.7),
        axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13))
dev.off()


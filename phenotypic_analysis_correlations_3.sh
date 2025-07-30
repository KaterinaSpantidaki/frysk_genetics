---
# author: "Katerina Spantidaki"
# title: "Phenotypic correlations"
# date: "2025-03-03"
---
# Script to make distribution and correlation plots for the different language learning phenotypic traits.
# Mainly in R #
---

# Project directories. 
cd /groups/umcg-lifelines/tmp02/projects/ov21_0398/Pheno_tables/frysk/
cd /home/umcg-aspantidaki/phenotype_data_copy


# Load R.
module load R/4.4.0-foss-2022a-bare  
R

# Load necessary libraries. 
library(dplyr)
library(data.table)
library(ggplot2)
library(corrplot)
library(Hmisc)
library(gridExtra)
library(DescTools)

# How to instal ggplot2 package:
#install.packages("pak")
#pak::pak("tidyverse/ggplot2")

# Read all traits:
language_type_only_fluent <- fread("language_type_only_fluent_cogn_filtered.csv", stringsAsFactors = FALSE)
speech_imitation <- fread("speech_imitation.csv", stringsAsFactors = FALSE)
speech_foreignaccent <- fread("speech_foreignaccent.csv", stringsAsFactors = FALSE)
speech_foreignlanguage <- fread("speech_foreignlanguage.csv", stringsAsFactors = FALSE)
age_of_acquisition <- fread("age_of_acquisition.csv", stringsAsFactors = FALSE)
educational_attainment <- fread("educational_attainment.csv", stringsAsFactors = FALSE)
people <- fread("people.csv", stringsAsFactors = FALSE)
rfft_results_na_labelled <- fread("rfft_results_na_labelled.csv", stringsAsFactors = FALSE)
situation <- fread("situation.csv", stringsAsFactors = FALSE)
mmse_scores_filtered <- fread("mmse_scores_filtered.csv", stringsAsFactors = FALSE)
aoa_filtered <- fread("aoa_filtered.csv", stringsAsFactors = FALSE) 


# Check the multilingualism trait individually for:
language_type_only_fluent <- fread("language_type_only_fluent_cogn_filtered.csv", stringsAsFactors = FALSE)  

#  ...the distribution of the number of languages spoken fluently by gender. 
png("n_fluent_languages_by_gender.png")
ggplot(language_type_only_fluent, aes(x = n_fluent_languages, fill = gender)) +
  geom_histogram(binwidth = 1, colour = "black", alpha = 0.5, position = "dodge") + 
  geom_vline(aes(xintercept = mean(language_type_only_fluent$n_fluent_languages)),
             color = "blue", linetype = "dashed", linewidth = 1) + 
  scale_x_continuous(breaks = seq(min(language_type_only_fluent$n_fluent_languages), max(language_type_only_fluent$n_fluent_languages), by = 1)) +
  labs(x = "Amount of fluent foreign languages", y = "Count", title = "Distribution of fluent foreign languages spoken by gender") 
dev.off()

#  ...the number of males and females speaking fluently foreign languages.
png("gender_distribution_fluent_languages.png") 
table(language_type_only_fluent$gender)
ggplot(language_type_only_fluent, aes(x = gender, fill = gender)) + 
  geom_bar(colour = "black", width = 0.6) +  
  scale_fill_manual(values = c("MALE" = "lightblue", "FEMALE" = "lightpink")) +  
  labs(x = "Gender", y = "Count", title = "Number of males and females for fluent languages") 
dev.off()

#  ...the distribution of age by gender for people speaking fluently foreign languages. 
png("age_distribution_by_gender_fluent_languages.png") 
ggplot(language_type_only_fluent, aes(x = gender, y = age, fill = gender)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("FEMALE" = "lightpink", "MALE" = "lightblue")) + 
  labs(title = "Age distribution by gender for fluent languages", x = "Gender", y = "Age")
dev.off()

# Check the speech_imitation trait individually for the distribution of:
speech_imitation <- fread("speech_imitation_results_na_labelled.csv", stringsAsFactors = FALSE)

#  ...the number of people finding imitating speach easy-difficult by gender.
# In the dataset the enumarations are ("Good", "Average", "Poor", "Very poor") but here we will change them to this:
score_labels <- c("Easy", "Average", "A bit difficult", "Difficult")
png("n_imitation_ability_by_gender.png") 
ggplot(speech_imitation, aes(x = speech_imitation_adu_q_1, fill = gender)) + 
  geom_histogram(bins = length(unique(speech_imitation$speech_imitation_adu_q_1)), alpha = 0.5, colour = "black", , position = "dodge") + 
  geom_vline(aes(xintercept = mean(speech_imitation$speech_imitation_adu_q_1)), color = "blue", linetype = "dashed", linewidth = 1) + 
  scale_x_continuous(breaks = 1:4, labels = score_labels) +
  labs(x = "Imitation ability", y = "Count", title = "Distribution of foreign speech imitation ability by gender")
dev.off()

#  ...the number of people finding imitating speach easy-difficult by gender. 
png("gender_distribution_speech_imitation.png") 
table(speech_imitation$gender)
ggplot(speech_imitation, aes(x = gender, fill = gender)) + 
  geom_bar(colour = "black", width = 0.6) +  
  scale_fill_manual(values = c("MALE" = "lightblue", "FEMALE" = "lightpink")) +  
  labs(x = "Gender", y = "Count", title = "Number of males and females for speech imitation") 
dev.off()

#  ...the number of people finding imitating speach easy-difficult by age per gender. 
png("age_distribution_by_gender_speech_imitation.png") 
ggplot(speech_imitation, aes(x = gender, y = age, fill = gender)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("FEMALE" = "lightpink", "MALE" = "lightblue")) + 
  labs(title = "Age distribution by gender for speech imitation", x = "Gender", y = "Age")
dev.off()

# Check the speech_foreignaccent trait individually for the distribution of:
speech_foreignaccent <- fread("speech_foreignaccent_results_na_labelled.csv", stringsAsFactors = FALSE)

#  ...the number of people finding foreign accents easy-difficult by gender.
score_labels <- c("Easy", "Average", "A bit difficult", "Difficult")
png("n_speech_foreignaccent_by_gender.png") 
ggplot(speech_foreignaccent, aes(x = speech_foreignaccent_adu_q_1, fill = gender)) + 
  geom_histogram(bins = length(unique(speech_foreignaccent$speech_foreignaccent_adu_q_1)), alpha = 0.5, colour = "black", , position = "dodge") + 
  geom_vline(aes(xintercept = mean(speech_foreignaccent$speech_foreignaccent_adu_q_1)), color = "blue", linetype = "dashed", linewidth = 1) + 
  scale_x_continuous(breaks = 1:4, labels = score_labels) +
  labs(x = "Foreign accent understanding ability", y = "Count", title = "Distribution of foreign accent understanding ability by gender")
dev.off()

#  ...the number of people finding foreign accents easy-difficult by gender. 
png("gender_distribution_speech_foreignaccent.png") 
table(speech_foreignaccent$gender)
ggplot(speech_foreignaccent, aes(x = gender, fill = gender)) + 
  geom_bar(colour = "black", width = 0.6) +  
  scale_fill_manual(values = c("MALE" = "lightblue", "FEMALE" = "lightpink")) +  
  labs(x = "Gender", y = "Count", title = "Number of males and females for foreign accents") 
dev.off()

#  ...the number of people finding finding foreign accents easy-difficult by age per gender. 
png("age_distribution_by_gender_speech_foreignaccent.png") 
ggplot(speech_foreignaccent, aes(x = gender, y = age, fill = gender)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("FEMALE" = "lightpink", "MALE" = "lightblue")) + 
  labs(title = "Age distribution by gender for foreign accents", x = "Gender", y = "Age")
dev.off()

# Check the speech_foreignlanguage trait individually for the distribution of:
speech_foreignlanguage <- fread("speech_foreignlanguage_results_na_labelled.csv", stringsAsFactors = FALSE)

#  ...the number of people finding foreign language learning easy-difficult by gender.
score_labels <- c("Easy", "Average", "A bit difficult", "Difficult")
png("n_speech_foreignlanguage_by_gender.png") 
ggplot(speech_foreignlanguage, aes(x = speech_foreignlanguage_adu_q_1, fill = gender)) + 
  geom_histogram(bins = length(unique(speech_foreignlanguage$speech_foreignlanguage_adu_q_1)), alpha = 0.5, colour = "black", , position = "dodge") + 
  geom_vline(aes(xintercept = mean(speech_foreignlanguage$speech_foreignlanguage_adu_q_1)), color = "blue", linetype = "dashed", linewidth = 1) + 
  scale_x_continuous(breaks = 1:4, labels = score_labels) +
  labs(x = "Foreign language learning ability", y = "Count", title = "Distribution of foreign language learning ability by gender")
dev.off()

#  ...the number of people finding foreign language learning easy-difficult by gender. 
png("gender_distribution_speech_foreignlanguage.png") 
table(speech_foreignlanguage$gender)
ggplot(speech_foreignlanguage, aes(x = gender, fill = gender)) + 
  geom_bar(colour = "black", width = 0.6) +  
  scale_fill_manual(values = c("MALE" = "lightblue", "FEMALE" = "lightpink")) +  
  labs(x = "Gender", y = "Count", title = "Number of males and females for foreign language learning") 
dev.off()

#  ...the number of people finding foreign language learning easy-difficult by age per gender. 
png("age_distribution_by_gender_speech_foreignlanguage.png") 
ggplot(speech_foreignlanguage, aes(x = gender, y = age, fill = gender)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c("FEMALE" = "lightpink", "MALE" = "lightblue")) + 
  labs(title = "Age distribution by gender for foreign language learning", x = "Gender", y = "Age")
dev.off()


# ------------ Check the age distribution in some of the traits (R) ------------ # 

mmse_3a <- fread("mmse_score_adu_m_01_3a_results_na_labbeled.csv", stringsAsFactors = FALSE)
png("age_distribution_mmse.png")
ggplot(mmse_3a, aes(x=age)) +
  geom_histogram(binwidth= 5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(age)),
             color="red", linetype="dashed", linewidth=1) + 
  labs(x = "Age", y = "Count", title = "Age distribution mmse scores") +
  scale_x_continuous(breaks = seq(min(mmse_3a$age), max(mmse_3a$age), by = 5)) 
dev.off()

language_type_only_fluent <- fread("language_type_only_fluent_cogn_filtered.csv", stringsAsFactors = FALSE)  
png("age_distribution_fluent_languages.png")
ggplot(language_type_only_fluent, aes(x=age)) +
  geom_histogram(binwidth= 5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(age)),
             color="red", linetype="dashed", linewidth=1) + 
  labs(x = "Age", y = "Count", title = "Age distribution fluent languages") +
  scale_x_continuous(breaks = seq(min(language_type_only_fluent$age), max(language_type_only_fluent$age), by = 5)) 
dev.off()

aoa <- fread("age_of_acquisition_results_na_labelled.csv", stringsAsFactors = FALSE)
png("age_distribution_aoa.png")
ggplot(aoa, aes(x=age)) +
  geom_histogram(binwidth= 5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(age)),
             color="red", linetype="dashed", linewidth=1) + 
  labs(x = "Age", y = "Count", title = "Age distribution aoa") +
  scale_x_continuous(breaks = seq(min(aoa$age), max(aoa$age), by = 5)) 
dev.off()

rfft <- fread("rfft_results_na_labelled.csv", stringsAsFactors = FALSE)
png("age_distribution_rfft.png")
ggplot(rfft, aes(x=age)) +
  geom_histogram(binwidth= 5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(age)),
             color="red", linetype="dashed", linewidth=1) + 
  labs(x = "Age", y = "Count", title = "Age distribution rfft scores") +
  scale_x_continuous(breaks = seq(min(rfft$age), max(rfft$age), by = 5)) 
dev.off()

speech_imitation <- fread("speech_imitation_results_na_labelled.csv", stringsAsFactors = FALSE)
png("age_distribution_speech_imitation.png")
ggplot(speech_imitation, aes(x=age)) +
  geom_histogram(binwidth= 5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(age)),
             color="red", linetype="dashed", linewidth=1) + 
  labs(x = "Age", y = "Count", title = "Age distribution speech imitation") +
  scale_x_continuous(breaks = seq(min(speech_imitation$age), max(speech_imitation$age), by = 5)) 
dev.off()

situation <- fread("situation_results_na_labelled.csv", stringsAsFactors = FALSE)
png("age_distribution_situation.png")
ggplot(situation, aes(x=age)) +
  geom_histogram(binwidth= 5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(age)),
             color="red", linetype="dashed", linewidth=1) + 
  labs(x = "Age", y = "Count", title = "Age distribution situation") +
  scale_x_continuous(breaks = seq(min(situation$age), max(situation$age), by = 5)) 
dev.off()

educational_attainment <- fread("educational_attainment_2a_results_na_labelled.csv", stringsAsFactors = FALSE)
png("age_distribution_educational_attainment.png")
ggplot(educational_attainment, aes(x=age)) +
  geom_histogram(binwidth= 5, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(age)),
             color="red", linetype="dashed", linewidth=1) + 
  labs(x = "Age", y = "Count", title = "Age distribution educational_attainment") +
  scale_x_continuous(breaks = seq(min(educational_attainment$age), max(educational_attainment$age), by = 5)) 
dev.off()


# ------------  Make a correlation matrix between the 4 main traits (R) ------------ # 

# Make a dataframe with the four main traits.
language_type_only_fluent <- fread("language_type_only_fluent_cogn_filtered.csv", stringsAsFactors = FALSE)
speech_imitation <- fread("speech_imitation_results_na_labelled.csv", stringsAsFactors = FALSE)
speech_foreignaccent <- fread("speech_foreignaccent_results_na_labelled.csv", stringsAsFactors = FALSE)
speech_foreignlanguage <- fread("speech_foreignlanguage_results_na_labelled.csv", stringsAsFactors = FALSE)
educational_attainment <- fread("educational_attainment_2a_results_na_labelled.csv", stringsAsFactors = FALSE)

# For the speech_imitation, speech_foreignaccent and speech_foreignlanguage traits flip the scores.
# So instead of 1 --> 4 meaning Easy --> Difficult, flip it so that it means Difficult --> Easy.
speech_imitation_fliped <- speech_imitation %>%
  mutate_at(7, ~case_when(
    . == 1 ~ 4,
    . == 2 ~ 3,
    . == 3 ~ 2,
    . == 4 ~ 1,
    TRUE ~ .  ))
speech_foreignaccent_fliped <- speech_foreignaccent %>%
  mutate_at(7, ~case_when(
    . == 1 ~ 4,
    . == 2 ~ 3,
    . == 3 ~ 2,
    . == 4 ~ 1,
    TRUE ~ .  ))
speech_foreignlanguage_fliped <- speech_foreignlanguage %>%
  mutate_at(7, ~case_when(
    . == 1 ~ 4,
    . == 2 ~ 3,
    . == 3 ~ 2,
    . == 4 ~ 1,
    TRUE ~ .  ))	

# Keep only the individuals that are not cognitively impaired (from 65+ yo). 
speech_imitation_fliped_65_plus <- speech_imitation_fliped %>%
  filter(age >= 65)
speech_imitation_fliped_65_plus_filtered <- speech_imitation_fliped_65_plus %>%
  filter(project_pseudo_id %in% mmse_scores_filtered$project_pseudo_id)
speech_imitation_fliped_65_minus <- speech_imitation_fliped %>%
  filter(age < 65)
speech_imitation <- bind_rows(speech_imitation_fliped_65_plus_filtered, speech_imitation_fliped_65_minus)

speech_foreignaccent_fliped_65_plus <- speech_foreignaccent_fliped %>%
  filter(age >= 65)
speech_foreignaccent_fliped_65_plus_filtered <- speech_foreignaccent_fliped_65_plus %>%
  filter(project_pseudo_id %in% mmse_scores_filtered$project_pseudo_id)
speech_foreignaccent_fliped_65_minus <- speech_foreignaccent_fliped %>%
  filter(age < 65)
speech_foreignaccent <- bind_rows(speech_foreignaccent_fliped_65_plus_filtered, speech_foreignaccent_fliped_65_minus)

speech_foreignlanguage_fliped_65_plus <- speech_foreignlanguage_fliped %>%
  filter(age >= 65)
speech_foreignlanguage_fliped_65_plus_filtered <- speech_foreignlanguage_fliped_65_plus %>%
  filter(project_pseudo_id %in% mmse_scores_filtered$project_pseudo_id)
speech_foreignlanguage_fliped_65_minus <- speech_foreignlanguage_fliped %>%
  filter(age < 65)
speech_foreignlanguage <- bind_rows(speech_foreignlanguage_fliped_65_plus_filtered, speech_foreignlanguage_fliped_65_minus)

# Save the cognitivelly filtered 3 main traits.
fwrite(speech_imitation, "speech_imitation_cogn_filtered.csv", quote = TRUE)
fwrite(speech_foreignaccent, "speech_foreignaccent_cogn_filtered.csv", quote = TRUE)
fwrite(speech_foreignlanguage, "speech_foreignlanguage_cogn_filtered.csv", quote = TRUE)

# Make a dataframe with the four main traits. 
four_main_traits  <- language_type_only_fluent %>%
  full_join(speech_imitation, by = "project_pseudo_id") %>%
  full_join(speech_foreignlanguage, by = "project_pseudo_id") %>%
  full_join(speech_foreignaccent, by = "project_pseudo_id") %>%
  select(project_pseudo_id, variant_id.x, date.x, age.x, gender.x, zip_code.x, 
         speech_imitation_adu_q_1, speech_foreignlanguage_adu_q_1, speech_foreignaccent_adu_q_1, n_fluent_languages) %>%
  rename(variant_id = variant_id.x, date = date.x, age = age.x, gender = gender.x, zip_code = zip_code.x)
  
# Compute the Spearman correlation matrix with the p-values for all the variables (traits).
corr_matrix_4_main_traits <- rcorr(as.matrix(four_main_traits [, 7:10]), type = "spearman")

# Adust the p-values (confirmed they are < 2e-16 with the corr.test()), 
# round the correlation coefficients to 6 decimals and save the three dataframes (p, n, corr).
corr_matrix_4_main_traits$P[corr_matrix_4_main_traits$P == 0] <- "< 2e-16"
p_values_4_main_traits <- as.data.frame(corr_matrix_4_main_traits$P) # They are all significant. 
p_values_4_main_traits$traits <- rownames(p_values_4_main_traits)
rownames(p_values_4_main_traits) <- NULL
p_values_4_main_traits <- p_values_4_main_traits[, c("traits", setdiff(names(p_values_4_main_traits), "traits"))]
p_values_4_main_traits_plotting <- p_values_4_main_traits[, c("traits", setdiff(names(p_values_4_main_traits), "traits"))] # Save this for later plotting. 

n_pseudoids_4_main_traits <- as.data.frame(corr_matrix_4_main_traits$n)
n_pseudoids_4_main_traits$traits <- rownames(n_pseudoids_4_main_traits)
rownames(n_pseudoids_4_main_traits) <- NULL
n_pseudoids_4_main_traits <- n_pseudoids_4_main_traits[, c("traits", setdiff(names(n_pseudoids_4_main_traits), "traits"))] 

corr_matrix_4_main_traits_plotting <- round(corr_matrix_4_main_traits$r, 6)	# Save this for later plotting. 
corr_matrix_4_main_traits <- as.data.frame(round(corr_matrix_4_main_traits$r, 6))						
corr_matrix_4_main_traits$traits <- rownames(corr_matrix_4_main_traits)
rownames(corr_matrix_4_main_traits) <- NULL
corr_matrix_4_main_traits <- corr_matrix_4_main_traits[, c("traits", setdiff(names(corr_matrix_4_main_traits), "traits"))]

# Save the three dataframes as data table images. 
corr_matrix_4_main_traits <- tableGrob(corr_matrix_4_main_traits)
n_pseudoids_4_main_traits <- tableGrob(n_pseudoids_4_main_traits)
p_values_4_main_traits <- tableGrob(p_values_4_main_traits)

ggsave("corr_matrix_4_main_traits.png", plot = corr_matrix_4_main_traits, width = 15, height = 2, dpi = 300)
ggsave("n_pseudoids_4_main_traits.png", plot = n_pseudoids_4_main_traits, width = 15, height = 2, dpi = 300)
ggsave("p_values_4_main_traits.png", plot = p_values_4_main_traits, width = 15, height = 2, dpi = 300)


# ------------ Make a correlation heatmap between the 4 main traits (R) ------------ # 

# Rename the 4 main traits with shorter and clearer names for the corr heatmap. 
colnames(corr_matrix_4_main_traits_plotting) <- c("foreign_speech_imitation", "foreign_language_learning", "foreign_accents_understanding", "n_fluent_languages")
rownames(corr_matrix_4_main_traits_plotting) <- c("foreign_speech_imitation", "foreign_language_learning", "foreign_accents_understanding", "n_fluent_languages")

p_values_4_main_traits_plotting[p_values_4_main_traits_plotting == "< 2e-16"] <- 0.01
p_values_4_main_traits_plotting <- data.frame(lapply(p_values_4_main_traits_plotting, as.numeric))
rownames(p_values_4_main_traits_plotting) <- c("foreign_speech_imitation", "foreign_language_learning", "foreign_accents_understanding", "n_fluent_languages")
p_values_4_main_traits_plotting <- p_values_4_main_traits_plotting[, -1]
colnames(p_values_4_main_traits_plotting) <- c("foreign_speech_imitation", "foreign_language_learning", "foreign_accents_understanding", "n_fluent_languages")
p_values_4_main_traits_plotting <- as.matrix(p_values_4_main_traits_plotting)
diag(p_values_4_main_traits_plotting) <- 0

png("corr_heatmap_4_main_traits.png", width = 2000, height = 1450, res = 150)
corrplot(corr_matrix_4_main_traits_plotting, 
         is.corr = TRUE, 
         cl.pos = 'b',
         method = 'circle',
         addgrid.col = 'grey', 
         addCoef.col = 'black',
         tl.col = "black",
         diag = FALSE,
         tl.srt = 45,            
         tl.cex = 1,
         p.mat = p_values_4_main_traits_plotting, 
         sig.level = 0.05,
         pch.col = '#636363',
         pch.cex = 6)
dev.off()



# ------------ Make correlation matrices between the 4 main traits and the additional 5 questions(R) ------------ # 

# Keep only the individuals that are not cognitively impaired (from 65+ yo). 
educational_attainment_65_plus <- educational_attainment %>%
  filter(age >= 65)
educational_attainment_65_plus_filtered <- educational_attainment_65_plus %>%
  filter(project_pseudo_id %in% mmse_scores_filtered$project_pseudo_id)
educational_attainment_65_minus <- educational_attainment %>%
  filter(age < 65)
educational_attainment <- bind_rows(educational_attainment_65_plus_filtered, educational_attainment_65_minus)

aoa_filtered <- aoa_filtered %>%
  filter(project_pseudo_id %in% mmse_scores_filtered$project_pseudo_id)

# Make a dataframe with the 4 main traits, the educational_attainment trait and the aoa trait. 
four_main_traits_edu_aoa  <- educational_attainment %>%
  full_join(speech_imitation, by = "project_pseudo_id") %>%
  full_join(speech_foreignlanguage, by = "project_pseudo_id") %>%
  full_join(speech_foreignaccent, by = "project_pseudo_id") %>%
  full_join(language_type_only_fluent, by = "project_pseudo_id") %>%
  full_join(aoa_filtered, by = "project_pseudo_id") %>%
  select(project_pseudo_id, variant_id.x, date.x, age.x, gender.x, zip_code.x, 
         speech_imitation_adu_q_1, speech_foreignlanguage_adu_q_1, speech_foreignaccent_adu_q_1, n_fluent_languages, educational_attainment_adu_c_1, age_fluent_language) %>%
  rename(variant_id = variant_id.x, date = date.x, age = age.x, gender = gender.x, zip_code = zip_code.x)

# Compute the Spearman correlation matrix with the p-values for all the variables (traits).
corr_matrix_4_main_traits_edu_aoa <- rcorr(as.matrix(four_main_traits_edu_aoa [, 7:12]), type = "spearman")

# Round the correlation coefficients to 6 decimals and save the three dataframes (p, n, corr).
p_values_4_main_traits_edu_aoa <- as.data.frame(corr_matrix_4_main_traits_edu_aoa$P) # They are all significant. 
p_values_4_main_traits_edu_aoa$traits <- rownames(p_values_4_main_traits_edu_aoa)
rownames(p_values_4_main_traits_edu_aoa) <- NULL
p_values_4_main_traits_edu_aoa <- p_values_4_main_traits_edu_aoa[, c("traits", setdiff(names(p_values_4_main_traits_edu_aoa), "traits"))]

n_pseudoids_4_main_traits_edu_aoa <- as.data.frame(corr_matrix_4_main_traits_edu_aoa$n)
n_pseudoids_4_main_traits_edu_aoa$traits <- rownames(n_pseudoids_4_main_traits_edu_aoa)
rownames(n_pseudoids_4_main_traits_edu_aoa) <- NULL
n_pseudoids_4_main_traits_edu_aoa <- n_pseudoids_4_main_traits_edu_aoa[, c("traits", setdiff(names(n_pseudoids_4_main_traits_edu_aoa), "traits"))]

corr_matrix_4_main_traits_edu_aoa_plotting <- round(corr_matrix_4_main_traits_edu_aoa$r, 6)	# Save this for later plotting. 
p_values_4_main_traits_edu_aoa_plotting <- p_values_4_main_traits_edu_aoa
corr_matrix_4_main_traits_edu_aoa <- as.data.frame(round(corr_matrix_4_main_traits_edu_aoa$r, 6))						
corr_matrix_4_main_traits_edu_aoa$traits <- rownames(corr_matrix_4_main_traits_edu_aoa)
rownames(corr_matrix_4_main_traits_edu_aoa) <- NULL
corr_matrix_4_main_traits_edu_aoa <- corr_matrix_4_main_traits_edu_aoa[, c("traits", setdiff(names(corr_matrix_4_main_traits_edu_aoa), "traits"))]

# Save the three dataframes as data table images. 
corr_matrix_4_main_traits_edu_aoa <- tableGrob(corr_matrix_4_main_traits_edu_aoa)
n_pseudoids_4_main_traits_edu_aoa <- tableGrob(n_pseudoids_4_main_traits_edu_aoa)
p_values_4_main_traits_edu_aoa <- tableGrob(p_values_4_main_traits_edu_aoa)

ggsave("corr_matrix_4_main_traits_edu_aoa.png", plot = corr_matrix_4_main_traits_edu_aoa, width = 20, height = 2.5, dpi = 300)
ggsave("n_pseudoids_4_main_traits_edu_aoa.png", plot = n_pseudoids_4_main_traits_edu_aoa, width = 20, height = 2.5, dpi = 300)
ggsave("p_values_4_main_traits_edu_aoa.png", plot = p_values_4_main_traits_edu_aoa, width = 20, height = 2.5, dpi = 300)

# Adjust the matrix to be plotted to include the educational_attainment and the aoa only as extra traits (not as rows).
corr_matrix_4_main_traits_edu_aoa_plotting <- corr_matrix_4_main_traits_edu_aoa_plotting[-c(5,6), ]

# Create a correlation matrix for the 4 main traits and the situation_results trait.
# Keep only the individuals that are not cognitively impaired (from 65+ yo). 
situation <- situation %>%
  filter(project_pseudo_id %in% mmse_scores_filtered$project_pseudo_id)

four_main_traits_situation  <- situation %>%
  full_join(speech_imitation, by = "project_pseudo_id") %>%
  full_join(speech_foreignlanguage, by = "project_pseudo_id") %>%
  full_join(speech_foreignaccent, by = "project_pseudo_id") %>%
  full_join(language_type_only_fluent, by = "project_pseudo_id") %>%
  select(project_pseudo_id, variant_id.x, date.x, age.x, gender.x, zip_code.x, 
         speech_imitation_adu_q_1, speech_foreignlanguage_adu_q_1, speech_foreignaccent_adu_q_1, n_fluent_languages, bilingual_situation_adu_q_1) %>%
  rename(variant_id = variant_id.x, date = date.x, age = age.x, gender = gender.x, zip_code = zip_code.x)

traits <- c("speech_imitation_adu_q_1", "speech_foreignlanguage_adu_q_1", "speech_foreignaccent_adu_q_1", "n_fluent_languages")

results_list <- list()
for (trait in traits) { 
  # Create the contingency table for the current trait and the situation trait.
  contingency_table <- table(four_main_traits_situation[[trait]], four_main_traits_situation$bilingual_situation_adu_q_1)
 
  # Perform the Chi-Square Test of Independence.
  chi_test <- chisq.test(contingency_table)
  
  # Calculate Cramér's V for the strength of the association.
  cramer_v <- CramerV(contingency_table)

  results_list[[trait]] <- data.frame(
    Trait = trait,
    Chi_Square_Statistic = chi_test$statistic,
    P_Value = chi_test$p.value,
    Cramers_V = cramer_v
  )
}

# Combine all results into a single data frame.
corr_four_main_traits_situation <- do.call(rbind, results_list)
rownames(corr_four_main_traits_situation) <- NULL
corr_four_main_traits_situation_plotting <- corr_four_main_traits_situation
corr_four_main_traits_situation <- tableGrob(corr_four_main_traits_situation)
ggsave("corr_four_main_traits_situation.png", plot = corr_four_main_traits_situation, width = 15, height = 2, dpi = 300)

# Make a correlation matrix dataframe and a p-value datafarme. 
corr_matrix_situation <- corr_four_main_traits_situation_plotting %>%
  select(Cramers_V) %>%
  rename(situation = Cramers_V) 
rownames(corr_matrix_situation) <- rownames(corr_matrix_4_main_traits_edu_aoa_plotting)

p_values_situation <- corr_four_main_traits_situation_plotting %>%
  select(P_Value) %>%
  rename(situation = P_Value) 
rownames(p_values_situation) <- rownames(corr_matrix_4_main_traits_edu_aoa_plotting)
  
# Create a correlation matrix for the 4 main traits and the people_results trait.
# Keep only the individuals that are not cognitively impaired (from 65+ yo). 
people <- people %>%
  filter(project_pseudo_id %in% mmse_scores_filtered$project_pseudo_id)

four_main_traits_people  <- people %>%
  full_join(speech_imitation, by = "project_pseudo_id") %>%
  full_join(speech_foreignlanguage, by = "project_pseudo_id") %>%
  full_join(speech_foreignaccent, by = "project_pseudo_id") %>%
  full_join(language_type_only_fluent, by = "project_pseudo_id") %>%
  select(project_pseudo_id, variant_id.x, date.x, age.x, gender.x, zip_code.x, 
         speech_imitation_adu_q_1, speech_foreignlanguage_adu_q_1, speech_foreignaccent_adu_q_1, n_fluent_languages, bilingual_people_adu_c_1_f1) %>%
  rename(variant_id = variant_id.x, date = date.x, age = age.x, gender = gender.x, zip_code = zip_code.x)
  
traits <- c("speech_imitation_adu_q_1", "speech_foreignlanguage_adu_q_1", "speech_foreignaccent_adu_q_1", "n_fluent_languages")

results_list <- list()
for (trait in traits) { 
  # Create the contingency table for the current trait and the people trait.
  contingency_table <- table(four_main_traits_people[[trait]], four_main_traits_people$bilingual_people_adu_c_1_f1)
  
  # Remove rows and columns where the total count is zero
  cleaned_table <- contingency_table[rowSums(contingency_table) > 0, colSums(contingency_table) > 0]

  # Perform the Chi-Square Test of Independence.
  chi_test <- chisq.test(cleaned_table)
  
  # Calculate Cramér's V for the strength of the association.
  cramer_v <- CramerV(cleaned_table)

  results_list[[trait]] <- data.frame(
    Trait = trait,
    Chi_Square_Statistic = chi_test$statistic,
    P_Value = chi_test$p.value,
    Cramers_V = cramer_v
  )
}

# Combine all results into a single data frame.
corr_four_main_traits_people <- do.call(rbind, results_list)
rownames(corr_four_main_traits_people) <- NULL
corr_four_main_traits_people_plotting <- corr_four_main_traits_people
corr_four_main_traits_people <- tableGrob(corr_four_main_traits_people)
ggsave("corr_four_main_traits_people.png", plot = corr_four_main_traits_people, width = 15, height = 2, dpi = 300)

# Make a correlation matrix dataframe and a p-value datafarme. 
corr_matrix_people <- corr_four_main_traits_people_plotting %>%
  select(Cramers_V) %>%
  rename(people = Cramers_V) 
rownames(corr_matrix_people) <- rownames(corr_matrix_4_main_traits_edu_aoa_plotting)

p_values_people <- corr_four_main_traits_people_plotting %>%
  select(P_Value) %>%
  rename(people = P_Value) 
rownames(p_values_people) <- rownames(corr_matrix_4_main_traits_edu_aoa_plotting)


# ------------ Make a correlation heatmap between the 4 main traits and the additional 5 questions(R) ------------ # 

# Combine the correlation and p-value matrices of the people, situation, education attaiment, aoa and the 4 main traits into one matrix.
corr_matrix_4_main_traits_edu_aoa_plotting <- as.data.frame(corr_matrix_4_main_traits_edu_aoa_plotting)
corr_matrix_4_main_traits_edu_aoa_plotting$people <- corr_matrix_people[rownames(corr_matrix_4_main_traits_edu_aoa_plotting), "people"]
corr_matrix_4_main_traits_edu_aoa_plotting$situation <- corr_matrix_situation[rownames(corr_matrix_4_main_traits_edu_aoa_plotting), "situation"]
corr_matrix_all_traits <- as.matrix(corr_matrix_4_main_traits_edu_aoa_plotting)

p_values_4_main_traits_edu_aoa_plotting <- as.data.frame(p_values_4_main_traits_edu_aoa_plotting)
p_values_4_main_traits_edu_aoa_plotting <- p_values_4_main_traits_edu_aoa_plotting[-c(5,6), ]
rownames(p_values_4_main_traits_edu_aoa_plotting) <- p_values_4_main_traits_edu_aoa_plotting$traits  
p_values_4_main_traits_edu_aoa_plotting <- p_values_4_main_traits_edu_aoa_plotting[, -which(names(p_values_4_main_traits_edu_aoa_plotting) == "traits")]
p_values_4_main_traits_edu_aoa_plotting$people <- p_values_people[rownames(p_values_4_main_traits_edu_aoa_plotting), "people"]
p_values_4_main_traits_edu_aoa_plotting$situation <- p_values_situation[rownames(p_values_4_main_traits_edu_aoa_plotting), "situation"]
p_values_all_traits <- as.matrix(p_values_4_main_traits_edu_aoa_plotting)

# Rename the trait names with shorter and clearer names for the corr heatmap. 
colnames(corr_matrix_all_traits) <- c("foreign_speech_imitation", "foreign_language_learning", "foreign_accents_understanding", "n_fluent_languages", "educational_attainment", "aoa", "people", "situation")
rownames(corr_matrix_all_traits) <- c("foreign_speech_imitation", "foreign_language_learning", "foreign_accents_understanding", "n_fluent_languages")

colnames(p_values_all_traits) <- c("foreign_speech_imitation", "foreign_language_learning", "foreign_accents_understanding", "n_fluent_languages", "educational_attainment", "aoa", "people", "situation")
rownames(p_values_all_traits) <- c("foreign_speech_imitation", "foreign_language_learning", "foreign_accents_understanding", "n_fluent_languages")

png("corr_heatmap_all_traits.png", width = 2000, height = 1450, res = 150)
corrplot(corr_matrix_all_traits, 
         is.corr = TRUE, 
         method = 'circle',
         addgrid.col = 'grey', 
         addCoef.col = 'black',
         tl.col = "black",
		 tl.srt = 45,
		 tl.cex = 1,  
		 cl.pos = 'b',   		 
         diag = FALSE,
         p.mat = p_values_all_traits, 
         sig.level = 0.05,
         pch.col = '#636363',
         pch.cex = 4)	 
dev.off()		 




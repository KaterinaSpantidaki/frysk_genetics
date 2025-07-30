---
# author: "Katerina Spantidaki"
# title: "Phenotypic correlations"
# date: "2025-03-03"
---
# Script to make adjust/define the different language learning phenotypic traits based on their enumerations. 
# Linux and R #
---

# Project directories. 
cd /groups/umcg-lifelines/tmp01/projects/ov21_0398/Pheno_tables/Katerina
cd /home/umcg-aspantidaki/phenotype_data_copy


# Load R.
module load R/4.2.2-foss-2022a-bare
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

# Output path.
output_path <- "/groups/umcg-lifelines/tmp01/projects/ov21_0398/Pheno_tables/Katerina/"


# ------------ Check the enumerations of the rest 3 main traits (Linux) ------------ # 
language_type_only_fluent <- fread("language_type_only_fluent.csv", stringsAsFactors = FALSE)
speech_imitation <- fread("speech_imitation_results_na_labelled.csv", stringsAsFactors = FALSE)
speech_foreignaccent <- fread("speech_foreignaccent_results_na_labelled.csv", stringsAsFactors = FALSE)
speech_foreignlanguage <- fread("speech_foreignlanguage_results_na_labelled.csv", stringsAsFactors = FALSE)
 
### Find the enumerations for the speech_foreignlanguage_adu_q_1 trait. 
awk -F, '$1 == "\"speech_foreignlanguage_adu_q_1\"" {print}' speq_q_1_enumerations.csv
# "speech_foreignlanguage_adu_q_1","1","easy","makkelijk"
# "speech_foreignlanguage_adu_q_1","2","average","gemiddeld"
# "speech_foreignlanguage_adu_q_1","3","a bit difficult","een beetje moeilijk"
# "speech_foreignlanguage_adu_q_1","4","difficult","moeilijk"

# Find the unique answers in the speech_foreignlanguage_adu_q_1 trait. 
cut -d',' -f7 speech_foreignlanguage_results_na_labelled.csv | tail -n +2 | sort | uniq -c
#  3663 "1"
# 18657 "2"
#  7813 "3"
#  4993 "4"
#     1 "5"

# Remove the person that answered "5" since it is not a part of the enumarations (answers).
awk -F',' '$7 != "\"5\"" {print}' speech_foreignlanguage_results_na_labelled.csv > temp_file && mv temp_file speech_foreignlanguage_results_na_labelled.csv
# Check if he was correctly removed. 
cut -d',' -f7 speech_foreignlanguage_results_na_labelled.csv | tail -n +2 | sort | uniq -c
#  3663 "1"
# 18657 "2"
#  7813 "3"
#  4993 "4"

### Find the enumerations for the speech_foreignaccent_adu_q_1 trait. 
awk -F, '$1 == "\"speech_foreignaccent_adu_q_1\"" {print}' speq_q_1_enumerations.csv
# "speech_foreignaccent_adu_q_1","1","easy","makkelijk"
# "speech_foreignaccent_adu_q_1","2","average","gemiddeld"
# "speech_foreignaccent_adu_q_1","3","a bit difficult","een beetje moeilijk"
# "speech_foreignaccent_adu_q_1","4","difficult","moeilijk"

# Find the unique answers in the speech_foreignaccent_adu_q_1 trait. 
cut -d',' -f7 speech_foreignaccent_results_na_labelled.csv | tail -n +2 | sort | uniq -c
#  7941 "1"
# 20109 "2"
#  5345 "3"
#  1720 "4"
#     1 "5"
#     3 "6"

# Remove the people that answered "5" and "6" since they are not a part of the enumarations (answers).
awk -F',' '$7 != "\"5\"" && $7 != "\"6\"" {print}' speech_foreignaccent_results_na_labelled.csv > temp_file && mv temp_file speech_foreignaccent_results_na_labelled.csv

# Check if he was correctly removed. 
cut -d',' -f7 speech_foreignaccent_results_na_labelled.csv | tail -n +2 | sort | uniq -c
#  7941 "1"
# 20109 "2"
#  5345 "3"
#  1720 "4"

### Find the enumerations for the speech_imitation_adu_q_1 trait. 
awk -F, '$1 == "\"speech_imitation_adu_q_1\"" {print}' speq_q_1_enumerations.csv
# "speech_imitation_adu_q_1","1","good","goed"
# "speech_imitation_adu_q_1","2","average","gemiddeld"
# "speech_imitation_adu_q_1","3","poor","niet zo goed"
# "speech_imitation_adu_q_1","4","very poor","slecht"

# Find the unique answers in the speech_imitation_adu_q_1 trait. 
cut -d',' -f7 speech_imitation_results_na_labelled.csv | tail -n +2 | sort | uniq -c
#  3352 "1"
# 13910 "2"
# 12989 "3"
#  4862 "4"


# ------------ Check the enumerations of the rest 5 additional questions (Linux) ------------ # 
rfft <- fread("rfft_results_na_labelled.csv", stringsAsFactors = FALSE)
educational_attainment <- fread("educational_attainment_2a_results_na_labelled.csv", stringsAsFactors = FALSE)
situation <- fread("situation_results_na_labelled.csv", stringsAsFactors = FALSE)
aoa <- fread("age_of_acquisition_results_na_labelled.csv", stringsAsFactors = FALSE)
people <- fread("people_results_na_labelled.csv", stringsAsFactors = FALSE)

### Find the enumerations for the educational_attainment_2a_results trait. 
awk -F, '$1 == "\"educational_attainment_adu_c_1\"" {print}' 2a_q_1_enumerations.csv
# "educational_attainment_adu_c_1","1","low","laag"
# "educational_attainment_adu_c_1","2","middle","midden"
# "educational_attainment_adu_c_1","3","high","hoog"

# Find the unique answers in the educational_attainment_2a_results trait. 
cut -d',' -f7 educational_attainment_2a_results_na_labelled.csv | tail -n +2 | sort | uniq -c

# Remove the people that asnwered anything but "1", "2" or "3" since it is not a part of the enumarations (answers).
(head -n 1 educational_attainment_2a_results_na_labelled.csv && tail -n +2 educational_attainment_2a_results_na_labelled.csv | awk -F',' '$7 == "\"1\"" || $7 == "\"2\"" || $7 == "\"3\""') > temp_file && mv temp_file educational_attainment_2a_results_na_labelled.csv

# Check if he was correctly removed. 
cut -d',' -f7 educational_attainment_2a_results_na_labelled.csv | tail -n +2 | sort | uniq -c
#  26555 "1"
# 33740 "2"
# 29475 "3"

### Find the enumerations for the situation_results trait. 
awk -F, '$1 == "\"bilingual_situation_adu_q_1_a\"" {print}' musq_q_1_enumerations.csv
awk -F, '$1 == "\"bilingual_situation_adu_q_1_b\"" {print}' musq_q_1_enumerations.csv
awk -F, '$1 == "\"bilingual_situation_adu_q_1_c\"" {print}' musq_q_1_enumerations.csv
awk -F, '$1 == "\"bilingual_situation_adu_q_1_d\"" {print}' musq_q_1_enumerations.csv
awk -F, '$1 == "\"bilingual_situation_adu_q_1_e\"" {print}' musq_q_1_enumerations.csv

# "bilingual_situation_adu_q_1_a","1","at home","thuis"
# "bilingual_situation_adu_q_1_b","1","at work","op het werk"
# "bilingual_situation_adu_q_1_c","1","at my association/club","bij mijn vereniging/club"
# "bilingual_situation_adu_q_1_d","1","professional","zakelijk"
# "bilingual_situation_adu_q_1_e","1","other, namely","anders, nl."

# Replace the numbers with the actual situation. 
situation$bilingual_situation_adu_q_1_a[situation$bilingual_situation_adu_q_1_a == 1] <- "at home"
situation$bilingual_situation_adu_q_1_b[situation$bilingual_situation_adu_q_1_b == 1] <- "at work"
situation$bilingual_situation_adu_q_1_c[situation$bilingual_situation_adu_q_1_c == 1] <- "at my association/club"
situation$bilingual_situation_adu_q_1_d[situation$bilingual_situation_adu_q_1_d == 1] <- "professional"
situation$bilingual_situation_adu_q_1_e[situation$bilingual_situation_adu_q_1_e == 1] <- "other"

# Find the unique answers in the situation_results trait. 
cut -d',' -f7 situation_results_na_labelled.csv | tail -n +2 | sort | uniq -c
# 3916 "at home"
# 2626 "at my association/club"
# 1202 "at work"
# 3125 "other"
# 675 "professional"

# Replace the numbers with the actual situation. 
situation$bilingual_situation_adu_q_1_a[situation$bilingual_situation_adu_q_1_a == 1] <- "at home"
situation$bilingual_situation_adu_q_1_b[situation$bilingual_situation_adu_q_1_b == 1] <- "at work"
situation$bilingual_situation_adu_q_1_c[situation$bilingual_situation_adu_q_1_c == 1] <- "at my association/club"
situation$bilingual_situation_adu_q_1_d[situation$bilingual_situation_adu_q_1_d == 1] <- "professional"
situation$bilingual_situation_adu_q_1_e[situation$bilingual_situation_adu_q_1_e == 1] <- "other"

# For those who have answered more than one situation, make sure to include them as many times as they gave an answer. 
situation$bilingual_situation_adu_q_1 <- apply(situation[, -c(1:6)], 1, function(row) {
  na.omit(row)
})
situation_expanded <- situation[rep(1:nrow(situation), times = sapply(situation$bilingual_situation_adu_q_1, length)),]
situation_expanded$bilingual_situation_adu_q_1 <- unlist(situation$bilingual_situation_adu_q_1, use.names = FALSE)
situation <- situation_expanded[, -c(7:11)]
file_path <- file.path(output_path, "situation_results_na_labelled.csv")
fwrite(situation, file_path, quote = TRUE)

### Find the enumerations for the people_results trait. 
awk -F, '$1 == "\"bilingual_people_adu_c_1_f1a\"" {print}' musq_q_1_enumerations.csv
awk -F, '$1 == "\"bilingual_people_adu_c_1_f1b\"" {print}' musq_q_1_enumerations.csv
awk -F, '$1 == "\"bilingual_people_adu_c_1_f1c\"" {print}' musq_q_1_enumerations.csv

# "bilingual_people_adu_c_1_f1a-c","1","relatives","familieleden"
# "bilingual_people_adu_c_1_f1a-c","2","friends/acquaintances","vrienden/kennissen"
# "bilingual_people_adu_c_1_f1a-c","3","nuclear family members","gezinsleden"
# "bilingual_people_adu_c_1_f1a-c","4","colleagues","collega's"
# "bilingual_people_adu_c_1_f1a-c","5","professional partners and business partners","zakenpartner"
# "bilingual_people_adu_c_1_f1a-c","6","hobby and association","hobby en vereniging"
# "bilingual_people_adu_c_1_f1a-c","7","holiday contacts","vakantiecontacten"
# "bilingual_people_adu_c_1_f1a-c","8","contacts at (volunteer)work","contacten op (vrijwilligers)werk"
# "bilingual_people_adu_c_1_f1a-c","9","community","gemeenschap"
# "bilingual_people_adu_c_1_f1a-c","10","with customers","met klanten"
# "bilingual_people_adu_c_1_f1a-c","11","contacts at training","contacten op de opleiding"
# "bilingual_people_adu_c_1_f1a-c","12","strangers","onbekenden"
# "bilingual_people_adu_c_1_f1a-c","13","foreign speakers (in the netherlands)","anderstaligen (in nederland)"
# "bilingual_people_adu_c_1_f1a-c","14","other","anders"
# "bilingual_people_adu_c_1_f1a-c","15","incomprehensible answer","onbegrijpelijk antwoord"

# Replace the numbers with the actual people. 
mapping <- c(
  "1" = "relatives",
  "2" = "friends/acquaintances",
  "3" = "nuclear family members",
  "4" = "colleagues",
  "5" = "professional partners and business partners",
  "6" = "hobby and association",
  "7" = "holiday contacts",
  "8" = "contacts at (volunteer)work",
  "9" = "community",
  "10" = "with customers",
  "11" = "contacts at training",
  "12" = "strangers",
  "13" = "foreign speakers (in the netherlands)",
  "14" = "other",
  "15" = "incomprehensible answer")

people$bilingual_people_adu_c_1_f1a <- as.character(people$bilingual_people_adu_c_1_f1a)
people$bilingual_people_adu_c_1_f1b <- as.character(people$bilingual_people_adu_c_1_f1b)
people$bilingual_people_adu_c_1_f1c <- as.character(people$bilingual_people_adu_c_1_f1c)

people$bilingual_people_adu_c_1_f1a <- mapping[people$bilingual_people_adu_c_1_f1a]
people$bilingual_people_adu_c_1_f1b <- mapping[people$bilingual_people_adu_c_1_f1b]
people$bilingual_people_adu_c_1_f1c <- mapping[people$bilingual_people_adu_c_1_f1c]

# Find the unique answers in the people_results trait. 
cut -d',' -f7 people_results_na_labelled.csv | tail -n +2 | sort | uniq -c
#  17 "colleagues"
# 260 "community"
#  29 "contacts at training"
# 167 "contacts at (volunteer)work"
# 234 "foreign speakers (in the netherlands)"
#  80 "friends/acquaintances"
# 103 "hobby and association"
# 912 "holiday contacts"
# 183 "incomprehensible answer"
#  14 "nuclear family members"
#  16 "other"
#  10 "professional partners and business partners"
#  39 "relatives"
# 137 "with customers"

# For those who have answered more than one "people", make sure to include them as many times as they gave an answer. 
people$bilingual_people_adu_c_1_f1 <- apply(people[, -c(1:6)], 1, function(row) {
  na.omit(row)
})
people_expanded <- people[rep(1:nrow(people), times = sapply(people$bilingual_people_adu_c_1_f1, length)),]
people_expanded$bilingual_people_adu_c_1_f1 <- unlist(people$bilingual_people_adu_c_1_f1, use.names = FALSE)
people <- people_expanded[, -c(7:9)]
file_path <- file.path(output_path, "people_results_na_labelled.csv")
fwrite(people, file_path, quote = TRUE)

# Adjust the age_of_acquisition (aoa) trait.
# When did you acquire a new language? 
# --> Do people acquire a new language when they are younger or when they are older? 
# --> Does that correlate with their language skills? So are they acquiring a language better if they acquire it earlier/later in life? 
aoa <- fread("age_of_acquisition.csv", stringsAsFactors = FALSE)

# Keep only the answers that belong to fluently spoken languages (based on the multilingualism trait). 
aoa_filtered <- merge(language_type_only_fluent, aoa[, -c(2:6)], by = "project_pseudo_id", all.x = TRUE)
for (i in 1:6) {
  type_col <- paste0("language_type_adu_q_1_0", i)
  startage_col <- paste0("language_startage_adu_q_1_0", i)
  if (type_col %in% names(aoa_filtered) & startage_col %in% names(aoa_filtered)) {
    aoa_filtered[[startage_col]] <- ifelse(is.na(aoa_filtered[[type_col]]), NA, aoa_filtered[[startage_col]])
  }
}


# ------------ AOA into two age groups (infant vs everything else) ------------ # 

# Replace the ages with age groups.
# 0-2 --> 1, 3-100 --> 2). 
age_grouping <- function(x) {
  ifelse(is.na(x), NA, 
         ifelse(x >= 0 & x <= 2, 1,
         ifelse(x >= 3 & x <= 100, 2, x)))
}
aoa_filtered[, 14:19] <- lapply(aoa_filtered[, 14:19], function(col) {
  if (is.numeric(col)) age_grouping(col) else col
})
# For those who have answered more than one "aoa", make sure to include them as many times as they gave an answer. 
aoa_filtered$age_fluent_language <- apply(aoa_filtered[, -c(1:13)], 1, function(row) {
  na.omit(row)
})
aoa_expanded <- aoa_filtered[rep(1:nrow(aoa_filtered), times = sapply(aoa_filtered$age_fluent_language, length)),]
aoa_expanded$age_fluent_language <- unlist(aoa_filtered$age_fluent_language, use.names = FALSE)
aoa_filtered  <- aoa_expanded[, -c(14:19)]
aoa_filtered <- aoa_filtered[, -c(7:12)]

png("distribution_aoa_2_groups.png", width = 2000, height = 1450, res = 150)
ggplot(aoa_filtered, aes(x = factor(age_fluent_language))) +
  geom_bar(fill = "steelblue") +
  labs(
    title = "Distribution of Age of Acquisition Groups",
    x = "Group (1 = Early, 2 = Late)",
    y = "Count"
  ) +
  theme_minimal()
dev.off()
fwrite(aoa_filtered, "aoa_filtered_2_age_groups.csv", quote = TRUE)


# ------------ AOA into four age groups ------------ # 

# Replace the ages with age groups.
# 0-10 --> 1, 11-20 --> 2, 21-40 --> 3, 41-100 --> 4). 
age_grouping <- function(x) {
  ifelse(is.na(x), NA, 
         ifelse(x >= 0 & x <= 10, 1,
         ifelse(x >= 11 & x <= 20, 2,
         ifelse(x >= 21 & x <= 40, 3,
         ifelse(x >= 41 & x <= 100, 4, x)))))
}
aoa_filtered[, 14:19] <- lapply(aoa_filtered[, 14:19], function(col) {
  if (is.numeric(col)) age_grouping(col) else col
})

# For those who have answered more than one "aoa", make sure to include them as many times as they gave an answer. 
aoa_filtered$age_fluent_language <- apply(aoa_filtered[, -c(1:13)], 1, function(row) {
  na.omit(row)
})
aoa_expanded <- aoa_filtered[rep(1:nrow(aoa_filtered), times = sapply(aoa_filtered$age_fluent_language, length)),]
aoa_expanded$age_fluent_language <- unlist(aoa_filtered$age_fluent_language, use.names = FALSE)
aoa_filtered  <- aoa_expanded[, -c(14:19)]
aoa_filtered <- aoa_filtered[, -c(7:12)]


png("distribution_aoa_4_groups.png", width = 2000, height = 1450, res = 150)
ggplot(aoa_filtered, aes(x = age_fluent_language)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  labs(
    title = "Distribution of Age of Acquisition",
    x = "Age of Acquisition",
    y = "Count"
  ) +
  theme_minimal()
dev.off()	



file_path <- file.path(output_path, "aoa_filtered.csv")
fwrite(aoa_filtered, file_path, quote = TRUE)


### Find the enumerations for the rfft_results trait. 
# As an index of planning efficiency, the error ratio can be calculated. 
# The error ratio can be calculated by summing the number of unique designs and the number of perseverative errors 
# across the five different parts, and then dividing the total number of perseverative errors by the total number of unique designs. 
# The error ratio assesses the degree to which the respondent is able to minimize repetition 
# while maximizing unique productions 1). 

rfft <- fread("rfft_results_na_labelled.csv", stringsAsFactors = FALSE)

rfft  <- rfft %>%
  select(project_pseudo_id, age.x, gender.x, n_fluent_languages, speech_imitation_adu_q_1,
         speech_foreignaccent_adu_q_1, speech_foreignlanguage_adu_q_1) %>%
  rename(age = age.x, gender = gender.x)

fwrite(main_pheno_scores, "main_pheno_scores.txt", sep = "\t", na = "NA")

 




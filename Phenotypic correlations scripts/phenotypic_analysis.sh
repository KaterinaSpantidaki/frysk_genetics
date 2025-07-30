---
# author: "Katerina Spantidaki"
# title: "Phenotypic correlations"
# date: "2025-03-03"
---
# Script to get descriptive statistics, distribution plots and correlations of language learning phenotypic data.
---

# Project directories. 
# Gearshift
cd /groups/umcg-lifelines/tmp01/projects/ov21_0398/
cd /home/umcg-aspantidaki/phenotype_data_copy

# Nibbler
cd /groups/umcg-lifelines/tmp02/projects/ov21_0398/Pheno_tables/frysk/

# Output path.
output_path <- "/groups/umcg-lifelines/tmp01/projects/ov21_0398/Pheno_tables/Katerina/"

# ------------------- Linux related ------------------- #

# Some basic commands to check the files:
# -- Check the file in a data table format. 
column -s, -t file.csv | less -S #If it is csv file.
column -t -s $'\t' file.txt | less -S #If it is txt file. 

# -- Check the column names of the file.
head -n 1 file.csv

# -- Check the number of rows of a file (excluding the header). 
tail -n +2 file.csv | wc -l


# Copy the phenotype data to my own directory. 
cp -r /groups/umcg-lifelines/tmp01/projects/ov21_0398/phenotype_data/ ~/phenotype_data_copy

# Print the variable names and questions.
cut -d',' -f1 musq_q_1_variables.csv | tr -d '"' 
cut -d',' -f4 musq_q_1_variables.csv | tr -d '"'

cut -d',' -f1 1a_v_1_variables.csv | tr -d '"' 
cut -d',' -f4 1a_v_1_variables.csv | tr -d '"'

cut -d',' -f1 3a_v_1_variables.csv | tr -d '"' 
cut -d',' -f4 3a_v_1_variables.csv | tr -d '"'

# Extract each variable/trait in a seperate csv file.
awk -F',' '
NR==1 {
    # Store first 7 columns as prefix and remove double quotes
    prefix = $1 "," $2 "," $3 "," $4 "," $5 "," $6;
    
    # Loop through columns starting from 7 to store file names
    for (i=7; i<=NF; i++) {
        gsub(/"/, "", $i);  # Remove double quotes from column names
        filename = $i "_results.csv";  # Name file based on column name
        files[i] = filename;
        print prefix "," $i > files[i];  # Add header to each new file
    }
    next;  # Skip header row
}
{
    # Append each row with first 7 columns + one additional column
    row_prefix = $1 "," $2 "," $3 "," $4 "," $5 "," $6;
    for (i=7; i<=NF; i++) {
        print row_prefix "," $i >> files[i];
    }
}
' musq_q_1_results.csv

# Check one of these files randomly to make sure of its content.
column -s, -t 1a_v_1_variables.csv | less -S

# Print the column names of the dataset to doublecheck all the variables are turned into files.
head -n 1 musq_q_1_results.csv | tr ',' '\n'

# Do the same for the variables in the 3a_v_1_results.csv dataset. #
awk -F',' '
NR==1 {
    # Store first 6 columns as prefix and remove double quotes
    prefix = $1 "," $2 "," $3 "," $4 "," $5 "," $6;
    
    # Loop through columns starting from 7 to store file names
    for (i=7; i<=NF; i++) {
        gsub(/"/, "", $i);  # Remove double quotes from column names
        filename = $i "_results.csv";  # Name file based on column name
        files[i] = filename;
        print prefix "," $i > files[i];  # Add header to each new file
    }
    next;  # Skip header row
}
{
    # Append each row with first 6 columns + one additional column
    row_prefix = $1 "," $2 "," $3 "," $4 "," $5 "," $6;
    for (i=7; i<=NF; i++) {
        print row_prefix "," $i >> files[i];
    }
}
' 3a_v_1_results.csv


# ------------------- R related ------------------- #
# Load R.
module load R/4.2.2-foss-2022a-bare
R

# Load necessary libraries. 
library(dplyr)
library(data.table)
library(gridExtra)
library(ggplot2)


# For each requested trait, merge all the variables that belong to it in a seperate dataframe and label the NAs. 

# Merge all the variables that will be used to define multilingualism into one dataframe. 
patterns <- c("language_other_adu_q_1_.*_results.csv",
              "language_speak_adu_q_1_.*_results.csv",
              "language_type_adu_q_1_.*_results.csv",
              "language_typespec_adu_c_1_.*_results.csv",
              "language_understand_adu_q_1_.*_results.csv")
files <- c()
for (pattern in patterns) {
  files <- c(files, list.files(pattern = pattern))
}
df1 <- fread(files[1], stringsAsFactors = FALSE)
df1[] <- lapply(df1, as.character)
for (i in 2:length(files)) {
  df2 <- fread(files[i], stringsAsFactors = FALSE)
  df2[] <- lapply(df2, as.character)
  df2_column_7 <- df2[[7]]
  column_name <- gsub("_results.csv", "", gsub(".csv", "", basename(files[i])))
  column_name <- as.character(column_name)
  df1[, (column_name) := df2_column_7]  
}
setnames(df1, old = names(df1), new = gsub("_results", "", names(df1)))
fwrite(df1, "multilingualism_results.csv", quote = TRUE)

# Label the NAs in the multilingualism dataframe.
multilingualism_results_na_labelled <- read.csv("multilingualism_results.csv", 
                                                  na.strings = c("$4", "$5", "$6", "$7"), 
                                                  colClasses = "character")
fwrite(multilingualism_results_na_labelled, "multilingualism_results_na_labelled.csv", quote = TRUE)

### Do the same for the rest of the variables. ###
# situation_results_na_labelled.csv: contains situation of multiple languages usage variables.
# people_results_na_labelled.csv: contains social circles for multiple languages usage variables.
# age_of_acquisition_results_na_labelled.csv: contains age of acquisition (AoA) of foreign language variables.
# rfft_results_na_labelled.csv: contains executive cognitive functioning variables.
# additional_traits_results_na_labelled.csv: contains additional traits variables.



# ------------ Continue with the in-house traits (Linux) ------------ # 

# Locate the in-house variables, save them as new datasets in your project folder and label the NAs. 

# educational_attainment_adu_c_1 --> dataset_order_202301 (1a and 2a study)
# speech_foreignaccent_adu_q_1 --> dataset_order_202301_speq
# speech_imitation_adu_q_1 --> dataset_order_202301_speq
# speech_foreignlanguage_adu_q_1 --> dataset_order_202301_speq

# Print for the variable name and question.
grep "educational_attainment_adu_c_1" 1a_q_1_variables.csv | cut -d',' -f1,4

# Find the column where this variable is present and write it together with the first 6 columns in a new file. 
awk -F, '{for(i=1;i<=NF;i++) if($i ~ /educational_attainment_adu_c_1/) print i, $i}' 1a_q_1_results.csv
awk -F, '{print $1 "," $2 "," $3 "," $4 "," $5 "," $6 "," $139}' 1a_q_1_results.csv > /home/umcg-aspantidaki/phenotype_data_copy/dataset_order_202502/results/educational_attainment_adu_c_1_1a_results.csv

# Label the NAs.
sed 's/\$4/NA/g; s/\$5/NA/g; s/\$6/NA/g; s/\$7/NA/g' educational_attainment_adu_c_1_1a_results.csv > educational_attainment_1a_results_na_labelled.csv

# Remove the NAs.
awk -F',' '!/NA/' educational_attainment_1a_results_na_labelled.csv > temp_file && mv temp_file educational_attainment_1a_results_na_labelled.csv

### Do the same for the rest of the variables. ###
# educational_attainment_1a_results_na_labelled.csv: contains educational attainment variables from 1a study.
# educational_attainment_2a_results_na_labelled.csv: contains educational attainment variables from 2a study.
# speech_foreignaccent_results_na_labelled.csv: contains foreign accents understanding ability variables.  
# speech_foreignlanguage_results_na_labelled.csv: contains foreign language learging ability variables.
# speech_imitation_results_na_labelled.csv: contains imitating/mimicking people talk variables.


# ------------ Define multilingualism(mainly R) ------------ # 
# This is under phenotypic_analysis_define_multilingualism_1.sh script. 


# ------------ Continue with the mmse_score_adu_m_01 trait (R) ------------ # 

# Because mmse_score_adu_m_01 exists in both studies, 3a_v_1_results.csv and 1a_v_1_results.csv, determie which one to keep. 
# Check if the pseudoIDs are the same between the 2 mmse_score_adu_m_01 datasets.
mmse_score_adu_m_01_1a_pseudoids <- tibble(ID = read.csv("mmse_score_adu_m_01_1a_results.csv", colClasses = "character")[, 1]) %>%
  mutate(Source = "mmse_score_adu_m_01_1a_pseudoids")
nrow(mmse_score_adu_m_01_1a_pseudoids) # 13.117

mmse_score_adu_m_01_3a_pseudoids <- tibble(ID = read.csv("mmse_score_adu_m_01_3a_results_no_na.csv", colClasses = "character")[, 1]) %>%
  mutate(Source = "mmse_score_adu_m_01_3a_pseudoids")
nrow(mmse_score_adu_m_01_3a_pseudoids) # 15.858

# Check how many pseudoIDs from both mmse dataframes combined are present in both mmse_score_adu_m_01_1a_pseudoids and mmse_score_adu_m_01_3a_pseudoids dataframes.
comparison_pseudoids_df <- bind_rows(mmse_score_adu_m_01_1a_pseudoids, mmse_score_adu_m_01_3a_pseudoids) %>%
  group_by(ID) %>%
  summarise(AppearsIn = paste(unique(Source), collapse = ", ")) %>%
  ungroup()
pseudoids_in_all_dfs <- comparison_pseudoids_df %>%
  filter(AppearsIn == paste(c("mmse_score_adu_m_01_1a_pseudoids", "mmse_score_adu_m_01_3a_pseudoids"), collapse = ", "))
nrow(pseudoids_in_all_dfs) # 4.680

# Check how many pseudoIDs from both mmse dataframes combined are not present in both mmse_score_adu_m_01_1a_pseudoids and mmse_score_adu_m_01_3a_pseudoids dataframes.
pseudoids_not_in_all_dfs <- comparison_pseudoids_df %>%
  filter(AppearsIn != paste(c("mmse_score_adu_m_01_1a_pseudoids", "mmse_score_adu_m_01_3a_pseudoids"), collapse = ", "))
nrow(pseudoids_not_in_all_dfs) # 19.615

# Label and remove the the NAs from the mmse dataframe. 
mmse_score_adu_m_01_3a_results_na_labbeled <- na.omit(read.csv("mmse_score_adu_m_01_3a_results_no_na.csv", 
                                                                  na.strings = c("$4", "$5", "$6", "$7"), 
                                                                  colClasses = "character"))
fwrite(mmse_score_adu_m_01_3a_results_na_labbeled, "mmse_score_adu_m_01_3a_results_na_labbeled.csv", quote = TRUE)

# Filter for the individuals that are not cognitively impaired.
mmse_scores <- fread("mmse_score_adu_m_01_3a_results_na_labbeled.csv", stringsAsFactors = FALSE)  

# Add a new column with labels based on MMSE scores.
mmse_scores_labelled <- mmse_scores %>%
  mutate(cognitive_impairment = case_when(
    mmse_score_adu_m_01 >= 24 & mmse_score_adu_m_01 <= 30 ~ "No cognitive impairment",
    mmse_score_adu_m_01 >= 18 & mmse_score_adu_m_01 <= 23 ~ "Mild cognitive impairment",
    mmse_score_adu_m_01 >= 0 & mmse_score_adu_m_01 <= 17 ~ "Severe cognitive impairment", 
  ))

# Filter for those that are not cognitively impaired.   
mmse_scores_filtered <- mmse_scores_labelled %>%
  filter(cognitive_impairment == "No cognitive impairment")
file_path <- file.path(output_path, "mmse_scores_filtered.csv")
fwrite(mmse_scores_filtered, file_path, quote = TRUE)

# Create an overlapping table with the overlapping pseudoID counts between all the traits and the mmse score. 
language_type_fluent_pseudoids <- data.frame(ID = read.csv("language_type_only_fluent.csv", colClasses = "character")[, 1])
speech_imitation_pseudoids <- data.frame(ID = read.csv("speech_imitation.csv", colClasses = "character")[, 1])
speech_foreignaccent_pseudoids <- data.frame(ID = read.csv("speech_foreignaccent.csv", colClasses = "character")[, 1])
speech_foreignlanguage_pseudoids <- data.frame(ID = read.csv("speech_foreignlanguage.csv", colClasses = "character")[, 1])
age_of_acquisition_pseudoids <- data.frame(ID = read.csv("age_of_acquisition.csv", colClasses = "character")[, 1])
educational_attainment_1a_pseudoids <- data.frame(ID = read.csv("educational_attainment_1a_results_na_labelled.csv", colClasses = "character")[, 1])
educational_attainment_2a_pseudoids <- data.frame(ID = read.csv("educational_attainment.csv", colClasses = "character")[, 1])
people_pseudoids <- data.frame(ID = read.csv("people.csv", colClasses = "character")[, 1])
rfft_pseudoids <- data.frame(ID = read.csv("rfft_results_na_labelled.csv", colClasses = "character")[, 1])
situation_pseudoids <- data.frame(ID = read.csv("situation.csv", colClasses = "character")[, 1])

datasets <- list(language_type_fluent_pseudoids, speech_imitation_pseudoids, speech_foreignaccent_pseudoids, speech_foreignlanguage_pseudoids,
age_of_acquisition_pseudoids, educational_attainment_1a_pseudoids, educational_attainment_2a_pseudoids, people_pseudoids, rfft_pseudoids, situation_pseudoids)
  
mmse_score_adu_m_01_3a_pseudoids <- data.frame(ID = read.csv("mmse_scores_filtered.csv", colClasses = "character")[, 1])
mmse_score_adu_m_01_3a_pseudoids <- mmse_score_adu_m_01_3a_pseudoids$ID

overlap_count <- sapply(datasets, function(df) sum(df$ID %in% mmse_score_adu_m_01_3a_pseudoids))
total_count <- sapply(datasets, nrow)
overlap_table <- data.frame(
  dataset = c("language_type_fluent_pseudoids", "speech_imitation_pseudoids", "speech_foreignaccent_pseudoids", "speech_foreignlanguage_pseudoids",
"age_of_acquisition_pseudoids", "educational_attainment_1a_pseudoids", "educational_attainment_2a_pseudoids", "people_pseudoids", "rfft_pseudoids", "situation_pseudoids"),
  all_pseudoids = total_count,
  mmse_pseudoid_overlap = overlap_count)

# Save it as a data table image. 
overlap_count_table <- tableGrob(overlap_table)
ggsave("overlap_count_table.png", plot = overlap_count_table, width = 8, height = 4, dpi = 300)


# ------------ Adjust/define the different language learning phenotypic traits based on their enumerations(Linux and R) ------------ # 
# This is under phenotypic_analysis_enumerations_2.sh script. 

# ------------ Make Distribution and Correlation Plots for the different traits (R) ------------ # 
# This is under phenotypic_analysis_correlations_3.sh script. 


# ------------ Merge the four main phenotype scores in a txt file for the GCTA analysis(R) ------------ # 

language_type_only_fluent <- fread("language_type_only_fluent_cogn_filtered.csv", stringsAsFactors = FALSE)
speech_imitation <- fread("speech_imitation_cogn_filtered.csv", stringsAsFactors = FALSE)
speech_foreignaccent <- fread("speech_foreignaccent_cogn_filtered.csv", stringsAsFactors = FALSE)
speech_foreignlanguage <- fread("speech_foreignlanguage_cogn_filtered.csv", stringsAsFactors = FALSE)

# Main phenotype scores + gender + age. 
main_pheno_scores_age_gender_cogn_filtered <- language_type_only_fluent %>%
  full_join(speech_imitation, by = "project_pseudo_id") %>%
  full_join(speech_foreignaccent, by = "project_pseudo_id") %>%
  full_join(speech_foreignlanguage, by = "project_pseudo_id") %>%
  select(project_pseudo_id, age.x, age.y, age.x.x, age.y.y, gender.x, gender.y, gender.x.x, gender.y.y,
		 n_fluent_languages, speech_imitation_adu_q_1, speech_foreignaccent_adu_q_1, speech_foreignlanguage_adu_q_1) %>%
  mutate(age_speech = coalesce(age.y, age.x.x, age.y.y)) %>% # For the age keep two seperate columns, one for the multilingualism age and one for the speech traits since they are from different questionairs. 
  mutate(gender = coalesce(gender.x, gender.y, gender.x.x, gender.y.y)) %>%
  select(project_pseudo_id, age_speech, age.x, gender, n_fluent_languages, speech_imitation_adu_q_1,
         speech_foreignaccent_adu_q_1, speech_foreignlanguage_adu_q_1) %>%
  rename(age_multilingualism = age.x, foreign_accents_understanding = speech_foreignaccent_adu_q_1, foreign_speech_imitation = speech_imitation_adu_q_1, foreign_language_learning = speech_foreignlanguage_adu_q_1)
  
fwrite(main_pheno_scores_age_gender_cogn_filtered, "main_pheno_scores_age_gender_cogn_filtered.txt", sep = "\t", na = "NA")
fwrite(main_pheno_scores_age_gender_cogn_filtered, "/groups/umcg-lifelines/tmp02/projects/ov21_0398/GCTA/frysk/heritability/input_files/main_pheno_scores_age_gender_cogn_filtered.txt", sep = "\t", na = "NA")
  

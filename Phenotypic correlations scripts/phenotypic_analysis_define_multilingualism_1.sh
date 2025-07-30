---
# author: "Katerina Spantidaki"
# title: "Phenotypic correlations"
# date: "2025-03-03"
---
# Script to get define the multilingualism phenotype of language learning phenotypic data.
# Mainly in R #
---

# Project directories. 
cd /groups/umcg-lifelines/tmp01/projects/ov21_0398/
cd /home/umcg-aspantidaki/phenotype_data_copy

# Load R.
module load R/4.2.2-foss-2022a-bare
R

# Load necessary libraries. 
library(dplyr)
library(data.table)
library(ggplot2)

# How to instal ggplot2 package:
#install.packages("pak")
#pak::pak("tidyverse/ggplot2")

# Make a dataframe with the variables in multilingualism that refer to the languges spoken by the individuals (In R). 
languages <- fread("multilingualism_results_na_labelled.csv", stringsAsFactors = FALSE)
languages[] <- lapply(languages, as.character)

languages <- languages %>%
  select(1:6, matches("language_type_adu_q_1_.*"), matches("language_typespec_adu_c_1_.*"))
type_columns <- grep("language_type_adu_q_1_.*", names(languages), value = TRUE)
typespec_columns <- grep("language_typespec_adu_c_1_.*", names(languages), value = TRUE)

# Replace the "10" in the language type columns with their respective language typespec column values.
# Appended "_new" to the typespec values used for this replacement.
# Loop through each type column and its corresponding typespec column
for (i in seq_along(type_columns)) {
  type_col <- type_columns[i]
  typespec_col <- typespec_columns[i]
  languages[[type_col]] <- ifelse(
    languages[[type_col]] == "10",
    paste0(languages[[typespec_col]], "_new"),
    languages[[type_col]]
  )
  languages[[typespec_col]] <- ifelse(
    languages[[type_col]] == paste0(languages[[typespec_col]], "_new"),
    paste0(languages[[typespec_col]], "_new"),
    languages[[typespec_col]]
  )
}
languages <- as.data.table(languages)
write.csv(languages, "language_type.csv", row.names = FALSE, quote = TRUE)

# Make a dataframe with the people (pseudoIDs) that answered in both type and typespec questions including their answers. 
results <- data.frame(Row = integer(), Column = character(), Value = character(), stringsAsFactors = FALSE)
for (col in typespec_columns) {
  invalid_rows <- which(!is.na(languages[[col]]) & !grepl("_new$", languages[[col]]))
  if (length(invalid_rows) > 0) {
    temp_results <- data.frame(
      Row = invalid_rows,
      Column = rep(col, length(invalid_rows)), # Use `rep()` to fill the column name
      Value = languages[invalid_rows, col, with = FALSE][[1]], # Extract the actual values
      stringsAsFactors = FALSE
    )
    results <- rbind(results, temp_results)
  }
}
column_mapping <- data.frame(
  Typespec = typespec_columns,
  Type = type_columns,
  stringsAsFactors = FALSE
)
comparison_results <- data.frame(PseudoID = character(), Typespec_Column = character(), 
                                 Type_Column = character(), Typespec_Value = character(), 
                                 Type_Value = character(), Match = logical(), 
                                 stringsAsFactors = FALSE)
for (i in 1:nrow(results)) {
  row_number <- results$Row[i]  
  typespec_column <- results$Column[i]  
  type_column <- column_mapping$Type[column_mapping$Typespec == typespec_column]
  pseudoid <- languages[row_number, 1][[1]]  
  typespec_value <- languages[row_number, ..typespec_column][[1]]
  type_value <- languages[row_number, ..type_column][[1]]
  comparison_results <- rbind(
    comparison_results,
    data.frame(
      PseudoID = pseudoid,
      Typespec_Column = typespec_column,
      Type_Column = type_column,
      Typespec_Value = typespec_value,
      Type_Value = type_value,
      Match = (typespec_value == type_value),  # Check if the values match
      stringsAsFactors = FALSE
    )
  )
}
write.csv(comparison_results, "answered_both_types_languages.csv", row.names = FALSE, quote = TRUE)

# There are 23 unique people (pseudoIDs) that answered in the same type and typespec question. 
# Exclude those from the language type dataframe.
cut -d',' -f1 answered_both_types_languages.csv | grep -Fv -f - language_type.csv > language_excluding_answered_both.csv	

language_excluding_answered_both <- fread("language_excluding_answered_both.csv", stringsAsFactors = FALSE)
language_type_filtered <- language_excluding_answered_both %>%
  select(-matches("language_typespec_adu_c_1_.*"))
	
# Replace the language numbers with the actual language type. 
language_types_enumerations <- fread("language_types_enumerations.csv", stringsAsFactors = FALSE)
for (i in 1:nrow(language_type_filtered)) {
  for (col in 7:ncol(language_type_filtered)) {  
    match_value <- language_type_filtered[i, ..col]  
    replacement <- language_types_enumerations[language_number == match_value, language_english]
    if (length(replacement) > 0) {
      language_type_filtered[i, (col) := replacement]  
    }
  }
}		 

# Remove the people (pseudoIDs) that have not answered (NA) any of the languages types. 
# 11043 --> 10987, 56 people were removed. 
cols_to_check <- 7:ncol(language_type_filtered)
na_row_numbers <- which(rowSums(is.na(language_type_filtered[, ..cols_to_check, with = FALSE])) == length(cols_to_check))
language_type_filtered <- language_type_filtered[-na_row_numbers, ]

write.csv(language_type_filtered, "language_type_not_with_numbers.csv", row.names = FALSE, quote = TRUE)	

# Make a dataframe for the language understanding ability variable and one for the language speaking ability.
# Keep the people (pseudoIDs) that have up to 6 languages (not the answered both ones). 
languages <- fread("multilingualism_results_na_labelled.csv", stringsAsFactors = FALSE)
languages[] <- lapply(languages, as.character)

language_understand <- languages %>%
  select(1:6, matches("language_understand_adu_q_1_.*")) %>%
  filter(project_pseudo_id %in% language_type_filtered$project_pseudo_id) 
language_speak <- languages %>%
  select(1:6, matches("language_speak_adu_q_1_.*")) %>%
  filter(project_pseudo_id %in% language_type_filtered$project_pseudo_id) 

# Merge the language type, understand and speak datasets into one. 
language_type_understand_speak <- language_type_filtered %>%
  inner_join(language_understand %>% select(project_pseudo_id, 7:ncol(language_understand)), by = "project_pseudo_id") %>%
  inner_join(language_speak %>% select(project_pseudo_id, 7:ncol(language_speak)), by = "project_pseudo_id")
  
write.csv(language_type_understand_speak, "language_type_understand_speak.csv", row.names = FALSE, quote = TRUE)	

# Replace all the scores below 7 as non-fluent and all the scores equal and above 7 as fluent. 
cols_to_update <- colnames(language_type_understand_speak[,13:22])
language_type_understand_speak <- language_type_understand_speak %>%
  mutate(across(all_of(cols_to_update), ~ ifelse(. < 7, "non-fluent", "fluent")))

# Add a fluent scores understand and speak columns for the native language type as well. 
# Find the people under language_type_adu_q_1_01 that have not answered to this question so that you can keep their scores also unanswered. 
na_row_numbers <- which(is.na(language_type_understand_speak$language_type_adu_q_1_01))
language_type_understand_speak <- language_type_understand_speak %>%
  mutate(
    language_understand_adu_q_1_01 = "fluent", 
    language_speak_adu_q_1_01 = "fluent"
  ) %>%
  mutate(
    language_understand_adu_q_1_01 = if_else(row_number() %in% na_row_numbers, NA_character_, language_understand_adu_q_1_01),
    language_speak_adu_q_1_01 = if_else(row_number() %in% na_row_numbers, NA_character_, language_speak_adu_q_1_01)
  ) %>%
  relocate(language_understand_adu_q_1_01, .before = "language_understand_adu_q_1_02") %>%
  relocate(language_speak_adu_q_1_01, .before = "language_speak_adu_q_1_02")

# For every person (pseudoID) if there is NA in language type but not NA in either of the language understand and language speak columns, replace those columns values with NAs.
# Do this for all 6 language type columns. 
# Remove the people (pseudoIDs) that have not answered (NA) any of the languages types. 
# From the 10987 people, 356 had their understand/speak scores replaced. 
understand_col <- paste0("language_understand_adu_q_1_0", 1:6)
speak_col <- paste0("language_speak_adu_q_1_0", 1:6)
type_col <- paste0("language_type_adu_q_1_0", 1:6)

replacements <- 0

for (i in 1:nrow(language_type_understand_speak)) {
  for (j in 1:6) {
    lt_col <- type_col[j]
    lu_col <- understand_col[j]
    ls_col <- speak_col[j]
    if (is.na(language_type_understand_speak[i, ..lt_col])) {
      if (!is.na(language_type_understand_speak[i, ..lu_col])) {
        set(language_type_understand_speak, i, lu_col, NA)
        replacements <- replacements + 1
      }
      if (!is.na(language_type_understand_speak[i, ..ls_col])) {
        set(language_type_understand_speak, i, ls_col, NA)
        replacements <- replacements + 1
      }
    }
  }
}

cat("Total replacements made:", replacements, "\n")

# For every person (pseudoID) keep only the language type that has a fluent score in both language understand and language speak columns. 
# Do this for all 6 languages. 
# 10987 --> 10947, 40 people were removed. 
language_type_only_fluent <- list()
for (i in 1:6) {
  understand_col <- paste0("language_understand_adu_q_1_0", i)
  speak_col <- paste0("language_speak_adu_q_1_0", i)
  type_col <- paste0("language_type_adu_q_1_0", i)
  language_type_fluent_df <- language_type_understand_speak %>%
    filter(!!sym(understand_col) == "fluent" & !!sym(speak_col) == "fluent") %>%
    select(project_pseudo_id, !!sym(type_col))
  colnames(language_type_fluent_df) <- c("project_pseudo_id", paste0("language_type_adu_q_1_0", i))
  language_type_only_fluent[[i]] <- language_type_fluent_df
}

# Merge the datasets for each of the 6 languages to one big dataset and add NAs in the empty languages per individual. 
language_type_only_fluent <- Reduce(function(x, y) merge(x, y, by = "project_pseudo_id", all = TRUE), language_type_only_fluent)

# Count the number of languages fluently spoken for each person and add them in a new column. 
language_type_only_fluent$n_fluent_languages <- rowSums(!is.na(language_type_only_fluent[, -1]))

# Add the additional information columns to the last dataset (language_type_only_fluent). 
language_type_only_fluent <- language_type_only_fluent %>%
  left_join(languages %>% select(1:6), by = "project_pseudo_id") %>%
  relocate(names(languages)[2:6], .after = project_pseudo_id)
write.csv(language_type_only_fluent, "language_type_only_fluent.csv", row.names = FALSE, quote = TRUE)

# Keep only the pseudoids that are not cognitively impaired.
# 10947 --> 8458, 2489 people were removed.
mmse_scores_filtered <- fread("mmse_scores_filtered.csv", stringsAsFactors = FALSE)
language_type_only_fluent_filtered <- language_type_only_fluent %>%
  filter(project_pseudo_id %in% mmse_scores_filtered$project_pseudo_id)
write.csv(language_type_only_fluent_filtered, "language_type_only_fluent_cogn_filtered.csv", row.names = FALSE, quote = TRUE)


# **README – Data analysis pipeline**

Study title: *Exploring the relationship between genetics and language acquisition.*

`Author: Katerina Spantidaki`

`Supervisors: Danielle Admiraal & Else Eising`

`Date: 28-07-2025`

——————————————————————————————————————————————

## Phenotypic correlations between the four main language acquisition traits and additional phenotypes.

We focus on four main traits: self-reported multilingualism, foreign speech mimicking ability, foreign language learning ability and foreign accents understanding ability. The latter three are ranked on a scale of 1-4, easy-difficult. For multilingualism, we collected more detailed data to support its definition, which we determined as the total number of languages an individual speaks fluently. Fluency was reported for speaking and understanding ability on a scale from 1 to 10. We considered an individual to be fluent in a language if they scored 7 or above on both measures. For the additional phenotypes we utilized these traits: educational attainment, type of situation of multiple language usage (home, work, social club, professional etc), type of people of multiple language usage (speaking with family, friends, colleagues etc.) and age of acquisition (AoA) of a language.

### [**Scripts & Steps Description**]{.smallcaps}

*\~\~\~\~\~\~\~ All scripts location: “/data/lag/workspaces/lg-lifelines/working_data/SCRIPTS/Phenotype_data/frysk"* *\~\~\~\~\~\~\~*

[Main script]{.underline}: *phenotypic_analysis.sh*

[Steps]{.underline}:

1.  Extract each trait, merge all the variables/questions that belong to it in a separate data frame and label the NAs.

2.  Define multilingualism (Script: *phenotypic_analysis_define_multilingualism_1.sh*):

    -   The language_type_adu_q_1_01-06 questions were asked and for each question, if the language the individual spoke was not present then they could chose another language from the language_typespec_adu_c_1_01-06 questions.

        –\> Remove the 23 people that answered in the same type and typespec question. We can only keep max 6 answers/languages per person since we have max 6 questions for the language_speak_adu_q_1_02-06 and the language_understand_adu_q_1_02-06 abilities.

    -   Replace all the scores below 7/10 in the language_speak_adu_q_1_02-06 and the language_understand_adu_q_1_02-06 abilities, as non-fluent and all the scores equal and above 7/10 as fluent.

    -   For every person keep only the language type that has a fluent score in both language understand and language speak columns.

    -   For each person count the number of languages fluently spoken.

    -   Keep only the people that are not cognitively impaired based on the mmse_score_adu_m_01 trait, next step, (8458/10947 --\> 2489 people were removed).

3.  Establish the non-cognitively impaired IDs trait:

    -   Because mmse_score_adu_m_01 exists in both studies, 3a_v_1_results.csv and 1a_v_1_results.csv, determine which one to keep. –\> mmse_score_adu_m_01_3a_results_na_labbeled.
    -   Filter for the non-cognitively impaired people (<https://wiki-lifelines.web.rug.nl/doku.php?id=mini_mental_state_examination&s>[[]=mmse&s[]=score&s[]=adu&s[]=01](https://wiki-lifelines.web.rug.nl/doku.php?id=mini_mental_state_examination&s%5B%5D=mmse&s%5B%5D=score&s%5B%5D=adu&s%5B%5D=01)):
        -   24-30 points: No cognitive impairment

        -   18-23 points: Mild to moderate cognitive impairment

        -   0-17 points: Severe cognitive impairment

4.  Adjust all the language learning phenotypic trait scores based on their enumerations (Script: *phenotypic_analysis_enumerations_2.sh*):

    -   For the age of acquisition (AoA) trait split it either in two groups: 1 –\> 0-2 yo & 2 –\> 3-100 yo or in four groups: 1 –\> 0-10 yo & 2 –\>11-20 yo & 3 –\> 21-40 yo & 4 –\> 41-100 yo. For the latter option make sure to include the people that answered more than one AoA (so more than one fluent languages), as many times as they gave an answer.

5.  Make distribution and correlation plots for the different traits (Script: *phenotypic_analysis_correlations_3.sh*):

    -   Keep only the individuals that are not cognitively impaired based on the mmse trait. Since the mmse questionnaire is available only for individuals aged 65 + and some of the traits include individuals across the full age range (0–100 years), apply this filter only to those aged 65+ in each trait, while keeping all individuals under 65 regardless of cognitive status.
    -   For the speech_imitation, speech_foreignaccent and speech_foreignlanguage traits flip the scores. So instead of 1 to 4 meaning Easy to Difficult, flip it so that it means Difficult to Easy.
    -   All the distribution and correlation plots are present here: *“/data/lag/workspaces/lg-lifelines/working_data/Phenotype_data/frysk/Phenotypes_plots"*.

——————————————————————————————————————————————

## Genetic analysis between the four main language acquisition traits.

### [**Scripts & Steps Description**]{.smallcaps}

1.  [**SNP-based heritability (SNP-h2) analysis**]{.underline} to determine weather the individual traits have a genetic liability. Measured SNP-h2 using the genome-based restricted maximum likelihood (GREML) analysis implemented in the genome-wide complex trait analysis (GCTA) (<https://yanglab.westlake.edu.cn/software/gcta/#GREMLanalysis>) software.

    *\~\~\~\~\~\~\~ All scripts location: “/data/lag/workspaces/lg-lifelines/working_data/SCRIPTS/*GCTA*/frysk/*SNP_h2*"* *\~\~\~\~\~\~\~*

    -   QC: MAF \> 0.01, HWE, pruning LD, r2 cut-off of 0.125

    -   We compared different genotyped data batches and evaluated the impact of including or excluding the cognitive impairment filter in our datasets to maximize the number of genotyped individuals.

        ![](images/SNPh2_bK02_language_acquisition_barplot-01.png)

    -   [Steps]{.underline}:

        1.  Prepare the genotype plink input files for the GRM construction GCTA (Script: *GCTA_prep_geno_gsa_UGLI2_Cyto_CHR.sh*).

        2.  Run the prep scripts of the plink input files for the GRM construction GCTA (Script: *GCTA_prep_geno_gsa_UGLI2_Cyto_make_and_run_all_chr.bash.sh*).

        3.  GRM generation per chromosome (to be submitted) (Script: *GCTA_calculate_GRM_CHR.sh*).

        4.  Run all GRM generation per chromosome scripts and create all GRMs (Script: *GCTA_calculate_GRM_make_and_run_all_chr.sh*).

        5.  Merge GRMs per chromosome, and create GRM based on IBS kernels (Script: *GCTA_calculate_GRM_merge.bK02.sh*).

        6.  Prepare pheno, geno and covar files for GCTA analysis (Script: *GCTA_prepare_covar_pheno_files_ugli_cyto_cogn_filtered.sh*).

        7.  Calculate SNP heritability using GCTA (a script per trait) (Script: *GCTA_calculate_SNP_h2_ugli_n_fluent_languages_cogn_filtered.sh*).

        8.  Plot SNP-h2 bK 0.02 (Script: *SNP_h2_plotting.sh*).

    -   All the SNP-h2 plots are present here: *“/data/lag/workspaces/lg-lifelines/working_data/GCTA/frysk/GCTA/heritability/Graphs"*.

2.  [**Genetic correlation (rG) analysis**]{.underline} to determine weather there are shared genetics across the traits. Measured rG using the genome-based restricted maximum likelihood (GREML) Bivariate analysis implemented in the genome-wide complex trait analysis (GCTA) (<https://yanglab.westlake.edu.cn/software/gcta/#BivariateGREMLanalysis>) software.

    *\~\~\~\~\~\~\~ All scripts location: “/data/lag/workspaces/lg-lifelines/working_data/SCRIPTS/*GCTA*/frysk/*gen_corr*"* *\~\~\~\~\~\~\~*

    -   QC: MAF \> 0.01, HWE, pruning LD, r2 cut-off of 0.125

    -   [Steps]{.underline}:

        1.  Prepare the input files for the pheno, geno and covar files for GCTA GREML Bivariate analysis (a script per combination of traits) (Script: *GCTA_bivar\_\_prepare_covar_pheno_files_ugli_cogn_filtered_LN.sh*).

        2.  Calculate rG using GCTA bivariate GREML (a script per combination of traits) (Script: *GCTA_GREML_genetic_correlations_BASE_LN.sh*).

            -   **LN**: foreign_language_learning vs n_fluent_language
            -   **IN**: foreign_language_learning vs foreign_accents_understanding
            -   **LI**: foreign_language_learning vs foreign_speech_imitation
            -   **UI**: foreign_accents_understanding vs foreign_speech_imitation
            -   **UL**: n_fluent_languages vs foreign_speech_imitation
            -   **UN**: n_fluent_languages vs foreign_accents_understanding

        3.  Plot rG (Script: *gen_corr_plotting.sh*).

    -   All the rG plots are present here: *“/data/lag/workspaces/lg-lifelines/working_data/GCTA/frysk/GCTA/gen_cor/Graphs"*.

3.  [**Polygenic score (PGS) analysis**]{.underline} to determine weather the trait outcomes can be predicted based on the trait-genetics. Measured PGS using the Polygenic Risk Scores via Continuous Shrinkage (PRS-CS) ([PRScs/README.md at master · getian107/PRScs · GitHub](https://github.com/getian107/PRScs/blob/master/README.md)) tool.

    *\~\~\~\~\~\~\~ All scripts location: “/data/lag/workspaces/lg-lifelines/working_data/SCRIPTS/PGS"* *\~\~\~\~\~\~\~*

    -   Reference paper: <https://pmc.ncbi.nlm.nih.gov/articles/PMC7612115/#S19>

    -   [Steps (]{.underline}<https://choishingwan.github.io/PRS-Tutorial/base/>) (Script: *PGS_steps.sh*):

        1.  `QC the base sumstats data to get the filtered SNP’s beta values` (schoolgrades.E4, <https://ipsych.dk/en/research/downloads/>) (*“/data/lag/workspaces/lg-lifelines/working_data/PGS/frysk/PGS_prep/sumstats"*) for the PGS analysis (Script: *Reformat_external_sumstats.sh*).

            | rsID | A1  | A2  | BETA   | SE     |
            |------|-----|-----|--------|--------|
            | rs…  | T   | A   | 0.0156 | 0.0201 |

        2.  `QC the target data to get the filtered SNP-list per language acquisition trait` for the PGS analysis. QC per SNPs and per individual (Script: *PGS_prep_geno_data_RUNALL.sh*).

            -   Merge UGLI1, UGLI2 and CYTO genotype platforms (Script: */data/lag/workspaces/lg-lifelines/working_data/SCRIPTS/*GCTA*/frysk/*SNP_h2/*GCTA_prep_geno_gsa_UGLI2_Cyto_CHR.sh*).

                -   Filter GSA, UGLI2 and CYTO genotype data for imputation quality and MAF.

                -   Convert them to PLINK best guess data.

                -   Merge GSA, UGLI2 and CYTO , keeping only overlapping SNPs

                -   Filter merged file on HWE 0.0001 and MAF 0.01.

            -   Filter out related individuals (Script: *PGS_individuals.extract.sh*).

            -   Filter on HWE 0.000001, MAF 0.01, SNPs with high missingness, individuals with high missingness and check and remove ambiguous SNPs (a script per chromosome) (Script: *PGS_prep_geno_data_filter.PLINK_CHR.sh*).

            -   Compile all chromosomal PLINK files into one (Script: *combine_chrom.sh*).

            -   Create PLINK files per trait with only unrelated individuals that have been genotyped and have pheno scores on that trait (Script: *plink_per.trait.sh*).

                | FID | IID | PAT | MAT | Sex | PHENO | rsID | rsID | rsID |
                |-----|-----|-----|-----|-----|-------|------|------|------|
                | 1   | …   | 0   | 0   | 0   | -9    | rs…  | rs…  | rs…  |

        3.  PGS_prep_geno_data_RUNALL.sh -\> `run PGS calculation analysis on the four main language acquisition traits.`

            -   Run PRScs - Adjust BETA’s (standardized posterior SNP effect sizes) based on LD structure genotype (Script: *PGS_RUN_base.sh*).

                | Chr | rsID | BPP | A1  | A2  | BETA_STD      |
                |-----|------|-----|-----|-----|---------------|
                | 1   | rs…  | …   | C   | A   | 1.062785e-03+ |

            -   Obtain individual-level polygenic scores by concatenating output files from all chromosomes with plink (Script: *PGS_plink_scores.sh*).

                | FID | IID | PHENO | CNT     | CNT2   | SCORE     |
                |-----|-----|-------|---------|--------|-----------|
                | 1   | …   | -9    | 1408796 | 325149 | 93405e-08 |

            -   Run regression between individual PGS scores and phenotype scores (Script: *PGS_steps.sh*).

            -   Plot PGS scores, pseudo R² and Odds Ratio (Script: *PGS_plotting.sh*).

    -   All the PGS analysis plots are present here: *“/data/lag/workspaces/lg-lifelines/working_data/PGS/frysk/PGS_plots"*.

——————————————————————————————————————————————

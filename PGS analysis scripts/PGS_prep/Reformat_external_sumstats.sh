#---
# author: Katerina Spantidaki.
# original author: Else Eising, 24 October 2023.
# title: "Script to prepare the external sumstats for PRScs".
# date: 14 May 2025.
# data: https://ipsych.dk/en/research/downloads/
#---

##-----------------------------------------------------
## Info on PRScs sources:
# Cleaning for PRScs: https://choishingwan.github.io/PRS-Tutorial/base/
# PRScs: see https://github.com/getian107/PRScs
# and PRS manual by Peyton Coleman: https://www.notion.so/PRS-Manual-Clapbeat-d644fc46e59b41229d178d5fb03fb73b
##-----------------------------------------------------

##-----------------------------------------------------
# Steps done in this script.
# STEP1:
# Totally 6391200 SNPs with MAF > 0.01 and INFO > 0.80 were included in the sumstats file. 
# To check (if possible): 
# N > Nmax/2.
# If chrX is present.
# Remove ambiguous SNPs (A/T or C/G).
# Remove duplicated SNPs.
# STEP2:
# Select and rename columns for PRScs
# Required columns: SNP, A1, A2, BETA, SE/p-value
# where SNP is the rs ID, A1 is the effect allele, A2 is the alternative allele, BETA is the effect of the A1 allele, SE is the standard error of the effect. 

##-----------------------------------------------------

##-----------------------------------------------------
# On gridmaster.
# Define paths.
##-----------------------------------------------------
WorkingDir="/data/workspaces/lag/workspaces/lg-lifelines/working_data/PGS/frysk/PGS_prep/sumstats/"
cd $WorkingDir/input/

##-----------------------------------------------------
## Sumstats: schoolgrades.E4.
## Obtained from here: https://ipsych.dk/en/research/downloads/
##-----------------------------------------------------
# Unzip the summary stats file.
gunzip schoolgrades.E4.sumstats.gz 

# Check header.
head -n 2 schoolgrades.E4.sumstats | column -t
#SNP         CHR  BP        A1  A2  BETA    SE      P       INFO      N
#rs10000000  4    40088896  T   A   0.0156  0.0201  0.4368  0.861565  30982

# Check for presence of X chromosome and for max and min N.
awk '$2 == "X"' schoolgrades.E4.sumstats | head # Not present
awk '{print $10}' schoolgrades.E4.sumstats | sort -n | head  # Nmin = 30982
awk '{print $10}' schoolgrades.E4.sumstats | sort -n | tail # Nmax = 30982

# Remove ambiguous and duplicated SNPs.
A1=4
A2=5
SNP_ID=1
awk -v A1="$A1" -v A2="$A2" -v SNP_ID="$SNP_ID" '{ if (($A1=="t" && $A2=="a") || ($A1=="a" && $A2=="t") || ($A1=="c" && $A2=="g") || ($A1=="g" && $A2=="c")) print $SNP_ID; else print $SNP_ID, "nonambig"; }' schoolgrades.E4.sumstats | grep -v nonambig > schoolgrades.E4._SNP.excl
wc -l schoolgrades.E4._SNP.excl # 0 SNPs are ambiguous.
awk -v SNP_ID="$SNP_ID" '{print $SNP_ID}' schoolgrades.E4.sumstats | sort | uniq -d > tmp1
wc -l tmp1 #22 duplicated rs ID's.
awk '{print $1}' schoolgrades.E4.sumstats | sort | uniq -c | awk '$1 > 1' | awk '{s+=$1} END {print s}' #these duplicated rs ID's are present 50197 times in these data. 
grep -vwf tmp1 schoolgrades.E4.sumstats > schoolgrades.E4.sumstats_cleaned.txt # 6341003 out of 6391200 SNPs remained.

# Obtain required columns for analyses, and adjust header name (SNP, A1, A2, BETA/OR, SE/p-value).
echo -e "rsID\tA1\tA2\tBETA\tSE" > schoolgrades.E4.sumstats_filtered.txt
awk '{ print $1, $4, $5, $6, $7 }' schoolgrades.E4.sumstats_cleaned.txt | tail -n +2 >> schoolgrades.E4.sumstats_filtered.txt

# Transform allele letters to uppercase.
awk '{ $2 = toupper($2); print }' schoolgrades.E4.sumstats_filtered.txt > schoolgrades.E4.sumstats_filtered1.txt
awk '{ $3 = toupper($3); print }' schoolgrades.E4.sumstats_filtered1.txt > schoolgrades.E4.sumstats_filtered2.txt

# Rename final output.
cp schoolgrades.E4.sumstats_filtered2.txt $WorkingDir/output/schoolgrades_E4_sumstats_PRScs_filtered.txt
# 6341003 out of 6391200 SNPs remained.

# Flip the effect sizes direction to have a positive correlation with the PGS of our traits. 
awk 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++) if($i=="BETA") col=i; print; next} {if(col){$col=-1*$col} print}' schoolgrades_E4_sumstats_PRScs_filtered.txt > schoolgrades_E4_sumstats_PRScs_filtered_flipped.txt

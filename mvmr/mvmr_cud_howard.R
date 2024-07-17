# The base script was created by Jasmine Khouja 27.04.22 and adapted by Chloe Burke April 2024.
# The script conducts multivariable MR exploring the effects of smoking initiation and cannabis use disorder on lifetime depression.
# The SNPs used in this script were selected based on the finding from GSCAN (Liu et al., 2019) and Levey et al., (2023)
# This script uses an alternative outcome GWAS by Howard et al., (2018), including UK Biobank

################################################################################
##### Load packages #####
################################################################################

library(usethis)
library(remotes)
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#devtools::install_github("mrcieu/ieugwasr") 
#install.packages("googleAuthR")
#ieugwasr::get_access_token() 
library(ieugwasr)
library(googleAuthR)
library(tidyverse)
library(stringr)
library(dplyr)
library(forestplot)
library(plyr)
library(gtable)
library(reshape)
library(gplots)
require(ggplot2)
library(ggplot2)
library(gridExtra)
library(grid)
library(extrafont)
library(plotly)
library(data.table)
library(curl)
#install.packages("MVMR")
#install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
library(MendelianRandomization)
#if (!require("devtools")) { install.packages("devtools") } else {}
#devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
library(simex)
library(writexl)
library(LDlinkR)


################################################################################
################################   MVMR  #######################################
################################################################################

######################################
##### Set SNP lists for SI + CUD  #### 
######################################
setwd(" ")
SI_SNPlist<- read.table("SI_SNPlist.txt", header=FALSE) 
SI_SNPlist<-rename(SI_SNPlist, c("SNP" = "V1"))
CUD_SNPlist<- read.table("CUD_SNPlist.txt", header=FALSE)
CUD_SNPlist<-rename(CUD_SNPlist, c("SNP" = "V1"))

########################################
##### Load full summary statistics #####
########################################

##############
##### SI #####
##############
setwd(" ")
si_dat <- read_exposure_data("SmokingInitiation.WithoutUKB.txt",
                             clump = FALSE,
                             sep = "\t",
                             snp_col = "RSID",
                             beta_col = "BETA",
                             se_col = "SE",
                             eaf_col = "AF",
                             effect_allele_col = "ALT",
                             other_allele_col = "REF",
                             pval_col = "PVALUE",
                             samplesize_col = "N",
                             min_pval = 1e-200,
                             log_pval = FALSE
)

si_dat_mr <- format_data(
  si_dat,
  type = "exposure",
  snps = SI_SNPlist$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

#######################################
#####   Clumping SI instrument    #####
#######################################
##Restrict the exposure data to independent SNPs only
si_dat_mr_clumped <- clump_data(si_dat_mr, clump_kb = 500, clump_r2 = 0.001) #316 SNPs

##################
## Saving files ##
##################
setwd(" ")
write.csv(si_dat_mr_clumped,"si_dat_mvmr.csv", row.names = FALSE)

###################
## Reading files ##
###################
#setwd(" ")
#si_dat_mr_clumped <- read.csv("si_dat_mvmr.csv", header = T)

###############
##### CUD #####
###############

# EAF column is missing as authors were unable to share this due to data sharing agreements
# Lead author confirmed that MAF was very similar to 1000G, so EAF are imputed for the analysis instruments using EAFs for the effect alleles available on: https://www.ensembl.org/index.html

setwd(" ")
cud_dat <- read_exposure_data("LOO_MVPmeta_iPSYCHremoved", 
                              clump = FALSE, 
                              sep = " ", 
                              snp_col = "rsid", 
                              beta_col = "logOR", 
                              se_col = "SE", 
                              effect_allele_col = "Allele1", 
                              other_allele_col = "Allele2", 
                              pval_col = "P", 
                              min_pval = 1e-200, 
                              log_pval = FALSE
)

cud_dat_mr <- format_data(
  cud_dat,
  type = "exposure",
  snps = CUD_SNPlist$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)


########################################
#####   Clumping CUD instrument    #####
########################################
##Restrict the exposure data to independent SNPs only
cud_dat_mr_clumped <- clump_data(cud_dat_mr, clump_kb = 500, clump_r2 = 0.001) #17 SNPs

##################
## Saving files ##
##################
setwd(" ")
write.csv(cud_dat_mr_clumped,"cud_dat_mvmr", row.names = FALSE)

###################
## Reading files ##
###################
#setwd(" ")
#cud_dat_mr_clumped <- read.csv("cud_dat_mvmr", header = T)

############################################
#####   Binding SI + CUD instruments    ####
############################################
SI_instr <- si_dat_mr_clumped[, c(
  "SNP",
  "effect_allele.exposure",
  "se.exposure",
  "pval.exposure",
  "beta.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure"
)]

CUD_instr <- cud_dat_mr_clumped[, c(
  "SNP",
  "effect_allele.exposure",
  "se.exposure",
  "pval.exposure",
  "beta.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure"
)]

INSTR_2 <- do.call("rbind", list(SI_instr, CUD_instr)) #333 SNPs

#######################################################################
##### Check for overlapping SNPs between the exposure instruments #####
#######################################################################
n_occur <- data.frame(table(INSTR_2$SNP))
n_occur[n_occur$Freq > 1, ]
INSTR_2[INSTR_2$SNP %in% n_occur$Var1[n_occur$Freq > 1], ] #0 overlapping SNPs

####################################################
##### Extract instruments from the SI data set #####
####################################################
si_dat_mvmr_2 <- format_data(
  si_dat,
  type = "exposure",
  snps = INSTR_2$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

##### Change name of GWAS and check n SNPs ####
si_dat_mvmr_2$id.exposure <- "1"
str(si_dat_mvmr_2) #missingSNPs = 0 (333/333)

#####################################################
##### Extract instruments from the CUD data set #####
#####################################################
cud_dat_mvmr <- format_data(
  cud_dat,
  type = "exposure",
  snps = INSTR_2$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  min_pval = 1e-200,
  log_pval = FALSE
)

# Adding in EAF from 1000G for SNPs in instrument (N = 333)
setwd(" ")
sicud_eaf <- read.csv("eaf1000g_mvmr_sicud.csv", header = T)
cud_dat_mvmr <- merge(cud_dat_mvmr, sicud_eaf, by = "SNP", all = TRUE)
cud_dat_mvmr$eaf.exposure <- cud_dat_mvmr$eaf_1000g
cud_dat_mvmr <- subset (cud_dat_mvmr, select = -eaf_1000g)

#Checking formatting of file
unique(cud_dat_mvmr$eaf.exposure)

#Removing trailing spaces
cud_dat_mvmr$eaf.exposure <- str_trim(cud_dat_mvmr$eaf.exposure)

#Convert to numeric
cud_dat_mvmr$eaf.exposure <- as.numeric(cud_dat_mvmr$eaf.exposure)

##### Change name of GWAS and check n SNPs ####
cud_dat_mvmr$id.exposure <- "2"
str(cud_dat_mvmr) #missingSNPs = 0 (333/333)

################################################################################
##### Merge them to extract from the other exposure as 'outcome' #####
################################################################################

si_dat_mvmr_3 <- subset(si_dat_mvmr_2, select = c(
  "SNP", "effect_allele.exposure",
  "other_allele.exposure",
  "beta.exposure", "se.exposure",
  "pval.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "eaf.exposure"
))

cud_dat_mvmr_1 <- subset(cud_dat_mvmr, select = c(
  "SNP",
  "effect_allele.exposure",
  "other_allele.exposure",
  "beta.exposure",
  "se.exposure",
  "pval.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "eaf.exposure"
))

##### check structure is the same #####
str(si_dat_mvmr_3)
str(cud_dat_mvmr_1)

################
##### Merge ####
################
si_cud_exp <- do.call("rbind", list(si_dat_mvmr_3, cud_dat_mvmr_1)) #666 SNPs

##########################
##### Add sample size ####
##########################

#Need something to differentiate exposure [1] and exposure [2] after binding to clump 
si_cud_exp$samplesize.exposure <- 0
si_cud_exp <- si_cud_exp %>%
  mutate(samplesize.exposure = case_when(
    id.exposure == 1 ~ samplesize.exposure + 249171,
    id.exposure == 2 ~ samplesize.exposure + 785635,
    TRUE ~ samplesize.exposure  
  ))

##### Save dataframe #####
setwd(" ")
write.csv(
  si_cud_exp,
  "SI_CUD_iPSYCH.csv",
  row.names = FALSE
)

##################################################################
##### Find proxies missing from either the SI or CUD data set #####
##################################################################
n_occur <- data.frame(table(si_cud_exp$SNP))
n_occur[n_occur$Freq < 2, ]
si_cud_exp[si_cud_exp$SNP %in% n_occur$Var1[n_occur$Freq < 2], ]

#0 missing SNPs

###################
##### Clumping ####
###################

##### Change all p-values for CUD to 1e-200 for clumping so that none are dropped #####
##### Save old p-values first #####
si_cud_exp$oldpvalues <- si_cud_exp$pval.exposure
si_cud_exp <- si_cud_exp %>%
  mutate(pval.exposure = if_else(
    si_cud_exp$SNP %in% CUD_instr$SNP,
    1e-201, pval.exposure
  ))

##### Clump the data #####
si_cud_exp$id.exposure[si_cud_exp$id.exposure == "2"] <- "1"
si_cud_exp <- clump_data(si_cud_exp, clump_kb = 500, clump_r2 = 0.001) 
#Removing 26 of 666 variants due to LD with other variants or absence from LD reference panel
str(si_cud_exp) #640 SNPs

##### Add ID's back #####
si_cud_exp$id.exposure[si_cud_exp$samplesize.exposure < 250000] <- "1"
si_cud_exp$id.exposure[si_cud_exp$samplesize.exposure > 250000] <- "2"

##### Revert all p-values for CUD from 1e-200 #####
si_cud_exp$pval.exposure <- si_cud_exp$oldpvalues
si_cud_exp <- select(si_cud_exp, -c(oldpvalues))

##### Split again to harmonise based on exposure id #####
SI_2 <- split(si_cud_exp, si_cud_exp$id.exposure)[["1"]]
CUD <- split(si_cud_exp, si_cud_exp$id.exposure)[["2"]]

##############################
##### Harmonise SI on CUD #####
##############################
names(SI_2) <- gsub("exposure", "outcome", names(SI_2))
SI_CUD <- harmonise_data(CUD, SI_2) 
#Removing the following SNPs for being palindromic with intermediate allele frequencies:
#rs10969352, rs1160685, rs13237637, rs1931431, rs3850736, rs4140932, rs4326350, rs56070621, rs6932350, rs7598402, rs7920501, rs7921378, rs986714

#####################################################################
##### Keep only snps that are present across both exposures ####
# Note: they would have frequency 1 if only available in one dataset
#####################################################################
n_occur <- data.frame(table(si_cud_exp$SNP))
n_occur[n_occur$Freq == 2, ]
si_cud_exp <- si_cud_exp[si_cud_exp$SNP %in% n_occur$Var1[n_occur$Freq == 2], ]
str(si_cud_exp) #640

############################
##### Format exposures #####
############################

##### Keep only SNPs where MrKeep = TRUE #####
SI_CUD <- SI_CUD[SI_CUD$mr_keep == TRUE, ]
str(SI_CUD) #307

##### Split the tables - CUD #####
CUD_H <- subset(
  SI_CUD,
  id.exposure == "2",
  select = c(
    SNP, exposure,
    id.exposure,
    effect_allele.exposure,
    other_allele.exposure,
    beta.exposure,
    se.exposure,
    pval.exposure,
    eaf.exposure
  )
)

##### Split the tables - SI #####
SI_H2 <- subset(SI_CUD,
                id.outcome == "1",
                select = c(
                  SNP,
                  outcome,
                  id.outcome,
                  effect_allele.outcome,
                  other_allele.outcome,
                  beta.outcome,
                  se.outcome,
                  pval.outcome,
                  eaf.outcome
                )
)

##### Turn SI from outcome to exposure to merge the data sets #####
names(SI_H2) <- gsub("outcome", "exposure", names(SI_H2))
SICUD_H <- merge(SI_H2, CUD_H, all = TRUE)
SICUD_H["Phenotype"] <- NA
SICUD_H$Phenotype[SICUD_H$id.exposure == 1] <- "SI"
SICUD_H$Phenotype[SICUD_H$id.exposure == 2] <- "CUD"
str(SICUD_H) #614 SNPs

#########################################
##### Extract outcome data for MVMR #####
#########################################

setwd(" ")

outcome_dat_sicud <- read_outcome_data(
  "howard_gwas_fixed.txt",
  snps = SICUD_H$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "OR",
  se_col = "SE",
  eaf_col = "FRQ_A_59851",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  min_pval = 1e-200,
  log_pval = FALSE
) #305

# x2 missing SNPs but this is a supplementary analysis so not going to proxy search for these

#################################
#####   Convert OR to log   #####
#################################
#Note. SE already for log odds
outcome_dat_sicud$beta.outcome<-log(outcome_dat_sicud$beta.outcome)

############################
##### Organise outcome #####
############################

outcome_dat_sicud["Phenotype"] <- NA
outcome_dat_sicud$Phenotype <- "Depression"

##################################
##### Harmonise with outcome #####
##################################

mvdat_sicud1 <- harmonise_data(SICUD_H, outcome_dat_sicud)
#Removing the following SNPs for being palindromic with intermediate allele frequencies: rs17554906

mvdat_sicud1 <- mvdat_sicud1[mvdat_sicud1$mr_keep == TRUE, ]
str(mvdat_sicud1) #SNPs = 608 

###########################
##### Save dataframes #####
###########################
setwd(" ")
write.csv(
  mvdat_sicud1,
  "mvdat_sicud_ipsych_howard.csv",
  row.names = FALSE
)

###########################
##### Reading in file #####
###########################
#setwd(" ")
#mvdat_sicud1 <- read.csv("mvdat_sicud_ipsych_howard.csv", header = T)

################################################################################
##### Run MVMR #####
################################################################################

##### IVW #####

bX1 <- c(mvdat_sicud1$beta.exposure[mvdat_sicud1$id.exposure == 1])
bX2 <- c(mvdat_sicud1$beta.exposure[mvdat_sicud1$id.exposure == 2])
bY <- c(mvdat_sicud1$beta.outcome[mvdat_sicud1$id.exposure == 1])
bYse <- c(mvdat_sicud1$se.outcome[mvdat_sicud1$id.exposure == 1])

set.seed(1234)
mod.MVMR_sicud <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(
  bY ~ bX1 + bX2 - 1,
  weights = bYse^-2
))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1,
                 weights = bYse^-2
  ))$sigma, 1)

mod_sicud <- summary(mod.MVMR_sicud)

mod_sicud_or <- coef(summary(mod.MVMR_sicud))
colnames(mod_sicud_or) <- c("b", "se", "t", "p")
mod_sicud_or <- as.data.frame(mod_sicud_or)
mod_sicud_or <- generate_odds_ratios(mod_sicud_or)

##### Orientation SI #####
##### As Egger analyses require the exposure betas to be positive,        #####
##### we first orient the betas to be positive for SI, and then           #####
##### orient the betas to be positive for CUD In the paper, we             #####
##### report the result for each exposure only with the right orientation #####

clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_sicud <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_sicud <- summary(mod.MVMRME_sicud)
mod_ME_sicud_or <- data.frame(mod.MVMRME_sicud[["coefficients"]])
colnames(mod_ME_sicud_or) <- c("b", "se", "t", "p")
mod_ME_sicud_or <- as.data.frame(mod_ME_sicud_or)
mod_ME_sicud_or <- generate_odds_ratios(mod_ME_sicud_or)

##### Orientation CUD #####
clist <- c("bX1", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX2>0,", var, ",", var, "*-1)")))
}
bX2 <- abs(bX2)

##### MVMR Egger #####
mod.MVMRME_sicud_2 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_sicud_2 <- summary(mod.MVMRME_sicud_2)

mod_ME_sicud_2_or <- data.frame(mod.MVMRME_sicud_2[["coefficients"]])
colnames(mod_ME_sicud_2_or) <- c("b", "se", "t", "p")
mod_ME_sicud_2_or <- as.data.frame(mod_ME_sicud_2_or)
mod_ME_sicud_2_or <- generate_odds_ratios(mod_ME_sicud_2_or)

##### Format to analyse / cross check with MVMR package #####
bX1 <- c(mvdat_sicud1$beta.exposure[mvdat_sicud1$id.exposure == 1])
bX2 <- c(mvdat_sicud1$beta.exposure[mvdat_sicud1$id.exposure == 2])
bY <- c(mvdat_sicud1$beta.outcome[mvdat_sicud1$id.exposure == 1])
bYse <- c(mvdat_sicud1$se.outcome[mvdat_sicud1$id.exposure == 1])
bXse1 <- c(mvdat_sicud1$se.exposure[mvdat_sicud1$id.exposure == 1])
bXse2 <- c(mvdat_sicud1$se.exposure[mvdat_sicud1$id.exposure == 2])
df_sicud <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

df_mvmr_sicud <- format_mvmr(
  df_sicud[, c(1, 3)],
  df_sicud[, 5],
  df_sicud[, c(2, 4)],
  df_sicud[, 6]
)

##### Cross check result with MVMR package #####
res_sicud <- ivw_mvmr(df_mvmr_sicud)

##################################################################################################
##### Calculate F-statistic and covariance #####
# Guidance here: https://wspiller.github.io/MVMR/articles/MVMR.html 
# Note: >10 is strong
# Where you have overlapping SNPs you can set the genetic covariance to 0:
# "If gene-exposure associations are estimated in seperate non-overlapping samples, then the covariances will be zero by design."
##### Robust estimates through Q-minimisation #####
# Following solution [2] to obtain pairwaise covariance i.e., correlation between the (phenotypic) exposures used by phenocov_mvmr() to approximate necessary covariance terms
# Correlation between SI and CUD in Levey (2023)  = 0.61 (Fig 1; doi: 10.1038/s41588-023-01563-z)
# At present, confidence intervals take a long time to generate, as such we report the OR compared to primary analysis
##################################################################################################
cov <- matrix(c(1, 0.61, 0.61, 1), nrow = 2, ncol = 2)

##### Calculate F-statistic #####
Xcovmat_sicud <- phenocov_mvmr(cov, df_mvmr_sicud[, c(6, 7)])
Fstat_sicud <- strength_mvmr(df_mvmr_sicud, Xcovmat_sicud)
#             exposure1 exposure2
#F-statistic   8.598132  4.158962

##### Calculate Cochran's Q #####
pres_sicud <- pleiotropy_mvmr(df_mvmr_sicud)
#Q = 501.4885  
#Qp = 3.217783e-12

##### Robust causal effect estimation #####
#"Where the MVMR assumptions are potentially violated, specifically where instruments are weak or exhibit pleiotropy...
#"it is possible to obtain more robust estimates through Q-statistic minimization. This can be performed using the qhet_mvmr() function."
cov_adj_mvmr_sicud <- qhet_mvmr(df_mvmr_sicud, cov, CI = FALSE, iterations = 1000)
# Converting into OR and LCI using exp(); first need to split 95% CI column by its separator (-), which needs to be adapted to consider negative values in 95%CI
cov_adj_mvmr_sicud2 <- separate(
  cov_adj_mvmr_sicud,
  col = "95% CI",
  into = c("LCI", "UCI"),
  sep = "(?<=\\d)-",
  remove = TRUE)

# Ensuring LCI and UCI are recognised as numeric values after splitting from original 95%CI column
cov_adj_mvmr_sicud2$LCI <- as.numeric(cov_adj_mvmr_sicud2$LCI)
cov_adj_mvmr_sicud2$UCI <- as.numeric(cov_adj_mvmr_sicud2$UCI)

# Converting values to OR and 95% CIs for tables/plots
cov_adj_mvmr_sicud_exp <- exp(cov_adj_mvmr_sicud2)

##### Calculate Cochran's Q #####
#With Q-minimisation
pres_sicud2 <- pleiotropy_mvmr(df_mvmr_sicud, gencov = Xcovmat_sicud)
#Q = 487.4888   
#Qp = 5.335574e-11
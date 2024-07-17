# The base script was created by Jasmine Khouja 27.04.22 and adapted by Chloe Burke April 2024.
# The script conducts multivariable MR exploring the effects of smoking initiation and cannabis initation on lifetime depression.
# The SNPs used in this script were selected based on the finding from GSCAN (Liu et al., 2019) and Pasman et al., (2018)

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
##### Set SNP lists for SI + CI  ##### 
######################################
setwd(" ")
SI_SNPlist<- read.table("SI_SNPlist.txt", header=FALSE) 
SI_SNPlist<-rename(SI_SNPlist, c("SNP" = "V1"))
CI_SNPlist<- read.table("CI_SNPlist.txt", header=FALSE)
CI_SNPlist<-rename(CI_SNPlist, c("SNP" = "V1"))

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

##############
##### CI #####
##############
setwd(" ")
ci_dat <- read_exposure_data("Cannabis_ICC_UKB_het.txt", 
                             clump = FALSE, 
                             sep = "\t", 
                             snp_col = "SNP", 
                             beta_col = "BETA", 
                             se_col = "SE", 
                             eaf_col = "FRQ",
                             effect_allele_col = "A1", 
                             other_allele_col = "A2", 
                             pval_col = "P", 
                             min_pval = 1e-200, 
                             log_pval = FALSE
)

ci_dat_mr <- format_data(
  ci_dat,
  type = "exposure",
  snps = CI_SNPlist$SNP,
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
#####   Clumping CI instrument    #####
#######################################
##Restrict the exposure data to independent SNPs only
ci_dat_mr_clumped <- clump_data(ci_dat_mr, clump_kb = 500, clump_r2 = 0.001) #6 SNPs

##################
## Saving files ##
##################
setwd(" ")
write.csv(ci_dat_mr_clumped,"ci_dat_mvmr.csv", row.names = FALSE)

###################
## Reading files ##
###################
#setwd(" ")
#ci_dat_mr_clumped <- read.csv("ci_dat_mvmr.csv", header = T)

############################################
#####   Binding SI + CI instruments    #####
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

CI_instr <- ci_dat_mr_clumped[, c(
  "SNP",
  "effect_allele.exposure",
  "se.exposure",
  "pval.exposure",
  "beta.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure"
)]

INSTR <- do.call("rbind", list(SI_instr, CI_instr))

#######################################################################
##### Check for overlapping SNPs between the exposure instruments #####
#######################################################################
n_occur <- data.frame(table(INSTR$SNP))
n_occur[n_occur$Freq > 1, ]
INSTR[INSTR$SNP %in% n_occur$Var1[n_occur$Freq > 1], ] #0 overlapping SNPs

####################################################
##### Extract instruments from the SI data set #####
####################################################
si_dat_mvmr <- format_data(
  si_dat,
  type = "exposure",
  snps = INSTR$SNP,
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
si_dat_mvmr$id.exposure <- "1"
str(si_dat_mvmr) #missingSNPs = 1 (321/322) = rs9773390
proxy_needed1 <- data.frame(setdiff(INSTR$SNP, si_dat_mvmr$SNP))

####################################################
##### Extract instruments from the CI data set #####
####################################################
ci_dat_mvmr <- format_data(
  ci_dat,
  type = "exposure",
  snps = INSTR$SNP,
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

##### Change name of GWAS and check n SNPs #####
ci_dat_mvmr$id.exposure <- "2"
str(ci_dat_mvmr) #missingSNPs = 1 (321/322) = rs2359180
proxy_needed2 <- data.frame(setdiff(INSTR$SNP, ci_dat_mvmr$SNP))

#########################################################################
##### Run proxy search and manually update SNP lists for extraction #####
#########################################################################
print(proxy_needed1) #rs9773390

#Proxy searching loop using LDLink
for (i in 1:nrow(proxy_needed1)) {
  x <- LDproxy(proxy_needed1[i , 1], pop = "EUR", r2d = "r2", token = "339b2b82e0ea")
  eligible <- x[x$R2 >= 0.8 , ]
  A <- si_dat %>%
    filter ( SNP %in% eligible$RS_Number)
  Common <- ci_dat %>%
    filter ( SNP %in% A$SNP)
  R2vals <- eligible %>%
    filter( RS_Number %in% Common$SNP)
  proxy_needed1$proxySNP[i] <- R2vals[1 , 1]
}

View(proxy_needed1) #now contains a column with a proxy SNP to use = rs73691962

print(proxy_needed2) #rs2359180

#Proxy searching loop using LDLink
for (i in 1:nrow(proxy_needed2)) {
  x <- LDproxy(proxy_needed2[i , 1], pop = "EUR", r2d = "r2", token = "339b2b82e0ea")
  eligible <- x[x$R2 >= 0.8 , ]
  A <- ci_dat %>%
    filter ( SNP %in% eligible$RS_Number)
  Common <- si_dat %>%
    filter ( SNP %in% A$SNP)
  R2vals <- eligible %>%
    filter( RS_Number %in% Common$SNP)
  proxy_needed2$proxySNP[i] <- R2vals[1 , 1]
}

View(proxy_needed2) # now contains a column with a proxy SNP to use = rs12607929

#Manually edit SNP lists swapping: rs2359180 in 'SI_SNPlist' for rs12607929, and rs9773390 in 'CI_SNPlist' for rs73691962

##########################################################################
##### Import new SNP lists and re-run initial format and clump steps #####
##########################################################################

###################################################
##### Set SNP lists with proxies for SI + CI  ##### 
###################################################
setwd(" ")
SI_SNPlist_proxy<- read.table("SI_SNPlist_Proxy.txt", header=FALSE) 
SI_SNPlist_proxy<-rename(SI_SNPlist_proxy, c("SNP" = "V1"))
CI_SNPlist_proxy<- read.table("CI_SNPlist_Proxy.txt", header=FALSE)
CI_SNPlist_proxy<-rename(CI_SNPlist_proxy, c("SNP" = "V1"))

#######################################
##### Extract SI SNPS with proxy  ##### 
#######################################

si_dat_mr <- format_data(
  si_dat,
  type = "exposure",
  snps = SI_SNPlist_proxy$SNP,
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
##### Extract CI SNPS with proxy  ##### 
#######################################

ci_dat_mr <- format_data(
  ci_dat,
  type = "exposure",
  snps = CI_SNPlist_proxy$SNP,
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

#############################################
##### Re-clump as using different SNPs  ##### 
#############################################
##Restrict the exposure data to independent SNPs only
si_dat_mr_proxy_clumped <- clump_data(si_dat_mr, clump_kb = 500, clump_r2 = 0.001) #316 SNPs

##Restrict the exposure data to independent SNPs only
ci_dat_mr_proxy_clumped <- clump_data(ci_dat_mr, clump_kb = 500, clump_r2 = 0.001) #6 SNPs

##################
## Saving files ##
##################
setwd(" ")
write.csv(si_dat_mr_proxy_clumped,"si_dat_mvmr_proxy.csv", row.names = FALSE)
write.csv(ci_dat_mr_proxy_clumped,"ci_dat_mvmr_proxy.csv", row.names = FALSE)

############################################
#####   Binding SI + CI instruments    #####
############################################
SI_instr <- si_dat_mr_proxy_clumped[, c(
  "SNP",
  "effect_allele.exposure",
  "se.exposure",
  "pval.exposure",
  "beta.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure"
)]

CI_instr <- ci_dat_mr_proxy_clumped[, c(
  "SNP",
  "effect_allele.exposure",
  "se.exposure",
  "pval.exposure",
  "beta.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure"
)]

INSTR <- do.call("rbind", list(SI_instr, CI_instr))

#######################################################################
##### Check for overlapping SNPs between the exposure instruments #####
#######################################################################
n_occur <- data.frame(table(INSTR$SNP))
n_occur[n_occur$Freq > 1, ]
INSTR[INSTR$SNP %in% n_occur$Var1[n_occur$Freq > 1], ] #0 overlapping SNPs

####################################################
##### Extract instruments from the SI data set #####
####################################################
si_dat_mvmr <- format_data(
  si_dat,
  type = "exposure",
  snps = INSTR$SNP,
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
si_dat_mvmr$id.exposure <- "1"
str(si_dat_mvmr) #missingSNPs = 0 (322/322)

####################################################
##### Extract instruments from the CI data set #####
####################################################
ci_dat_mvmr <- format_data(
  ci_dat,
  type = "exposure",
  snps = INSTR$SNP,
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

##### Change name of GWAS and check n SNPs #####
ci_dat_mvmr$id.exposure <- "2"
str(ci_dat_mvmr) #missingSNPs = 0 (322/322)

######################################################################
##### Merge them to extract from the other exposure as 'outcome' #####
######################################################################
si_dat_mvmr_1 <- subset(si_dat_mvmr, select = c(
  "SNP", "effect_allele.exposure",
  "other_allele.exposure",
  "beta.exposure", "se.exposure",
  "pval.exposure",
  "exposure",
  "mr_keep.exposure",
  "id.exposure",
  "eaf.exposure"
))

ci_dat_mvmr_1 <- subset(ci_dat_mvmr, select = c(
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
str(si_dat_mvmr_1)
str(ci_dat_mvmr_1)

################
##### Merge ####
################
si_ci_exp <- do.call("rbind", list(si_dat_mvmr_1, ci_dat_mvmr_1))

##########################
##### Add sample size ####
##########################
#Need something to differentiate exposure [1] and exposure [2] after binding to clump 
si_ci_exp$samplesize.exposure <- 0
si_ci_exp <- si_ci_exp %>%
  mutate(samplesize.exposure = case_when(
    id.exposure == 1 ~ samplesize.exposure + 249171,
    id.exposure == 2 ~ samplesize.exposure + 162082,
    TRUE ~ samplesize.exposure  # Handle other cases if needed
  ))

##### Save dataframe #####
setwd(" ")
write.csv(
  si_ci_exp,
  "SI_CI+PROXY.csv",
  row.names = FALSE
)

###################
##### Clumping ####
###################

##### Change all p-values for CI to 1e-200 for clumping so that none are dropped #####
##### Save old p-values first #####
si_ci_exp$oldpvalues <- si_ci_exp$pval.exposure
si_ci_exp <- si_ci_exp %>%
  mutate(pval.exposure = if_else(
    si_ci_exp$SNP %in% CI_instr$SNP,
    1e-201, pval.exposure
  ))

##### Clump the data #####
si_ci_exp$id.exposure[si_ci_exp$id.exposure == "2"] <- "1"
si_ci_exp <- clump_data(si_ci_exp, clump_kb = 500, clump_r2 = 0.001) 
#Removing 8 of 644 variants due to LD with other variants or absence from LD reference panel
str(si_ci_exp)

##### Add ID's back #####
si_ci_exp$id.exposure[si_ci_exp$samplesize.exposure > 240000] <- "1"
si_ci_exp$id.exposure[si_ci_exp$samplesize.exposure < 240000] <- "2"

##### Revert all p-values for CI from 1e-200 #####
si_ci_exp$pval.exposure <- si_ci_exp$oldpvalues
si_ci_exp <- select(si_ci_exp, -c(oldpvalues))

##### Split again to harmonise based on exposure id #####
SI <- split(si_ci_exp, si_ci_exp$id.exposure)[["1"]]
CI <- split(si_ci_exp, si_ci_exp$id.exposure)[["2"]]

##############################
##### Harmonise SI on CI #####
##############################
names(SI) <- gsub("exposure", "outcome", names(SI))
SI_CI <- harmonise_data(CI, SI)
#Removing the following SNPs for being palindromic with intermediate allele frequencies:
  #rs10969352, rs1160685, rs17554906, rs1931431, rs1937443, rs3850736, rs4140932, rs4326350, rs6932350, rs7598402, rs7920501, rs7921378, rs986714

#####################################################################
##### Keep only snps that are present across both exposures ####
# Note: they would have frequency 1 if only available in one dataset
#####################################################################
n_occur <- data.frame(table(si_ci_exp$SNP))
n_occur[n_occur$Freq == 2, ]
si_ci_exp <- si_ci_exp[si_ci_exp$SNP %in% n_occur$Var1[n_occur$Freq == 2], ]
str(si_ci_exp) #636

############################
##### Format exposures #####
############################

##### Keep onlySNPs where MrKeep = TRUE #####
SI_CI <- SI_CI[SI_CI$mr_keep == TRUE, ]
str(SI_CI) #305

##### Split the tables - CI #####
CI_H <- subset(
  SI_CI,
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
SI_H <- subset(SI_CI,
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
names(SI_H) <- gsub("outcome", "exposure", names(SI_H))
SICI_H <- merge(SI_H, CI_H, all = TRUE)
SICI_H["Phenotype"] <- NA
SICI_H$Phenotype[SICI_H$id.exposure == 1] <- "SI"
SICI_H$Phenotype[SICI_H$id.exposure == 2] <- "CI"
str(SICI_H)

#########################################
##### Extract outcome data for MVMR #####
#########################################

setwd(" ")
outcome_dat_sici <- read_outcome_data(
  "dep_full_imputed.txt",
  snps = SICI_H$SNP,
  sep = "\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM_INF",
  min_pval = 1e-200,
  log_pval = FALSE
) 
 
###########################################
##### Convert odds ratios to log odds #####
###########################################
outcome_dat_sici$beta.outcome <- as.numeric(
  as.character(
    outcome_dat_sici$beta.outcome
  )
)
outcome_dat_sici["beta.outcome"] <- (outcome_dat_sici$beta.outcome/(0.22652*((1-0.22652))))
outcome_dat_sici["se.outcome"] <- (outcome_dat_sici$se.outcome/(0.22652*((1-0.22652))))

############################
##### Organise outcome #####
############################
outcome_dat_sici["Phenotype"] <- NA
outcome_dat_sici$Phenotype <- "Depression"

##################################
##### Harmonise with outcome #####
##################################
mvdat_sici1 <- harmonise_data(SICI_H, outcome_dat_sici)
mvdat_sici1 <- mvdat_sici1[mvdat_sici1$mr_keep == TRUE, ]
str(mvdat_sici1)

###########################
##### Save dataframes #####
###########################
setwd(" ")
write.csv(
  mvdat_sici1,
  "mvdat_sici_mdd.csv",
  row.names = FALSE
)

###########################
##### Read dataframes #####
###########################
#setwd(" ")
#mvdat_sici1 <- read.csv ("mvdat_sici_mdd.csv", header = T)

################################################################################
##### Run MVMR #####
################################################################################

##### IVW #####

bX1 <- c(mvdat_sici1$beta.exposure[mvdat_sici1$id.exposure == 1])
bX2 <- c(mvdat_sici1$beta.exposure[mvdat_sici1$id.exposure == 2])
bY <- c(mvdat_sici1$beta.outcome[mvdat_sici1$id.exposure == 1])
bYse <- c(mvdat_sici1$se.outcome[mvdat_sici1$id.exposure == 1])

set.seed(1234)
mod.MVMR_sici <- lm(bY ~ bX1 + bX2 - 1, weights = bYse^-2)
se_theta1MI.random <- summary(lm(
  bY ~ bX1 + bX2 - 1,
  weights = bYse^-2
))$coef[1, 2] /
  min(summary(lm(bY ~ bX1 + bX2 - 1,
                 weights = bYse^-2
  ))$sigma, 1)

mod_sici <- summary(mod.MVMR_sici)

mod_sici_or <- coef(summary(mod.MVMR_sici))
colnames(mod_sici_or) <- c("b", "se", "t", "p")
mod_sici_or <- as.data.frame(mod_sici_or)
mod_sici_or <- generate_odds_ratios(mod_sici_or)

##### Orientation SI #####
##### As Egger analyses require the exposure betas to be positive,        #####
##### we first orient the betas to be positive for SI, and then           #####
##### orient the betas to be positive for CI In the paper, we             #####
##### report the result for each exposure only with the right orientation #####

clist <- c("bX2", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX1>0,", var, ",", var, "*-1)")))
}
bX1 <- abs(bX1)

##### MVMR Egger #####
mod.MVMRME_sici <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(
  bY ~ bX1 + bX2,
  weights = bYse^-2
))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_sici <- summary(mod.MVMRME_sici)
mod_ME_sici_or <- data.frame(mod.MVMRME_sici[["coefficients"]])
colnames(mod_ME_sici_or) <- c("b", "se", "t", "p")
mod_ME_sici_or <- as.data.frame(mod_ME_sici_or)
mod_ME_sici_or <- generate_odds_ratios(mod_ME_sici_or)

##### Orientation CI #####
clist <- c("bX1", "bY")
for (var in clist) {
  eval(parse(text = paste0(var, "<-ifelse(bX2>0,", var, ",", var, "*-1)")))
}
bX2 <- abs(bX2)

##### MVMR Egger #####
mod.MVMRME_sici_2 <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))
se_theta1ME.random <- summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$coef[2, 2] /
  min(summary(lm(bY ~ bX1 + bX2, weights = bYse^-2))$sigma, 1)
mod_ME_sici_2 <- summary(mod.MVMRME_sici_2)

mod_ME_sici_2_or <- data.frame(mod.MVMRME_sici_2[["coefficients"]])
colnames(mod_ME_sici_2_or) <- c("b", "se", "t", "p")
mod_ME_sici_2_or <- as.data.frame(mod_ME_sici_2_or)
mod_ME_sici_2_or <- generate_odds_ratios(mod_ME_sici_2_or)

##### Format to analyse / cross check with MVMR package #####
bX1 <- c(mvdat_sici1$beta.exposure[mvdat_sici1$id.exposure == 1])
bX2 <- c(mvdat_sici1$beta.exposure[mvdat_sici1$id.exposure == 2])
bY <- c(mvdat_sici1$beta.outcome[mvdat_sici1$id.exposure == 1])
bYse <- c(mvdat_sici1$se.outcome[mvdat_sici1$id.exposure == 1])
bXse1 <- c(mvdat_sici1$se.exposure[mvdat_sici1$id.exposure == 1])
bXse2 <- c(mvdat_sici1$se.exposure[mvdat_sici1$id.exposure == 2])
df_sici <- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

df_mvmr_sici <- format_mvmr(
  df_sici[, c(1, 3)],
  df_sici[, 5],
  df_sici[, c(2, 4)],
  df_sici[, 6]
)

##### Cross check result with MVMR package #####
res_sici <- ivw_mvmr(df_mvmr_sici)

##################################################################################################
##### Calculate F-statistic and covariance #####
# Guidance here: https://wspiller.github.io/MVMR/articles/MVMR.html 
# Note: >10 is strong
# Where you have overlapping SNPs you can set the genetic covariance to 0:
# "If gene-exposure associations are estimated in seperate non-overlapping samples, then the covariances will be zero by design."
# This is not the case for your analysis 
##### Robust estimates through Q-minimisation #####
# Following solution [2] to obtain pairwaise covariance 
# i.e., correlation between the (phenotypic) exposures used by phenocov_mvmr() to approximate necessary covariance terms
# Correlation between SI and CI in GSCAN  = 0.60 (Fig 1; doi: 10.1038/s41588-018-0307-5)
# At present, confidence intervals take a long time to generate
##################################################################################################
cov <- matrix(c(1, 0.60, 0.60, 1), nrow = 2, ncol = 2)

##### Calculate F-statistic #####
Xcovmat_sici <- phenocov_mvmr(cov, df_mvmr_sici[, c(6, 7)])
Fstat_sici <- strength_mvmr(df_mvmr_sici, gencov = Xcovmat_sici)
#              exposure1    exposure2
#F-statistic   9.137364     2.849941

##### Calculate Cochran's Q #####
pres_sici <- pleiotropy_mvmr(df_mvmr_sici)
#Q = 710.8835
#Qp = 5.352388e-35

##### Robust causal effect estimation #####
#"Where the MVMR assumptions are potentially violated, specifically where instruments are weak or exhibit pleiotropy...
#"it is possible to obtain more robust estimates through Q-statistic minimization. This can be performed using the qhet_mvmr() function."
cov_adj_mvmr_sici <- qhet_mvmr(df_mvmr_sici, cov, CI = TRUE, iterations = 1000)
# Converting into OR and LCI using exp(); first need to split 95% CI column by its separator (-), which needs to be adapted to consider negative values in 95%CI
cov_adj_mvmr_sici2 <- separate(
  cov_adj_mvmr_sici,
  col = "95% CI",
  into = c("LCI", "UCI"),
  sep = "(?<=\\d)-",
  remove = TRUE)

# Ensuring LCI and UCI are recognised as numeric values after splitting from original 95%CI column
cov_adj_mvmr_sici2$LCI <- as.numeric(cov_adj_mvmr_sici2$LCI)
cov_adj_mvmr_sici2$UCI <- as.numeric(cov_adj_mvmr_sici2$UCI)

# Converting values to OR and 95% CIs for tables/plots
cov_adj_mvmr_sici_exp <- exp(cov_adj_mvmr_sici2)

##### Calculate Cochran's Q #####
#With Q-minimisation
pres_sici2 <- pleiotropy_mvmr(df_mvmr_sici, gencov = Xcovmat_sici)
#Q = 721.0061 
#Qp = 2.800037e-36


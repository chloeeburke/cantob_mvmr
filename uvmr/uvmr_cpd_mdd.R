# Script was created by Chloe Burke January 2024.
# The script conducts univariable MR exploring the effects of smoking heaviness on lifetime depression.
# The SNPs used in this script were selected based on the finding from GSCAN (Liu et al., 2019)

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
library(simex)
library(writexl)
library(LDlinkR)

################################
##### Set SNP list for CPD ##### 
################################
setwd(" ")
CPD_SNPlist<- read.table("CPD_SNPlist.txt", header=FALSE) 
CPD_SNPlist<-rename(CPD_SNPlist, c("SNP" = "V1"))

##############################################################
##### Extract exposure data for MR of CPD and depression #####
##############################################################
setwd(" ")
cpd_dat <- read_exposure_data("CigarettesPerDay.WithoutUKB.txt",
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

cpd_dat_mr <- format_data(
  cpd_dat,
  type = "exposure",
  snps = CPD_SNPlist$SNP,
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
#####   Clumping CPD instrument    #####
#######################################
##Restrict the exposure data to independent SNPs only
cpd_dat_mr_clumped <- clump_data(cpd_dat_mr, clump_kb = 500, clump_r2 = 0.001) #41 SNPs

################################################################################
##### Extract outcome data for MR #####
################################################################################
#Set working directory to location of depression summary stats
setwd(" ")

#########
## CPD ##
#########
#Never smokers
mdd_dat_never_cpd_clumped <- read_outcome_data(
  "dep_never_imputed.txt",
  snps = cpd_dat_mr_clumped$SNP,
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

#Ever smokers
mdd_dat_ever_cpd_clumped <- read_outcome_data(
  "dep_ever_imputed.txt",
  snps = cpd_dat_mr_clumped$SNP,
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

##################
## Saving files ##
##################
setwd(" ")
write.csv(cpd_dat_mr_clumped,"cpd_dat_uvmr.csv", row.names = FALSE)
write.csv(mdd_dat_never_cpd_clumped,"cpd_mdd_uvmr_never.csv", row.names = FALSE)
write.csv(mdd_dat_ever_cpd_clumped,"cpd_mdd_uvmr_never_ever.csv", row.names = FALSE)

###################################################
##### Convert odds ratios to log odds #####
# Note: the exposures are already on the log scale.
###################################################
# To convert binary outcomes from BOLT-LMM files to log odds, the following calculation needs to be made: logOR<-beta(caseprevelance*(1-caseprevelance))
# Case prevalence is found in the log file by searching Phenotype and looking at the mean value.
# As the phenotype was coded as 1 or 2, the case prevalence (%) is the mean -1*100
# Standard errors of SNP effect size estimates should also be divided by (case prevalence * (1 - case prevalence)) when applying the above transformation
# Mean value MDD never smokers = 1.20349, case prevalence = 20.349% of 194881: caseN = 39656, ctrlN = 155225
# Mean value MDD ever smokers = 1.25484, case prevalence = 25.484% of 160248: caseN = 40838, ctrlN = 119410

mdd_dat_never_cpd_clumped$beta.outcome <- as.numeric(as.character(
  mdd_dat_never_cpd_clumped$beta.outcome
))

mdd_dat_ever_cpd_clumped$beta.outcome <- as.numeric(as.character(
  mdd_dat_ever_cpd_clumped$beta.outcome
))

mdd_dat_never_cpd_clumped["beta.outcome"] <- (mdd_dat_never_cpd_clumped$beta.outcome/(0.20349*((1-0.20349))))
mdd_dat_never_cpd_clumped["se.outcome"] <- (mdd_dat_never_cpd_clumped$se.outcome/(0.20349*((0.20349))))

mdd_dat_ever_cpd_clumped["beta.outcome"] <- (mdd_dat_ever_cpd_clumped$beta.outcome/(0.25484*((1-0.25484))))
mdd_dat_ever_cpd_clumped["se.outcome"] <- (mdd_dat_ever_cpd_clumped$se.outcome/(0.25484*((1-0.25484))))

################################################################################
##### Harmonising #####
################################################################################

##CPD - rs1737894 palindromic with intermediate allele frequencies
outcome_mdd_dat_cpd_never_clumped <-harmonise_data(cpd_dat_mr_clumped, mdd_dat_never_cpd_clumped, action = 2)
outcome_mdd_dat_cpd_ever_clumped <-harmonise_data(cpd_dat_mr_clumped, mdd_dat_ever_cpd_clumped, action = 2)

##################
## Saving files ##
##################
setwd(" ")
write.csv(outcome_mdd_dat_cpd_never_clumped,"outcome_mdd_cpd_uvmr_never.csv", row.names = FALSE)
write.csv(outcome_mdd_dat_cpd_ever_clumped,"outcome_mdd_cpd_uvmr_ever.csv", row.names = FALSE)

###################
## Reading files ##
###################
#setwd(" ")
#outcome_mdd_dat_cpd_never_clumped <- read.csv("outcome_mdd_cpd_uvmr_never.csv", header = T)
#outcome_mdd_dat_cpd_ever_clumped <- read.csv("outcome_mdd_cpd_uvmr_ever.csv", header = T)

################################################################################
################################## MR ##########################################
################################################################################

# CPD never smokers
result_cpd_never_clumped <- mr(
  outcome_mdd_dat_cpd_never_clumped,
  method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
)
result_cpd_never_clumped <- generate_odds_ratios(result_cpd_never_clumped)
outcome_mdd_dat_cpd_never_clumped <- subset(outcome_mdd_dat_cpd_never_clumped, mr_keep)
mr_n_cpd_clmp <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_cpd_never_clumped$beta.exposure,
  outcome_mdd_dat_cpd_never_clumped$beta.outcome,
  outcome_mdd_dat_cpd_never_clumped$se.exposure,
  outcome_mdd_dat_cpd_never_clumped$se.outcome,
  parameters = default_parameters()
)
ptr_n_cpd_clmp <- data.frame(mr_n_cpd_clmp["Q"])
egger_n_cpd_clmp <- mr_egger_regression(
  outcome_mdd_dat_cpd_never_clumped$beta.exposure,
  outcome_mdd_dat_cpd_never_clumped$beta.outcome,
  outcome_mdd_dat_cpd_never_clumped$se.exposure,
  outcome_mdd_dat_cpd_never_clumped$se.outcome,
  parameters
)
ptr_n_cpd_clmp[2, 1] <- egger_n_cpd_clmp["Q"]
F <- abs(outcome_mdd_dat_cpd_never_clumped$beta.exposure)^2 / outcome_mdd_dat_cpd_never_clumped$se.exposure^2
mF <- mean(F)
ptr_n_cpd_clmp[3, 1] <- mF

cpd_n_mdd_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                             SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                             OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = outcome_mdd_dat_cpd_never_clumped, NbDistribution = 1000,  SignifThreshold = 0.05)

cpd_n_mdd_presso_mainres <- cpd_n_mdd_presso$`Main MR results`
cpd_n_mdd_presso_mainres <- cpd_n_mdd_presso_mainres %>%
  rename(b = `Causal Estimate`)
cpd_n_mdd_presso_mainres <- cpd_n_mdd_presso_mainres %>%
  rename(se = `Sd`)
cpd_n_mdd_presso_mainres <- generate_odds_ratios(cpd_n_mdd_presso_mainres)
cpd_n_mdd_presso_mainres$Distortion_P <- cpd_n_mdd_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
cpd_n_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT <- ifelse(cpd_n_mdd_presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, TRUE, FALSE)
summary(cpd_n_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT) #No outliers identified
ptr_n_cpd_clmp[4, 1] <- "NA"

# CPD ever smokers
result_cpd_ever_clumped <- mr(
  outcome_mdd_dat_cpd_ever_clumped,
  method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
)
result_cpd_ever_clumped <- generate_odds_ratios(result_cpd_ever_clumped)
outcome_mdd_dat_cpd_ever_clumped <- subset(outcome_mdd_dat_cpd_ever_clumped, mr_keep)
mr_e_cpd_clmp <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_cpd_ever_clumped$beta.exposure,
  outcome_mdd_dat_cpd_ever_clumped$beta.outcome,
  outcome_mdd_dat_cpd_ever_clumped$se.exposure,
  outcome_mdd_dat_cpd_ever_clumped$se.outcome,
  parameters = default_parameters()
)
ptr_e_cpd_clmp <- data.frame(mr_e_cpd_clmp["Q"])
egger_e_cpd_clmp <- mr_egger_regression(
  outcome_mdd_dat_cpd_ever_clumped$beta.exposure,
  outcome_mdd_dat_cpd_ever_clumped$beta.outcome,
  outcome_mdd_dat_cpd_ever_clumped$se.exposure,
  outcome_mdd_dat_cpd_ever_clumped$se.outcome,
  parameters
)
ptr_e_cpd_clmp[2, 1] <- egger_e_cpd_clmp["Q"]
F <- abs(outcome_mdd_dat_cpd_ever_clumped$beta.exposure)^2 / outcome_mdd_dat_cpd_ever_clumped$se.exposure^2
mF <- mean(F)
ptr_e_cpd_clmp[3, 1] <- mF

cpd_e_mdd_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                             SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                             OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = outcome_mdd_dat_cpd_ever_clumped, NbDistribution = 1000,  SignifThreshold = 0.05)

cpd_e_mdd_presso_mainres <- cpd_e_mdd_presso$`Main MR results`
cpd_e_mdd_presso_mainres <- cpd_e_mdd_presso_mainres %>%
  rename(b = `Causal Estimate`)
cpd_e_mdd_presso_mainres <- cpd_e_mdd_presso_mainres %>%
  rename(se = `Sd`)
cpd_e_mdd_presso_mainres <- generate_odds_ratios(cpd_e_mdd_presso_mainres)
cpd_e_mdd_presso_mainres$Distortion_P <- cpd_e_mdd_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
cpd_e_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT <- ifelse(cpd_e_mdd_presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, TRUE, FALSE)
summary(cpd_e_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT) #nOUT TRUE = 0 SNPs  
ptr_e_cpd_clmp[4, 1] <- "NA"

setwd(" ")
write.csv(cpd_n_mdd_presso_mainres, "presso_cpd_never.csv", row.names=FALSE, quote=FALSE)
write.csv(cpd_e_mdd_presso_mainres, "presso_cpd_ever.csv", row.names=FALSE, quote=FALSE)

################################################################################
##### Calculate I2GX #####
# Isq(y, s) where y = vector of effects and s = vector of standard errors
# Apply simulation extrapolation SIMEX corrections to MR-Egger analysis where
# I2GX estimates < 0.9
# This indicates the effect estimate is biased by 10% due to measurement error
# (Bowden, Del Greco, et al., 2016)
# (Lederer & K?chenhoff, 2006)
################################################################################

##Manual ISQ function:
ISQ <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  ISQ        = (Q - (k-1))/Q
  ISQ        = max(0,ISQ)
  return(ISQ)
}
##Never smokers
#Calculate the regression dilution Isq. 
BetaXG   = outcome_mdd_dat_cpd_never_clumped$beta.exposure
seBetaXG = outcome_mdd_dat_cpd_never_clumped$se.exposure 
seBetaYG <- outcome_mdd_dat_cpd_never_clumped$se.outcome

BXG  = abs(BetaXG)         # gene--exposure estimates are positive  

Isq_unweighted_cpd_n <- ISQ(BXG,seBetaXG) #unweighted = 0.91
Isq_weighted_cpd_n <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.91

##Ever smokers
BetaXG   = outcome_mdd_dat_cpd_ever_clumped$beta.exposure
seBetaXG = outcome_mdd_dat_cpd_ever_clumped$se.exposure 
seBetaYG <- outcome_mdd_dat_cpd_ever_clumped$se.outcome

BXG  = abs(BetaXG)         # gene--exposure estimates are positive  

Isq_unweighted_cpd_e <- ISQ(BXG,seBetaXG) #unweighted = 0.91
Isq_weighted_cpd_e <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.91

#################################################
#########   Binding results tables   ############
#################################################
##Never smokers
#Add other results into main results data frame (i.e., PRESSO)
#PRESSO:
#No outliers identified so this is the "raw" estimate (i.e., not outlier corrected)
result_cpd_never_clumped[5, 1] <- result_cpd_never_clumped[1, "id.exposure"]
result_cpd_never_clumped[5,2] <- result_cpd_never_clumped[1, "id.outcome"]
result_cpd_never_clumped[5,3] <- result_cpd_never_clumped[1, "outcome"]
result_cpd_never_clumped[5,4] <- result_cpd_never_clumped[1, "exposure"]
result_cpd_never_clumped[5,5] <- "MR PRESSO"
result_cpd_never_clumped[5,6] <- result_cpd_never_clumped[1, "nsnp"] #No outliers removed so same as other analyses
result_cpd_never_clumped[5,7] <- cpd_n_mdd_presso_mainres[1, "b"]
result_cpd_never_clumped[5,8] <- cpd_n_mdd_presso_mainres[1, "se"]
result_cpd_never_clumped[5,9] <- cpd_n_mdd_presso_mainres[1, "P-value"]
result_cpd_never_clumped[5,10] <- cpd_n_mdd_presso_mainres[1, "lo_ci"]
result_cpd_never_clumped[5,11] <- cpd_n_mdd_presso_mainres[1, "up_ci"]
result_cpd_never_clumped[5,12] <- cpd_n_mdd_presso_mainres[1, "or"]
result_cpd_never_clumped[5,13] <- cpd_n_mdd_presso_mainres[1, "or_lci95"]
result_cpd_never_clumped[5,14] <- cpd_n_mdd_presso_mainres[1, "or_uci95"]

##Ever smokers
#Add other results into main results data frame (i.e., PRESSO)
#PRESSO:
#Outliers not identified so this is the raw estimate 
result_cpd_ever_clumped[5, 1] <- result_cpd_ever_clumped[1, "id.exposure"]
result_cpd_ever_clumped[5,2] <- result_cpd_ever_clumped[1, "id.outcome"]
result_cpd_ever_clumped[5,3] <- result_cpd_ever_clumped[1, "outcome"]
result_cpd_ever_clumped[5,4] <- result_cpd_ever_clumped[1, "exposure"]
result_cpd_ever_clumped[5,5] <- "MR PRESSO"
result_cpd_ever_clumped[5,6] <- result_cpd_ever_clumped[1, "nsnp"] #No outliers removed so same as other analyses
result_cpd_ever_clumped[5,7] <- cpd_e_mdd_presso_mainres[1, "b"]
result_cpd_ever_clumped[5,8] <- cpd_e_mdd_presso_mainres[1, "se"]
result_cpd_ever_clumped[5,9] <- cpd_e_mdd_presso_mainres[1, "P-value"]
result_cpd_ever_clumped[5,10] <- cpd_e_mdd_presso_mainres[1, "lo_ci"]
result_cpd_ever_clumped[5,11] <- cpd_e_mdd_presso_mainres[1, "up_ci"]
result_cpd_ever_clumped[5,12] <- cpd_e_mdd_presso_mainres[1, "or"]
result_cpd_ever_clumped[5,13] <- cpd_e_mdd_presso_mainres[1, "or_lci95"]
result_cpd_ever_clumped[5,14] <- cpd_e_mdd_presso_mainres[1, "or_uci95"]

setwd(" ")
write.csv(result_cpd_never_clumped, "results_cpd_never.csv", row.names=FALSE, quote=FALSE)
write.csv(result_cpd_ever_clumped, "results_cpd_ever.csv", row.names=FALSE, quote=FALSE)

#####################################################
#########   Test for reverse causation   ############
#####################################################

### Steiger filtering ###
# Binary exposures/outcomes must contain: samplesize, ncase, ncontrol, prevalence and units="log odds"
# For smoking heaviness there were were 337,334 but UKB and 23+Me removed = 143,210
# Mean value MDD never smokers = 1.20349, case prevalence = 20.349% of 194881: caseN = 39656, ctrlN = 155225
# Mean value MDD ever smokers = 1.25484, case prevalence = 25.484% of 160248: caseN = 40838, ctrlN = 119410

##Never smokers
outcome_mdd_dat_cpd_never_clumped$samplesize.exposure<-143210
outcome_mdd_dat_cpd_never_clumped$samplesize.outcome<- 194881
outcome_mdd_dat_cpd_never_clumped$ncase.outcome<-39656
outcome_mdd_dat_cpd_never_clumped$ncontrol.outcome<- 155225

# Adding a unit column labelled "log odds" and add a prevalence column
# Average prevalence of current-smoking amongst cohorts (-23&Me, -UKBiobank) = 40%; Supplementary Table 7 GSCAN (2019) paper
# Prevalence of depression in UK Biobank GWAS sample = 20%
outcome_mdd_dat_cpd_never_clumped$prevalence.outcome<-0.20
outcome_mdd_dat_cpd_never_clumped$units.exposure <-"SD"
outcome_mdd_dat_cpd_never_clumped$units.outcome <-"log odds"

#Run Steiger filtering
steiger <- steiger_filtering(outcome_mdd_dat_cpd_never_clumped)

#How many SNPs explain more variance in the exposure than in the outcome?
table(steiger$steiger_dir) #0 are FALSE
true_cpd_n_clmp<-subset(steiger, steiger$steiger_dir==TRUE)
str(true_cpd_n_clmp)

##Ever smokers
outcome_mdd_dat_cpd_ever_clumped$samplesize.exposure<-143851
outcome_mdd_dat_cpd_ever_clumped$samplesize.outcome<- 160248
outcome_mdd_dat_cpd_ever_clumped$ncase.outcome<-40838
outcome_mdd_dat_cpd_ever_clumped$ncontrol.outcome<- 119410

# Adding a unit column labelled "log odds" and add a prevalence column
# Average prevalence of current-smoking amongst cohorts (-23&Me, -UKBiobank) = 40%; Supplementary Table 7 GSCAN (2019) paper
# Prevalence of depression in UK Biobank GWAS sample = 25%
outcome_mdd_dat_cpd_ever_clumped$prevalence.outcome<-0.25
outcome_mdd_dat_cpd_ever_clumped$units.exposure <-"SD"
outcome_mdd_dat_cpd_ever_clumped$units.outcome <-"log odds"

#Run Steiger filtering
steiger <- steiger_filtering(outcome_mdd_dat_cpd_ever_clumped)

#How many SNPs explain more variance in the exposure than in the outcome?
table(steiger$steiger_dir) #0 are FALSE
true_cpd_e_clmp<-subset(steiger, steiger$steiger_dir==TRUE)
str(true_cpd_e_clmp)

################################
######### Plot results #########
################################
#https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#plots 

setwd(" ")
#Renaming exposure and outcome for plots
outcome_mdd_dat_cpd_never_clumped$exposure <- gsub("exposure", "Smoking Heaviness", outcome_mdd_dat_cpd_never_clumped$exposure)
outcome_mdd_dat_cpd_never_clumped$outcome<- gsub("outcome", "MDD", outcome_mdd_dat_cpd_never_clumped$outcome)
outcome_mdd_dat_cpd_ever_clumped$exposure <- gsub("exposure", "Smoking Heaviness", outcome_mdd_dat_cpd_ever_clumped$exposure)
outcome_mdd_dat_cpd_ever_clumped$outcome<- gsub("outcome", "MDD", outcome_mdd_dat_cpd_ever_clumped$outcome)

# single SNP analyses
cpd_never_clmp_single <- mr_singlesnp(outcome_mdd_dat_cpd_never_clumped, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
cpd_never_clmp_single <- mr_forest_plot(cpd_never_clmp_single)
ggsave(cpd_never_clmp_single[[1]], file="single_forest_cpd_never_clumped.png", width=7, height=10)

# single SNP analyses
cpd_ever_clmp_single <- mr_singlesnp(outcome_mdd_dat_cpd_ever_clumped, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
cpd_ever_clmp_single <- mr_forest_plot(cpd_ever_clmp_single)
ggsave(cpd_ever_clmp_single[[1]], file="single_forest_cpd_ever_clumped.png", width=7, height=10)

# leave one out analyses
cpd_never_clmp_loo <- mr_leaveoneout(outcome_mdd_dat_cpd_never_clumped)
cpd_never_clmp_loo <- mr_leaveoneout_plot(cpd_never_clmp_loo)
ggsave(cpd_never_clmp_loo[[1]], file="loo_cpd_never_clumped.png", width=7, height=8)

# leave one out analyses
cpd_ever_clmp_loo <- mr_leaveoneout(outcome_mdd_dat_cpd_ever_clumped)
cpd_ever_clmp_loo <- mr_leaveoneout_plot(cpd_ever_clmp_loo)
ggsave(cpd_ever_clmp_loo[[1]], file="loo_cpd_ever_clumped.png", width=7, height=8)

# scatter plot
plot_result_cpd_never_clumped <- subset(result_cpd_never_clumped, method != "MR PRESSO")
scatter_cpd_never_clmp <- mr_scatter_plot(plot_result_cpd_never_clumped, outcome_mdd_dat_cpd_never_clumped)
ggsave(scatter_cpd_never_clmp[[1]], file="scatter_cpd_never_clumped.png", width=7, height=7)

# scatter plot
scatter_cpd_ever_clmp <- mr_scatter_plot(result_cpd_ever_clumped, outcome_mdd_dat_cpd_ever_clumped)
ggsave(scatter_cpd_ever_clmp[[1]], file="scatter_cpd_ever_clumped.png", width=7, height=7)

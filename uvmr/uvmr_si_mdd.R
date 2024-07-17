# Script was created by Chloe Burke January 2024.
# The script conducts univariable MR exploring the effects of smoking initiation on lifetime depression.
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
#if (!require("devtools")) { install.packages("devtools") } else {}
#devtools::install_github("rondolab/MR-PRESSO")
library(MRPRESSO)
library(MendelianRandomization)
library(simex)
library(writexl)
library(LDlinkR)

##################################
#####  Set SNP list for SI   ##### 
##################################
setwd(" ")
SI_SNPlist<- read.table("SI_SNPlist.txt", header=FALSE) 
SI_SNPlist<-rename(SI_SNPlist, c("SNP" = "V1"))

############################################################
#####Extract exposure data for MR of SI and depression #####
############################################################
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

# Only 377/378 SNPs which is not an error
# Appears that rs111842178 is not available in the summary statistics without UK Biobank and 23&Me, given high number of SNPs - will not search for proxy

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

#######################################
##### Extract outcome data for MR #####
#######################################

setwd(" ")

########
## SI ##
########

#Smoking initiation is run in GSCAN in ever and never smokers, t/f we use the GWAS from the non-stratified Biobank sample
#316 of 316 SNPs available (i.e., 0 missing from outcome set)
mdd_dat_si_clumped <- read_outcome_data(
  "dep_full_imputed.txt",
  snps = si_dat_mr_clumped$SNP,
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
write.csv(si_dat_mr_clumped,"si_dat_uvmr.csv", row.names = FALSE)
write.csv(mdd_dat_si_clumped,"si_mdd_uvmr.csv", row.names = FALSE)

###################
## Reading files ##
###################
#setwd(" ")
#si_dat_mr_clumped <- read.csv("si_dat_uvmr.csv", header = T)
#mdd_dat_si_clumped <- read.csv("si_mdd_uvmr.csv", header = T)

############################################
##### Convert odds ratios to log odds #####
# Note: exposure already on the log scale.
############################################
# To convert binary outcomes from BOLT-LMM files to log odds, the following calculation needs to be made: logOR<-beta(caseprevelance*(1-caseprevelance))
# Case prevalence is found in the log file by searching Phenotype and looking at the mean value.
# As the phenotype was coded as 1 or 2, the case prevalence (%) is the mean -1*100
# Standard errors of SNP effect size estimates should also be divided by (case prevalence * (1 - case prevalence)) when applying the above transformation
# Mean value MDD all sample = 1.22652, case prevalence = 22.652% of 356641: caseN = 80786, ctrlN = 275855 

mdd_dat_si_clumped$beta.outcome <- as.numeric(as.character(
  mdd_dat_si_clumped$beta.outcome
))

mdd_dat_si_clumped["beta.outcome"] <- (mdd_dat_si_clumped$beta.outcome/(0.22652*((1-0.22652))))
mdd_dat_si_clumped["se.outcome"] <- (mdd_dat_si_clumped$se.outcome/(0.22652*((1-0.22652))))

#######################
##### Harmonising #####
#######################
##SI
# Palindromic with intermediate allele frequencies:
# rs10969352, rs1160685, rs13237637, rs17554906, rs1931431, rs1937443, rs3850736, 
# rs4140932, rs4326350, rs6932350, rs7598402, rs7920501, rs7921378, rs986714
# TRUE = 302/FALSE = 14
outcome_mdd_dat_si_clumped <-harmonise_data(si_dat_mr_clumped, mdd_dat_si_clumped, action = 2)
summary(outcome_mdd_dat_si_clumped$mr_keep == TRUE)

##################
## Saving files ##
##################
setwd(" ")
write.csv(outcome_mdd_dat_si_clumped,"outcome_mdd_si_uvmr.csv", row.names = FALSE)

###################
## Reading files ##
###################
#setwd(" ")
#outcome_mdd_dat_si_clumped <- read.csv("outcome_mdd_si_uvmr.csv", header = T)

################################################################################
################################## MR ##########################################
################################################################################

#IVW, Egger, Weighted Median + Weighted Mode
result_si_clumped <- mr(
  outcome_mdd_dat_si_clumped,
  method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
)
result_si_clumped <- generate_odds_ratios(result_si_clumped)
outcome_mdd_dat_si_clumped <- subset(outcome_mdd_dat_si_clumped, mr_keep)
mr_si_clmp <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_si_clumped$beta.exposure,
  outcome_mdd_dat_si_clumped$beta.outcome,
  outcome_mdd_dat_si_clumped$se.exposure,
  outcome_mdd_dat_si_clumped$se.outcome,
  parameters = default_parameters()
)
ptr_si_clmp <- data.frame(mr_si_clmp["Q"])
egger_si_clmp <- mr_egger_regression(
  outcome_mdd_dat_si_clumped$beta.exposure,
  outcome_mdd_dat_si_clumped$beta.outcome,
  outcome_mdd_dat_si_clumped$se.exposure,
  outcome_mdd_dat_si_clumped$se.outcome,
  parameters
)
ptr_si_clmp[2, 1] <- egger_si_clmp["Q"]
F <- abs(outcome_mdd_dat_si_clumped$beta.exposure)^2 / outcome_mdd_dat_si_clumped$se.exposure^2
mF <- mean(F)
ptr_si_clmp[3, 1] <- mF

# Note = at NbDistribution of 1000 warning message of:
# "Outlier test unstable. The significance threshold of 0.05 for the outlier test is not achievable with only 1000 to compute the null distribution." "The current precision is <0.302. Increase NbDistribution."
# Note = at NbDistribution of 2000 warning message of:
# "Outlier test unstable. The significance threshold of 0.05 for the outlier test is not achievable with only 2000 to compute the null distribution." "The current precision is <0.151. Increase NbDistribution."
# Upped until Nb = 7000 to meet p<0.05 threshold; upping the NDistribution increases the run time substantially (i.e., 45-60 mins); 
# Result interpretation advice here: https://github.com/rondolab/MR-PRESSO/issues/9 

si_mdd_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
          SdOutcome = "se.outcome", SdExposure = "se.exposure", 
          OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = outcome_mdd_dat_si_clumped, NbDistribution = 7000,  SignifThreshold = 0.05)

si_mdd_presso_mainres <- si_mdd_presso$`Main MR results`
si_mdd_presso_mainres <- si_mdd_presso_mainres %>%
  rename(b = `Causal Estimate`)
si_mdd_presso_mainres <- si_mdd_presso_mainres %>%
  rename(se = `Sd`)
si_mdd_presso_mainres <- generate_odds_ratios(si_mdd_presso_mainres)
si_mdd_presso_mainres$Distortion_P <- si_mdd_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
si_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT <- ifelse(si_mdd_presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, TRUE, FALSE)
summary(si_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT) #nOUT TRUE = 11 SNPs  
si_mdd_presso_mainres$nsnp <- c(".","291")
si_mdd_presso_mainres$nOUT_true <- c(".","11")
ptr_si_clmp[4, 1] <- si_mdd_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
  
#Save results
setwd(" ")
write.csv(si_mdd_presso_mainres, "presso_si_mdd_uvmr.csv", row.names=FALSE, quote=FALSE)

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
#Calculate the regression dilution Isq. 
BetaXG   = outcome_mdd_dat_si_clumped$beta.exposure
seBetaXG = outcome_mdd_dat_si_clumped$se.exposure 
seBetaYG <- outcome_mdd_dat_si_clumped$se.outcome

BXG  = abs(BetaXG)         # gene--exposure estimates are positive  

Isq_unweighted_si <- ISQ(BXG,seBetaXG) #unweighted = 0.16
Isq_weighted_si <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.07

#SIMEX correction not appropriate; don't report MR-Egger

#################################################
#########   Binding results tables   ############
#################################################

#Add other results into main results data frame (i.e.,  PRESSO)
#PRESSO:
#Outliers identified so this is the "outlier-corrected" estimate 
result_si_clumped[5, 1] <- result_si_clumped[1, "id.exposure"]
result_si_clumped[5,2] <- result_si_clumped[1, "id.outcome"]
result_si_clumped[5,3] <- result_si_clumped[1, "outcome"]
result_si_clumped[5,4] <- result_si_clumped[1, "exposure"]
result_si_clumped[5,5] <- "MR PRESSO"
result_si_clumped[5,6] <- si_mdd_presso_mainres[2, "nsnp"] 
result_si_clumped[5,7] <- si_mdd_presso_mainres[2, "b"]
result_si_clumped[5,8] <- si_mdd_presso_mainres[2, "se"]
result_si_clumped[5,9] <- si_mdd_presso_mainres[2, "P-value"]
result_si_clumped[5,10] <- si_mdd_presso_mainres[2, "lo_ci"]
result_si_clumped[5,11] <- si_mdd_presso_mainres[2, "up_ci"]
result_si_clumped[5,12] <- si_mdd_presso_mainres[2, "or"]
result_si_clumped[5,13] <- si_mdd_presso_mainres[2, "or_lci95"]
result_si_clumped[5,14] <- si_mdd_presso_mainres[2, "or_uci95"]
#Remove MR-Egger
result_si_clumped <- subset(result_si_clumped, method != "MR Egger")

#Save results
setwd(" ")
write.csv(result_si_clumped, "results_si_mdd_uvmr.csv", row.names=FALSE, quote=FALSE)

#####################################################
#########   Test for reverse causation   ############
#####################################################

### Steiger filtering ###
# Binary exposures/outcomes must contain: samplesize, ncase, ncontrol, prevalence and units="log odds"
# For smoking initiation there were 1,232,091 but now 599,289+ 73,331+ 39,480+ 270,820 were removed (23andMe and UKBB)
# Mean value MDD all sample = 1.22652, case prevalence = 22.652% of 356641: caseN = 80786, ctrlN = 275855 
outcome_mdd_dat_si_clumped$samplesize.exposure<-249171
outcome_mdd_dat_si_clumped$ncase.exposure<-140856
outcome_mdd_dat_si_clumped$ncontrol.exposure<-108314
outcome_mdd_dat_si_clumped$samplesize.outcome<- 356641
outcome_mdd_dat_si_clumped$ncase.outcome<-80786
outcome_mdd_dat_si_clumped$ncontrol.outcome<- 275855

# Adding a unit column labelled "log odds" and add a prevalence column
# Prevalence of depression in UK Biobank GWAS sample = 23%
outcome_mdd_dat_si_clumped$prevalence.exposure<-0.56
outcome_mdd_dat_si_clumped$prevalence.outcome<-0.23
outcome_mdd_dat_si_clumped$units.exposure <-"log odds"
outcome_mdd_dat_si_clumped$units.outcome <-"log odds"

#Run Steiger filtering
steiger <- steiger_filtering(outcome_mdd_dat_si_clumped)

#How many SNPs explain more variance in the exposure than in the outcome?
#Generates steiger_dir (which is `TRUE` if the rsq.exposure is larger than rsq.outcome)
table(steiger$steiger_dir) #40 are FALSE
true_si_clumped <-subset(steiger, steiger$steiger_dir==TRUE)
str(true_si_clumped)

#Re-running analysis restricted to SNPs that explain more variance in exposure than the outcome
result_si_clumped_true <- mr(
  true_si_clumped,
  method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
)
result_si_clumped_true <- generate_odds_ratios(result_si_clumped_true)
true_si_clumped <- subset(true_si_clumped, mr_keep)
mr_si_clumped_true <- TwoSampleMR::mr_ivw(
  true_si_clumped$beta.exposure,
  true_si_clumped$beta.outcome,
  true_si_clumped$se.exposure,
  true_si_clumped$se.outcome,
  parameters = default_parameters()
)
ptr_si_clmp_true <- data.frame(mr_si_clumped_true["Q"])
egger_si_clmp_true <- mr_egger_regression(
  true_si_clumped$beta.exposure,
  true_si_clumped$beta.outcome,
  true_si_clumped$se.exposure,
  true_si_clumped$se.outcome,
  parameters
)
ptr_si_clmp_true[2, 1] <- egger_si_clmp_true["Q"]
F <- abs(true_si_clumped$beta.exposure)^2 / true_si_clumped$se.exposure^2
mF <- mean(F)
ptr_si_clmp_true[3, 1] <- mF

si_mdd_true_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                           SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                           OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = true_si_clumped, NbDistribution = 7000,  SignifThreshold = 0.05)

si_mdd_true_presso_mainres <- si_mdd_true_presso$`Main MR results`
si_mdd_true_presso_mainres <- si_mdd_true_presso_mainres %>%
  rename(b = `Causal Estimate`)
si_mdd_true_presso_mainres <- si_mdd_true_presso_mainres %>%
  rename(se = `Sd`)
si_mdd_true_presso_mainres <- generate_odds_ratios(si_mdd_true_presso_mainres)
si_mdd_true_presso_mainres$Distortion_P <- si_mdd_true_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
si_mdd_true_presso$`MR-PRESSO results`$`Outlier Test`$nOUT <- ifelse(si_mdd_true_presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, TRUE, FALSE)
summary(si_mdd_true_presso$`MR-PRESSO results`$`Outlier Test`$nOUT) #nOUT TRUE = 3 SNPs  
si_mdd_true_presso_mainres$nsnp <- c(".","259")
si_mdd_true_presso_mainres$nOUT_true <- c(".","3")
ptr_si_clmp_true[4, 1] <- si_mdd_true_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue

#Calculate the regression dilution Isq. 
BetaXG   = true_si_clumped$beta.exposure
seBetaXG = true_si_clumped$se.exposure 
seBetaYG <- true_si_clumped$se.outcome

BXG  = abs(BetaXG)         # gene--exposure estimates are positive  

Isq_unweighted_si_true <- ISQ(BXG,seBetaXG) #unweighted = 0.10
Isq_weighted_si_true <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.004

#################################################
#########   Binding results tables   ############
#################################################

#Add other results into main results data frame (i.e.,  PRESSO)
#PRESSO:
#Outliers identified so this is the "outlier-corrected" estimate 
result_si_clumped_true[5, 1] <- result_si_clumped_true[1, "id.exposure"]
result_si_clumped_true[5,2] <- result_si_clumped_true[1, "id.outcome"]
result_si_clumped_true[5,3] <- result_si_clumped_true[1, "outcome"]
result_si_clumped_true[5,4] <- result_si_clumped_true[1, "exposure"]
result_si_clumped_true[5,5] <- "MR PRESSO(OC)"
result_si_clumped_true[5,6] <- si_mdd_true_presso_mainres[2, "nsnp"] 
result_si_clumped_true[5,7] <- si_mdd_true_presso_mainres[2, "b"]
result_si_clumped_true[5,8] <- si_mdd_true_presso_mainres[2, "se"]
result_si_clumped_true[5,9] <- si_mdd_true_presso_mainres[2, "P-value"]
result_si_clumped_true[5,10] <- si_mdd_true_presso_mainres[2, "lo_ci"]
result_si_clumped_true[5,11] <- si_mdd_true_presso_mainres[2, "up_ci"]
result_si_clumped_true[5,12] <- si_mdd_true_presso_mainres[2, "or"]
result_si_clumped_true[5,13] <- si_mdd_true_presso_mainres[2, "or_lci95"]
result_si_clumped_true[5,14] <- si_mdd_true_presso_mainres[2, "or_uci95"]
#Remove MR-Egger
result_si_clumped_true <- subset(result_si_clumped_true, method != "MR Egger")

#Save results
setwd(" ")
write.csv(result_si_clumped_true, "results_si_mdd_steiger.csv", row.names=FALSE, quote=FALSE)

################################
######### Plot results #########
################################
#https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#plots 

setwd(" ")
#Renaming exposure and outcome for plots
outcome_mdd_dat_si_clumped$exposure <- gsub("exposure", "Smoking Initiation", outcome_mdd_dat_si_clumped$exposure)
outcome_mdd_dat_si_clumped$outcome<- gsub("outcome", "MDD", outcome_mdd_dat_si_clumped$outcome)

# single SNP analyses
si_clmp_single <- mr_singlesnp(outcome_mdd_dat_si_clumped, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
si_clmp_single <- mr_forest_plot(si_clmp_single)
ggsave(si_clmp_single[[1]], file="single_forest_si_clumped.png", width=7, height=20)

# leave one out analyses
si_clmp_loo <- mr_leaveoneout(outcome_mdd_dat_si_clumped)
si_clmp_loo <- mr_leaveoneout_plot(si_clmp_loo)

ggsave(si_clmp_loo[[1]], file="loo_si_clumped.png", width=7, height=10)

# scatter plot
scatter_si_clmp <- mr_scatter_plot(result_si_clumped, outcome_mdd_dat_si_clumped)
ggsave(scatter_si_clmp[[1]], file="scatter_si_clumped.png", width=7, height=7)
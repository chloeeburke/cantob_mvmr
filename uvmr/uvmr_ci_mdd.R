# Script was created by Chloe Burke January 2024.
# The script conducts univariable MR exploring the effects of lifetime cannabis use on lifetime depression.
# The SNPs used in this script were selected based on the findings from a GWAS of lifetime cannabis use (Pasman et al, 2018)

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

################################
##### Set SNP list for CI  ##### 
################################
setwd(" ")
CI_SNPlist<- read.table("CI_SNPlist.txt", header=FALSE)
CI_SNPlist<-rename(CI_SNPlist, c("SNP" = "V1"))

##############################################################
##### Extract exposure data for MR of CUD and depression #####
##############################################################
setwd(" ")
ci_dat <- read_exposure_data("MA_ICC_23andme.txt", 
                              clump = FALSE, 
                              sep = " ", 
                              snp_col = "SNPID", 
                              beta_col = "Effect", 
                              se_col = "StdErr", 
                              eaf_col = "Freq1",
                              effect_allele_col = "Allele1", 
                              other_allele_col = "Allele2", 
                              pval_col = "PVAL", 
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

################################################################################
##### Extract outcome data for MR #####
################################################################################

setwd(" ")

########
## CI ##
########

#All sample
#No proxies needed, all 6 SNPs available in the outcome data
mdd_dat_ci_clumped <- read_outcome_data(
  "dep_full_imputed.txt",
  snps = ci_dat_mr_clumped$SNP,
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
write.csv(ci_dat_mr_clumped,"ci_dat_uvmr.csv", row.names = FALSE)
write.csv(mdd_dat_ci_clumped,"ci_mdd_uvmr.csv", row.names = FALSE)

###################################################
##### Convert odds ratios to log odds #####
# Note: the exposures are already on the log scale.
###################################################
# To convert binary outcomes from BOLT-LMM files to log odds, the following calculation needs to be made: logOR<-beta(caseprevelance*(1-caseprevelance))
# Case prevalence is found in the log file by searching Phenotype and looking at the mean value.
# As the phenotype was coded as 1 or 2, the case prevalence (%) is the mean -1*100
# Standard errors of SNP effect size estimates should also be divided by (case prevalence * (1 - case prevalence)) when applying the above transformation
# Mean value MDD all sample = 1.22652, case prevalence = 22.652% of 356641: caseN = 80786, ctrlN = 275855 

mdd_dat_ci_clumped$beta.outcome <- as.numeric(as.character(
  mdd_dat_ci_clumped$beta.outcome
))
mdd_dat_ci_clumped["beta.outcome"] <- (mdd_dat_ci_clumped$beta.outcome/(0.22652*((1-0.22652))))
mdd_dat_ci_clumped["se.outcome"] <- (mdd_dat_ci_clumped$se.outcome/(0.22652*((1-0.22652))))

#######################
##### Harmonising #####
#######################
outcome_mdd_dat_ci_clumped <-harmonise_data(ci_dat_mr_clumped, mdd_dat_ci_clumped, action = 2)

##################
## Saving files ##
##################
setwd(" ")
write.csv(outcome_mdd_dat_ci_clumped,"outcome_mdd_ci_uvmr.csv", row.names = FALSE)

###################
## Reading files ##
###################
#setwd(" ")
#outcome_mdd_dat_ci_clumped <- read.csv("outcome_mdd_ci_uvmr.csv", header = T)

################################################################################
################################## MR ##########################################
################################################################################

#IVW, Egger, Weighted Median + Weighted Mode
result_ci_clumped <- mr(
  outcome_mdd_dat_ci_clumped,
  method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
)
result_ci_clumped <- generate_odds_ratios(result_ci_clumped)
outcome_mdd_dat_ci_clumped <- subset(outcome_mdd_dat_ci_clumped, mr_keep)
mr_ci_clumped <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_ci_clumped$beta.exposure,
  outcome_mdd_dat_ci_clumped$beta.outcome,
  outcome_mdd_dat_ci_clumped$se.exposure,
  outcome_mdd_dat_ci_clumped$se.outcome,
  parameters = default_parameters()
)
ptr_ci_clumped <- data.frame(mr_ci_clumped["Q"])
egger_ci_clumped <- mr_egger_regression(
  outcome_mdd_dat_ci_clumped$beta.exposure,
  outcome_mdd_dat_ci_clumped$beta.outcome,
  outcome_mdd_dat_ci_clumped$se.exposure,
  outcome_mdd_dat_ci_clumped$se.outcome,
  parameters
)
ptr_ci_clumped[2, 1] <- egger_ci_clumped ["Q"]
F <- outcome_mdd_dat_ci_clumped$beta.exposure^2 / outcome_mdd_dat_ci_clumped$se.exposure^2
mF <- mean(F)
ptr_ci_clumped [3, 1] <- mF

#MR-PRESSO
ci_mdd_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                            SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                            OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = outcome_mdd_dat_ci_clumped, NbDistribution = 1000,  SignifThreshold = 0.05)

ci_mdd_presso_mainres <- ci_mdd_presso$`Main MR results`
ci_mdd_presso_mainres <- ci_mdd_presso_mainres %>%
  rename(b = `Causal Estimate`)
ci_mdd_presso_mainres <- ci_mdd_presso_mainres %>%
  rename(se = `Sd`)
ci_mdd_presso_mainres <- generate_odds_ratios(ci_mdd_presso_mainres)
ci_mdd_presso_mainres$Distortion_P <- ci_mdd_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
ci_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT <- ifelse(ci_mdd_presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, TRUE, FALSE)
summary(ci_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT) #No outliers identified
ptr_ci_clumped [4, 1] <- NA

#Save results
setwd(" ")
write.csv(ci_mdd_presso_mainres, "presso_ci_mdd_uvmr.csv", row.names=FALSE, quote=FALSE)

################################################################################
##### Calculate I2GX #####
# Isq(y, s) where y = vector of effects and s = vector of standard errors
# Apply simulation extrapolation SIMEX corrections to MR-Egger analysis where
# I2GX estimates < 0.9
# This indicates the effect estimate is biased by 10% due to measurement error
# (Bowden, Del Greco, et al., 2016)
# (Lederer & K?chenhoff, 2006)
################################################################################
##Two-sample MR package:
#Doesn't have gene-exposure estimates as positive
#Isq(outcome_mdd_dat_ci_clumped$beta.exposure, outcome_mdd_dat_ci_clumped$se.exposure) #0.95

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
BetaXG   = outcome_mdd_dat_ci_clumped$beta.exposure
seBetaXG = outcome_mdd_dat_ci_clumped$se.exposure 
seBetaYG <- outcome_mdd_dat_ci_clumped$se.outcome

BXG  = abs(BetaXG)         # gene--exposure estimates are positive  

Isq_unweighted_ci <- ISQ(BXG,seBetaXG) #unweighted = 0.69
Isq_weighted_ci <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.05

## SIMEX Correction
# If the Isquared (weighted or unweighted) is between 0.6 and 0.9, conduct a SIMEX correction
BetaXG   = outcome_mdd_dat_ci_clumped$beta.exposure
BetaYG <- outcome_mdd_dat_ci_clumped$beta.outcome
seBetaXG = outcome_mdd_dat_ci_clumped$se.exposure
seBetaYG <- outcome_mdd_dat_ci_clumped$se.outcome

BYG <- BetaYG*sign(BetaXG) # Pre-processing steps to ensure all gene--exposure estimates are positive
BXG <- abs(BetaXG)  # gene--exposure estimates are positive  

# MR-Egger regression (weighted) 
#Fit1 <- lm(BYG ~ BXG,weights=1/seBetaYG^2,x=TRUE,y=TRUE)

# MR-Egger regression (unweighted)
Fit2 <- lm(BYG~BXG,x=TRUE,y=TRUE) 

# Simulation extrapolation 
#mod.sim1 <- simex(Fit1,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
mod.sim2 <- simex(Fit2,B=10000, measurement.error = seBetaXG, SIMEXvariable="BXG",fitting.method ="quad",asymptotic="FALSE") 
#mod1<-summary(mod.sim1)
mod2<-summary(mod.sim2)
simex_ci <- data.frame(mod2$coefficients)
simex_ci <- simex_ci %>%
  rename(b = `jackknife.Estimate`, se = 'jackknife.Std..Error', pval = 'jackknife.Pr...t..')
simex_ci <- generate_odds_ratios(simex_ci)

#Save results
setwd(" ")
write.csv(simex_ci, "simex_ci_mdd_uvmr.csv", row.names=FALSE, quote=FALSE)

#################################################
#########   Binding results tables   ############
#################################################

#Add other results into main results data frame (i.e., SIMEX + PRESSO)
#SIMEX:
result_ci_clumped[5, 1] <- result_ci_clumped[1, "id.exposure"]
result_ci_clumped[5,2] <- result_ci_clumped[1, "id.outcome"]
result_ci_clumped[5,3] <- result_ci_clumped[1, "outcome"]
result_ci_clumped[5,4] <- result_ci_clumped[1, "exposure"]
result_ci_clumped[5,5] <- "MR Egger (SIMEX)"
result_ci_clumped[5,6] <- result_ci_clumped[1, "nsnp"]
result_ci_clumped[5,7] <- simex_ci[2, "b"]
result_ci_clumped[5,8] <- simex_ci[2, "se"]
result_ci_clumped[5,9] <- simex_ci[2, "pval"]
result_ci_clumped[5,10] <- simex_ci[2, "lo_ci"]
result_ci_clumped[5,11] <- simex_ci[2, "up_ci"]
result_ci_clumped[5,12] <- simex_ci[2, "or"]
result_ci_clumped[5,13] <- simex_ci[2, "or_lci95"]
result_ci_clumped[5,14] <- simex_ci[2, "or_uci95"]
#PRESSO:
#No outliers identified so this is the "raw" estimate (i.e., not outlier corrected)
result_ci_clumped[6, 1] <- result_ci_clumped[1, "id.exposure"]
result_ci_clumped[6,2] <- result_ci_clumped[1, "id.outcome"]
result_ci_clumped[6,3] <- result_ci_clumped[1, "outcome"]
result_ci_clumped[6,4] <- result_ci_clumped[1, "exposure"]
result_ci_clumped[6,5] <- "MR PRESSO"
result_ci_clumped[6,6] <- result_ci_clumped[1, "nsnp"] #No outliers removed so same as other analyses
result_ci_clumped[6,7] <- ci_mdd_presso_mainres[1, "b"]
result_ci_clumped[6,8] <- ci_mdd_presso_mainres[1, "se"]
result_ci_clumped[6,9] <- ci_mdd_presso_mainres[1, "P-value"]
result_ci_clumped[6,10] <- ci_mdd_presso_mainres[1, "lo_ci"]
result_ci_clumped[6,11] <- ci_mdd_presso_mainres[1, "up_ci"]
result_ci_clumped[6,12] <- ci_mdd_presso_mainres[1, "or"]
result_ci_clumped[6,13] <- ci_mdd_presso_mainres[1, "or_lci95"]
result_ci_clumped[6,14] <- ci_mdd_presso_mainres[1, "or_uci95"]
#Remove MR-Egger
result_ci_clumped <- subset(result_ci_clumped, method != "MR Egger")

#Save results
setwd(" ")
write.csv(result_ci_clumped, "results_ci_mdd_uvmr.csv", row.names=FALSE, quote=FALSE)

#####################################################
#########   Test for reverse causation   ############
#####################################################
### Steiger filtering ###
# Binary exposures/outcomes must contain: samplesize, ncase, ncontrol, prevalence and units="log odds"
# ICC N = 35297 (42.8% cases i.e., 15107), 23+Me N = 22683 (43.2% cases, i.e., 9799), total = 57980, cases = 24906, ctrls = 33074
# Mean value MDD all sample = 1.22652, case prevalence = 22.652% of 356641: caseN = 80786, ctrlN = 275855 
outcome_mdd_dat_ci_clumped$samplesize.exposure<-57980
outcome_mdd_dat_ci_clumped$ncase.exposure<-24906
outcome_mdd_dat_ci_clumped$ncontrol.exposure<- 33074
outcome_mdd_dat_ci_clumped$samplesize.outcome<- 356641
outcome_mdd_dat_ci_clumped$ncase.outcome<-80786
outcome_mdd_dat_ci_clumped$ncontrol.outcome<- 275855

# Adding a unit column labelled "log odds" and add a prevalence column
outcome_mdd_dat_ci_clumped$prevalence.exposure<-0.43
outcome_mdd_dat_ci_clumped$prevalence.outcome<-0.23
outcome_mdd_dat_ci_clumped$units.exposure <- "log odds"
outcome_mdd_dat_ci_clumped$units.outcome <- "log odds"

#Run Steiger filtering
steiger <- steiger_filtering(outcome_mdd_dat_ci_clumped)

#How many SNPs explain more variance in the exposure than in the outcome?
#Generates steiger_dir (which is `TRUE` if the rsq.exposure is larger than rsq.outcome)
table(steiger$steiger_dir) #0 are FALSE
true_ci_clumped <-subset(steiger, steiger$steiger_dir==TRUE)
str(true_ci_clumped)

################################
######### Plot results #########
################################
#https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#plots 
setwd(" ")

#Renaming exposure and outcome for plots
outcome_mdd_dat_ci_clumped$exposure <- gsub("exposure", "CI", outcome_mdd_dat_ci_clumped$exposure)
outcome_mdd_dat_ci_clumped$outcome<- gsub("outcome", "MDD", outcome_mdd_dat_ci_clumped$outcome)

# single SNP analyses
ci_clmp_single <- mr_singlesnp(outcome_mdd_dat_ci_clumped, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
ci_clmp_single <- mr_forest_plot(ci_clmp_single)
ggsave(ci_clmp_single[[1]], file="single_forest_ci_clumped.png", width=7, height=4)

# leave one out analyses
ci_clmp_loo <- mr_leaveoneout(outcome_mdd_dat_ci_clumped)
ci_clmp_loo <- mr_leaveoneout_plot(ci_clmp_loo)
ggsave(ci_clmp_loo[[1]], file="loo_ci_clumped.png", width=7, height=4)

# scatter plot
scatter_ci_clmp <- mr_scatter_plot(result_ci_clumped, outcome_mdd_dat_ci_clumped)
ggsave(scatter_ci_clmp[[1]], file="scatter_ci_clumped.png", width=7, height=7)
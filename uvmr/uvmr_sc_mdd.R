# Script was created by Chloe Burke January 2024.
# The script conducts univariable MR exploring the effects of smoking continuation (i.e. current vs former smokers) on lifetime depression.
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
library(MRPRESSO)

##################################
#####  Set SNP list for SC   ##### 
##################################
setwd(" ")
SC_SNPlist<- read.table("SC_SNPlist.txt", header=FALSE) 
SC_SNPlist<-rename(SC_SNPlist, c("SNP" = "V1"))

################################################################################
##### Extract exposure data for MR of SC and depression #####
################################################################################
setwd(" ")
sc_dat <- read_exposure_data("SmokingCessation.WithoutUKB.txt",
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

sc_dat_mr <- format_data(
  sc_dat,
  type = "exposure",
  snps = SC_SNPlist$SNP,
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
#####   Clumping SC instrument    #####
#######################################
##Restrict the exposure data to independent SNPs only
sc_dat_mr_clumped <- clump_data(sc_dat_mr, clump_kb = 500, clump_r2 = 0.001) #17 SNPs

################################################################################
##### Extract outcome data for MR #####
################################################################################
setwd(" ")

#########
## SC ##
#########
#Smoking cessation is run in GSCAN in ever smokers, t/f we use the GWAS from the stratified UK Biobank sample
#17 of 17 SNPs available (i.e., 0 missing from outcome set)
#Never smokers
mdd_dat_never_sc_clumped <- read_outcome_data(
  "dep_never_imputed.txt",
  snps = sc_dat_mr_clumped$SNP,
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
#17 of 17 SNPs available (i.e., 0 missing from outcome set)
mdd_dat_ever_sc_clumped <- read_outcome_data(
  "dep_ever_imputed.txt",
  snps = sc_dat_mr_clumped$SNP,
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
write.csv(sc_dat_mr_clumped,"sc_dat_uvmr.csv", row.names = FALSE)
write.csv(mdd_dat_never_sc_clumped,"sc_mdd_uvmr_never.csv", row.names = FALSE)
write.csv(mdd_dat_ever_sc_clumped,"sc_mdd_uvmr_ever.csv", row.names = FALSE)

###################
## Reading files ##
###################
#setwd(" ")
#sc_dat_mr_clumped <- read.csv("sc_dat_uvmr.csv", header = T)
#mdd_dat_never_sc_clumped <- read.csv("sc_mdd_uvmr_never.csv", header = T)
#mdd_dat_ever_sc_clumped <- read.csv("sc_mdd_uvmr_ever.csv", header = T)

################################################################################
##### Convert odds ratios to log odds #####
# Note: the exposures are already on the log scale.
################################################################################
# To convert binary outcomes from BOLT-LMM files to log odds, the following calculation needs to be made: logOR<-beta(caseprevelance*(1-caseprevelance))
# Case prevalence is found in the log file by searching Phenotype and looking at the mean value.
# As the phenotype was coded as 1 or 2, the case prevalence (%) is the mean -1*100
# Standard errors of SNP effect size estimates should also be divided by (case prevalence * (1 - case prevalence)) when applying the above transformation
# Mean value MDD never smokers = 1.20349, case prevalence = 20.349% of 194881: caseN = 39656, ctrlN = 155225
# Mean value MDD ever smokers = 1.25484, case prevalence = 25.484% of 160248: caseN = 40838, ctrlN = 119410

mdd_dat_never_sc_clumped$beta.outcome <- as.numeric(as.character(
  mdd_dat_never_sc_clumped$beta.outcome
))
mdd_dat_ever_sc_clumped$beta.outcome <- as.numeric(as.character(
  mdd_dat_ever_sc_clumped$beta.outcome
))

mdd_dat_never_sc_clumped["beta.outcome"] <- (mdd_dat_never_sc_clumped$beta.outcome/(0.20349*((1-0.20349))))
mdd_dat_never_sc_clumped["se.outcome"] <- (mdd_dat_never_sc_clumped$se.outcome/(0.20349*((0.20349))))
mdd_dat_ever_sc_clumped["beta.outcome"] <- (mdd_dat_ever_sc_clumped$beta.outcome/(0.25484*((1-0.25484))))
mdd_dat_ever_sc_clumped["se.outcome"] <- (mdd_dat_ever_sc_clumped$se.outcome/(0.25484*((1-0.25484))))

################################################################################
##### Harmonising #####
################################################################################
##SC - none removed
outcome_mdd_dat_sc_never_clumped <-harmonise_data(sc_dat_mr_clumped, mdd_dat_never_sc_clumped, action = 2)
outcome_mdd_dat_sc_ever_clumped <-harmonise_data(sc_dat_mr_clumped, mdd_dat_ever_sc_clumped, action = 2)

##################
## Saving files ##
##################
setwd(" ")
write.csv(outcome_mdd_dat_sc_never_clumped,"outcome_mdd_sc_uvmr_never.csv", row.names = FALSE)
write.csv(outcome_mdd_dat_sc_ever_clumped,"outcome_mdd_sc_uvmr_ever.csv", row.names = FALSE)

###################
## Reading files ##
###################
#setwd("//rdsfcifs.acrc.bris.ac.uk/UK_Biobank_Project_9142/Working folder Chloe/Data/2sMR Cannabis Tobacco")
#outcome_mdd_dat_sc_never_clumped <- read.csv("outcome_mdd_sc_uvmr_never.csv", header = T)
#outcome_mdd_dat_sc_ever_clumped <- read.csv("outcome_mdd_sc_uvmr_ever.csv", header = T)

################################################################################
################################## MR ##########################################
################################################################################

# SC never smokers
# IVW, Egger, Weighted Median + Weighted Mode
result_sc_never_clumped <- mr(
  outcome_mdd_dat_sc_never_clumped,
  method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
)
result_sc_never_clumped <- generate_odds_ratios(result_sc_never_clumped)
outcome_mdd_dat_sc_never_clumped <- subset(outcome_mdd_dat_sc_never_clumped, mr_keep)
mr_n_sc_clmp <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_sc_never_clumped$beta.exposure,
  outcome_mdd_dat_sc_never_clumped$beta.outcome,
  outcome_mdd_dat_sc_never_clumped$se.exposure,
  outcome_mdd_dat_sc_never_clumped$se.outcome,
  parameters = default_parameters()
)
ptr_n_sc_clmp <- data.frame(mr_n_sc_clmp["Q"])
egger_n_sc_clmp <- mr_egger_regression(
  outcome_mdd_dat_sc_never_clumped$beta.exposure,
  outcome_mdd_dat_sc_never_clumped$beta.outcome,
  outcome_mdd_dat_sc_never_clumped$se.exposure,
  outcome_mdd_dat_sc_never_clumped$se.outcome,
  parameters
)
ptr_n_sc_clmp[2, 1] <- egger_n_sc_clmp["Q"]
F <- abs(outcome_mdd_dat_sc_never_clumped$beta.exposure)^2 / outcome_mdd_dat_sc_never_clumped$se.exposure^2
mF <- mean(F)
ptr_n_sc_clmp[3, 1] <- mF

#MR-PRESSO
sc_n_mdd_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                             SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                             OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = outcome_mdd_dat_sc_never_clumped, NbDistribution = 1000,  SignifThreshold = 0.05)

sc_n_mdd_presso_mainres <- sc_n_mdd_presso$`Main MR results`
sc_n_mdd_presso_mainres <- sc_n_mdd_presso_mainres %>%
  rename(b = `Causal Estimate`)
sc_n_mdd_presso_mainres <- sc_n_mdd_presso_mainres %>%
  rename(se = `Sd`)
sc_n_mdd_presso_mainres <- generate_odds_ratios(sc_n_mdd_presso_mainres)
sc_n_mdd_presso_mainres$Distortion_P <- sc_n_mdd_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
sc_n_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT <- ifelse(sc_n_mdd_presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, TRUE, FALSE)
summary(sc_n_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT) #0 outliers identified 
ptr_n_sc_clmp[4, 1] <- NA


# SC ever smokers
# IVW, Egger, Weighted Median + Weighted Mode
result_sc_ever_clumped <- mr(
  outcome_mdd_dat_sc_ever_clumped,
  method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
)
result_sc_ever_clumped <- generate_odds_ratios(result_sc_ever_clumped)
outcome_mdd_dat_sc_ever_clumped <- subset(outcome_mdd_dat_sc_ever_clumped, mr_keep)
mr_e_sc_clmp <- TwoSampleMR::mr_ivw(
  outcome_mdd_dat_sc_ever_clumped$beta.exposure,
  outcome_mdd_dat_sc_ever_clumped$beta.outcome,
  outcome_mdd_dat_sc_ever_clumped$se.exposure,
  outcome_mdd_dat_sc_ever_clumped$se.outcome,
  parameters = default_parameters()
)
ptr_e_sc_clmp <- data.frame(mr_e_sc_clmp["Q"])
egger_e_sc_clmp <- mr_egger_regression(
  outcome_mdd_dat_sc_ever_clumped$beta.exposure,
  outcome_mdd_dat_sc_ever_clumped$beta.outcome,
  outcome_mdd_dat_sc_ever_clumped$se.exposure,
  outcome_mdd_dat_sc_ever_clumped$se.outcome,
  parameters
)
ptr_e_sc_clmp[2, 1] <- egger_e_sc_clmp["Q"]
F <- abs(outcome_mdd_dat_sc_ever_clumped$beta.exposure)^2 / outcome_mdd_dat_sc_ever_clumped$se.exposure^2
mF <- mean(F)
ptr_e_sc_clmp[3, 1] <- mF

#MR-PRESSO
sc_e_mdd_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                           SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                           OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = outcome_mdd_dat_sc_ever_clumped, NbDistribution = 1000,  SignifThreshold = 0.05)

sc_e_mdd_presso_mainres <- sc_e_mdd_presso$`Main MR results`
sc_e_mdd_presso_mainres <- sc_e_mdd_presso_mainres %>%
  rename(b = `Causal Estimate`)
sc_e_mdd_presso_mainres <- sc_e_mdd_presso_mainres %>%
  rename(se = `Sd`)
sc_e_mdd_presso_mainres <- generate_odds_ratios(sc_e_mdd_presso_mainres)
sc_e_mdd_presso_mainres$Distortion_P <- sc_e_mdd_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
sc_e_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT <- ifelse(sc_e_mdd_presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, TRUE, FALSE)
summary(sc_e_mdd_presso$`MR-PRESSO results`$`Outlier Test`$nOUT) #nOUT TRUE = 1 SNPs  
sc_e_mdd_presso_mainres$nsnp <- c(".","16")
sc_e_mdd_presso_mainres$nOUT_true <- c(".","1")
ptr_e_sc_clmp[4, 1] <- sc_e_mdd_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue

setwd(" ")
write.csv(sc_n_mdd_presso_mainres, "presso_sc_never.csv", row.names=FALSE, quote=FALSE)
write.csv(sc_e_mdd_presso_mainres, "presso_sc_ever.csv", row.names=FALSE, quote=FALSE)

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
BetaXG   = outcome_mdd_dat_sc_never_clumped$beta.exposure
seBetaXG = outcome_mdd_dat_sc_never_clumped$se.exposure 
seBetaYG <- outcome_mdd_dat_sc_never_clumped$se.outcome

BXG  = abs(BetaXG)         # gene--exposure estimates are positive  

Isq_unweighted_sc_n <- ISQ(BXG,seBetaXG) #unweighted = 0.69
Isq_weighted_sc_n <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.58

##Ever smokers
#Calculate the regression dilution Isq. 
BetaXG   = outcome_mdd_dat_sc_ever_clumped$beta.exposure
seBetaXG = outcome_mdd_dat_sc_ever_clumped$se.exposure 
seBetaYG <- outcome_mdd_dat_sc_ever_clumped$se.outcome

BXG  = abs(BetaXG)         # gene--exposure estimates are positive  

Isq_unweighted_sc_e <- ISQ(BXG,seBetaXG) #unweighted = 0.69
Isq_weighted_sc_e <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.58

## SIMEX Correction - Never Smokers
# If the Isquared (weighted or unweighted) is between 0.6 and 0.9, conduct a SIMEX correction
BetaXG   = outcome_mdd_dat_sc_never_clumped$beta.exposure
BetaYG <- outcome_mdd_dat_sc_never_clumped$beta.outcome
seBetaXG = outcome_mdd_dat_sc_never_clumped$se.exposure
seBetaYG <- outcome_mdd_dat_sc_never_clumped$se.outcome

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
simex_sc_n <- data.frame(mod2$coefficients)
simex_sc_n <- simex_sc_n %>%
  rename(b = `jackknife.Estimate`, se = 'jackknife.Std..Error', pval = 'jackknife.Pr...t..')
simex_sc_n <- generate_odds_ratios(simex_sc_n)

## SIMEX Correction - Ever Smokers
# If the Isquared (weighted or unweighted) is between 0.6 and 0.9, conduct a SIMEX correction
BetaXG   = outcome_mdd_dat_sc_ever_clumped$beta.exposure
BetaYG <- outcome_mdd_dat_sc_ever_clumped$beta.outcome
seBetaXG = outcome_mdd_dat_sc_ever_clumped$se.exposure
seBetaYG <- outcome_mdd_dat_sc_ever_clumped$se.outcome

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
simex_sc_e <- data.frame(mod2$coefficients)
simex_sc_e <- simex_sc_e %>%
  rename(b = `jackknife.Estimate`, se = 'jackknife.Std..Error', pval = 'jackknife.Pr...t..')
simex_sc_e <- generate_odds_ratios(simex_sc_e)

#Save results
setwd(" ")
write.csv(simex_sc_n, "simex_sc_never.csv", row.names=FALSE, quote=FALSE)
write.csv(simex_sc_e, "simex_sc_ever.csv", row.names=FALSE, quote=FALSE)

#################################################
#########   Binding results tables   ############
#################################################

##Never smokers
#Add other results into main results data frame (i.e., SIMEX + PRESSO)
#SIMEX:
result_sc_never_clumped[5, 1] <- result_sc_never_clumped[1, "id.exposure"]
result_sc_never_clumped[5,2] <- result_sc_never_clumped[1, "id.outcome"]
result_sc_never_clumped[5,3] <- result_sc_never_clumped[1, "outcome"]
result_sc_never_clumped[5,4] <- result_sc_never_clumped[1, "exposure"]
result_sc_never_clumped[5,5] <- "MR Egger (SIMEX)"
result_sc_never_clumped[5,6] <- result_sc_never_clumped[1, "nsnp"]
result_sc_never_clumped[5,7] <- simex_sc_n[2, "b"]
result_sc_never_clumped[5,8] <- simex_sc_n[2, "se"]
result_sc_never_clumped[5,9] <- simex_sc_n[2, "pval"]
result_sc_never_clumped[5,10] <- simex_sc_n[2, "lo_ci"]
result_sc_never_clumped[5,11] <- simex_sc_n[2, "up_ci"]
result_sc_never_clumped[5,12] <- simex_sc_n[2, "or"]
result_sc_never_clumped[5,13] <- simex_sc_n[2, "or_lci95"]
result_sc_never_clumped[5,14] <- simex_sc_n[2, "or_uci95"]
#PRESSO:
#No outliers identified so this is the "raw" estimate (i.e., not outlier corrected)
result_sc_never_clumped[6, 1] <- result_sc_never_clumped[1, "id.exposure"]
result_sc_never_clumped[6,2] <- result_sc_never_clumped[1, "id.outcome"]
result_sc_never_clumped[6,3] <- result_sc_never_clumped[1, "outcome"]
result_sc_never_clumped[6,4] <- result_sc_never_clumped[1, "exposure"]
result_sc_never_clumped[6,5] <- "MR PRESSO"
result_sc_never_clumped[6,6] <- result_sc_never_clumped[1, "nsnp"] 
result_sc_never_clumped[6,7] <- sc_n_mdd_presso_mainres[1, "b"]
result_sc_never_clumped[6,8] <- sc_n_mdd_presso_mainres[1, "se"]
result_sc_never_clumped[6,9] <- sc_n_mdd_presso_mainres[1, "P-value"]
result_sc_never_clumped[6,10] <- sc_n_mdd_presso_mainres[1, "lo_ci"]
result_sc_never_clumped[6,11] <- sc_n_mdd_presso_mainres[1, "up_ci"]
result_sc_never_clumped[6,12] <- sc_n_mdd_presso_mainres[1, "or"]
result_sc_never_clumped[6,13] <- sc_n_mdd_presso_mainres[1, "or_lci95"]
result_sc_never_clumped[6,14] <- sc_n_mdd_presso_mainres[1, "or_uci95"]
#Remove MR-Egger
result_sc_never_clumped <- subset(result_sc_never_clumped, method != "MR Egger")

##Ever smokers
#Add other results into main results data frame (i.e., SIMEX + PRESSO)
#SIMEX:
result_sc_ever_clumped[5, 1] <- result_sc_ever_clumped[1, "id.exposure"]
result_sc_ever_clumped[5,2] <- result_sc_ever_clumped[1, "id.outcome"]
result_sc_ever_clumped[5,3] <- result_sc_ever_clumped[1, "outcome"]
result_sc_ever_clumped[5,4] <- result_sc_ever_clumped[1, "exposure"]
result_sc_ever_clumped[5,5] <- "MR Egger (SIMEX)"
result_sc_ever_clumped[5,6] <- result_sc_ever_clumped[1, "nsnp"]
result_sc_ever_clumped[5,7] <- simex_sc_e[2, "b"]
result_sc_ever_clumped[5,8] <- simex_sc_e[2, "se"]
result_sc_ever_clumped[5,9] <- simex_sc_e[2, "pval"]
result_sc_ever_clumped[5,10] <- simex_sc_e[2, "lo_ci"]
result_sc_ever_clumped[5,11] <- simex_sc_e[2, "up_ci"]
result_sc_ever_clumped[5,12] <- simex_sc_e[2, "or"]
result_sc_ever_clumped[5,13] <- simex_sc_e[2, "or_lci95"]
result_sc_ever_clumped[5,14] <- simex_sc_e[2, "or_uci95"]
#PRESSO:
#Outliers identified so this is the "OC" estimate (i.e., outlier corrected)
result_sc_ever_clumped[6, 1] <- result_sc_ever_clumped[1, "id.exposure"]
result_sc_ever_clumped[6,2] <- result_sc_ever_clumped[1, "id.outcome"]
result_sc_ever_clumped[6,3] <- result_sc_ever_clumped[1, "outcome"]
result_sc_ever_clumped[6,4] <- result_sc_ever_clumped[1, "exposure"]
result_sc_ever_clumped[6,5] <- "MR PRESSO"
result_sc_ever_clumped[6,6] <- sc_e_mdd_presso_mainres[2, "nsnp"] 
result_sc_ever_clumped[6,7] <- sc_e_mdd_presso_mainres[2, "b"]
result_sc_ever_clumped[6,8] <- sc_e_mdd_presso_mainres[2, "se"]
result_sc_ever_clumped[6,9] <- sc_e_mdd_presso_mainres[2, "P-value"]
result_sc_ever_clumped[6,10] <- sc_e_mdd_presso_mainres[2, "lo_ci"]
result_sc_ever_clumped[6,11] <- sc_e_mdd_presso_mainres[2, "up_ci"]
result_sc_ever_clumped[6,12] <- sc_e_mdd_presso_mainres[2, "or"]
result_sc_ever_clumped[6,13] <- sc_e_mdd_presso_mainres[2, "or_lci95"]
result_sc_ever_clumped[6,14] <- sc_e_mdd_presso_mainres[2, "or_uci95"]
#Remove MR-Egger
result_sc_ever_clumped <- subset(result_sc_ever_clumped, method != "MR Egger")

#Save tables
setwd(" ")
write.csv(result_sc_never_clumped, "results_sc_never.csv", row.names=FALSE, quote=FALSE)
write.csv(result_sc_ever_clumped, "results_sc_ever.csv", row.names=FALSE, quote=FALSE)


#####################################################
#########   Test for reverse causation   ############
#####################################################

### Steiger filtering ###
# Binary exposures/outcomes must contain: samplesize, ncase, ncontrol, prevalence and units="log odds"
# For smoking cessation there were were 547,219 but UKB and 23+Me removed = 143851
# Mean value MDD never smokers = 1.20349, case prevalence = 20.349% of 194881: caseN = 39656, ctrlN = 155225
# Mean value MDD ever smokers = 1.25484, case prevalence = 25.484% of 160248: caseN = 40838, ctrlN = 119410

##Never smokers
outcome_mdd_dat_sc_never_clumped$samplesize.exposure<-143851
outcome_mdd_dat_sc_never_clumped$ncase.exposure <- 65251
outcome_mdd_dat_sc_never_clumped$ncontrol.exposure <- 78599
outcome_mdd_dat_sc_never_clumped$samplesize.outcome<- 194881
outcome_mdd_dat_sc_never_clumped$ncase.outcome<-39656
outcome_mdd_dat_sc_never_clumped$ncontrol.outcome<- 155225

# Adding a unit column labelled "log odds" and add a prevalence column
# Average prevalence of current-smoking amongst cohorts (-23&Me, -UKBiobank) = 40%; Supplementary Table 7 GSCAN (2019) paper
# Prevalence of depression in UK Biobank GWAS sample = 20%
outcome_mdd_dat_sc_never_clumped$prevalence.exposure<-0.45
outcome_mdd_dat_sc_never_clumped$prevalence.outcome<-0.20
outcome_mdd_dat_sc_never_clumped$units.exposure <-"log odds"
outcome_mdd_dat_sc_never_clumped$units.outcome <-"log odds"

#Run Steiger filtering
steiger <- steiger_filtering(outcome_mdd_dat_sc_never_clumped)

#How many SNPs explain more variance in the exposure than in the outcome?
table(steiger$steiger_dir) #1 is FALSE
true_sc_n_clmp<-subset(steiger, steiger$steiger_dir==TRUE)
str(true_sc_n_clmp)

##Ever smokers
outcome_mdd_dat_sc_ever_clumped$samplesize.exposure<-143851
outcome_mdd_dat_sc_ever_clumped$ncase.exposure <- 65251
outcome_mdd_dat_sc_ever_clumped$ncontrol.exposure <- 78599
outcome_mdd_dat_sc_ever_clumped$samplesize.outcome<- 160248
outcome_mdd_dat_sc_ever_clumped$ncase.outcome<-40838
outcome_mdd_dat_sc_ever_clumped$ncontrol.outcome<- 119410

# Adding a unit column labelled "log odds" and add a prevalence column
# Prevalence of depression in UK Biobank GWAS sample = 25%
outcome_mdd_dat_sc_ever_clumped$prevalence.exposure<-0.45
outcome_mdd_dat_sc_ever_clumped$prevalence.outcome<-0.25
outcome_mdd_dat_sc_ever_clumped$units.exposure <-"log odds"
outcome_mdd_dat_sc_ever_clumped$units.outcome <-"log odds"

#Run Steiger filtering
steiger <- steiger_filtering(outcome_mdd_dat_sc_ever_clumped)

#How many SNPs explain more variance in the exposure than in the outcome?
table(steiger$steiger_dir) #1 is FALSE
true_sc_e_clmp<-subset(steiger, steiger$steiger_dir==TRUE)
str(true_sc_e_clmp)

#Re-running analysis for ever smokers restricted to SNPs that explain more variance in exposure than the outcome
result_sc_ever_clumped_true <- mr(
  true_sc_e_clmp,
  method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
)
result_sc_ever_clumped_true <- generate_odds_ratios(result_sc_ever_clumped_true)
true_sc_e_clmp <- subset(true_sc_e_clmp, mr_keep)
mr_e_sc_clmp_true <- TwoSampleMR::mr_ivw(
  true_sc_e_clmp$beta.exposure,
  true_sc_e_clmp$beta.outcome,
  true_sc_e_clmp$se.exposure,
  true_sc_e_clmp$se.outcome,
  parameters = default_parameters()
)
ptr_e_sc_clmp_true <- data.frame(mr_e_sc_clmp_true["Q"])
egger_e_sc_clmp_true <- mr_egger_regression(
  true_sc_e_clmp$beta.exposure,
  true_sc_e_clmp$beta.outcome,
  true_sc_e_clmp$se.exposure,
  true_sc_e_clmp$se.outcome,
  parameters
)
ptr_e_sc_clmp_true[2, 1] <- egger_e_sc_clmp_true["Q"]
F <- abs(true_sc_e_clmp$beta.exposure)^2 / true_sc_e_clmp$se.exposure^2
mF <- mean(F)
ptr_e_sc_clmp_true[3, 1] <- mF

#MR-PRESSO
sc_e_mdd__true_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                             SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                             OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = true_sc_e_clmp, NbDistribution = 1000,  SignifThreshold = 0.05)

sc_e_mdd_true_presso_mainres <- sc_e_mdd__true_presso$`Main MR results`
sc_e_mdd_true_presso_mainres <- sc_e_mdd_true_presso_mainres %>%
  rename(b = `Causal Estimate`)
sc_e_mdd_true_presso_mainres <- sc_e_mdd_true_presso_mainres %>%
  rename(se = `Sd`)
sc_e_mdd_true_presso_mainres <- generate_odds_ratios(sc_e_mdd_true_presso_mainres)
sc_e_mdd_true_presso_mainres$Distortion_P <- sc_e_mdd__true_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
sc_e_mdd__true_presso$`MR-PRESSO results`$`Outlier Test`$nOUT <- ifelse(sc_e_mdd__true_presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, TRUE, FALSE)
summary(sc_e_mdd__true_presso$`MR-PRESSO results`$`Outlier Test`$nOUT) #nOUT TRUE = 0 SNPs  
ptr_e_sc_clmp_true[4, 1] <- NA

##Ever smokers
#Calculate the regression dilution Isq. 
BetaXG   = true_sc_e_clmp$beta.exposure
seBetaXG = true_sc_e_clmp$se.exposure 
seBetaYG <- true_sc_e_clmp$se.outcome

BXG  = abs(BetaXG)         # gene--exposure estimates are positive  

Isq_unweighted_sc_e_true <- ISQ(BXG,seBetaXG) #unweighted = 0.66
Isq_weighted_sc_e_true <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0.53

## SIMEX Correction - Ever Smokers
# If the Isquared (weighted or unweighted) is between 0.6 and 0.9, conduct a SIMEX correction
BetaXG   = true_sc_e_clmp$beta.exposure
BetaYG <- true_sc_e_clmp$beta.outcome
seBetaXG = true_sc_e_clmp$se.exposure
seBetaYG <- true_sc_e_clmp$se.outcome

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
simex_sc_e_true <- data.frame(mod2$coefficients)
simex_sc_e_true <- simex_sc_e_true %>%
  rename(b = `jackknife.Estimate`, se = 'jackknife.Std..Error', pval = 'jackknife.Pr...t..')
simex_sc_e_true <- generate_odds_ratios(simex_sc_e_true)

#Save results
setwd(" ")
write.csv(simex_sc_e_true, "simex_sc_ever_steiger.csv", row.names=FALSE, quote=FALSE)

#################################################
#########   Binding results tables   ############
#################################################

##Ever
#Add other results into main results data frame (i.e., SIMEX + PRESSO)
#SIMEX:
result_sc_ever_clumped_true[5, 1] <- result_sc_ever_clumped_true[1, "id.exposure"]
result_sc_ever_clumped_true[5,2] <- result_sc_ever_clumped_true[1, "id.outcome"]
result_sc_ever_clumped_true[5,3] <- result_sc_ever_clumped_true[1, "outcome"]
result_sc_ever_clumped_true[5,4] <- result_sc_ever_clumped_true[1, "exposure"]
result_sc_ever_clumped_true[5,5] <- "MR Egger (SIMEX)"
result_sc_ever_clumped_true[5,6] <- result_sc_ever_clumped_true[1, "nsnp"]
result_sc_ever_clumped_true[5,7] <- simex_sc_e_true[2, "b"]
result_sc_ever_clumped_true[5,8] <- simex_sc_e_true[2, "se"]
result_sc_ever_clumped_true[5,9] <- simex_sc_e_true[2, "pval"]
result_sc_ever_clumped_true[5,10] <- simex_sc_e_true[2, "lo_ci"]
result_sc_ever_clumped_true[5,11] <- simex_sc_e_true[2, "up_ci"]
result_sc_ever_clumped_true[5,12] <- simex_sc_e_true[2, "or"]
result_sc_ever_clumped_true[5,13] <- simex_sc_e_true[2, "or_lci95"]
result_sc_ever_clumped_true[5,14] <- simex_sc_e_true[2, "or_uci95"]
#PRESSO:
#No outliers identified so this is the "raw" estimate (i.e., not outlier corrected)
result_sc_ever_clumped_true[6, 1] <- result_sc_ever_clumped_true[1, "id.exposure"]
result_sc_ever_clumped_true[6,2] <- result_sc_ever_clumped_true[1, "id.outcome"]
result_sc_ever_clumped_true[6,3] <- result_sc_ever_clumped_true[1, "outcome"]
result_sc_ever_clumped_true[6,4] <- result_sc_ever_clumped_true[1, "exposure"]
result_sc_ever_clumped_true[6,5] <- "MR PRESSO"
result_sc_ever_clumped_true[6,6] <- result_sc_ever_clumped_true[1, "nsnp"] 
result_sc_ever_clumped_true[6,7] <- sc_e_mdd_true_presso_mainres[1, "b"]
result_sc_ever_clumped_true[6,8] <- sc_e_mdd_true_presso_mainres[1, "se"]
result_sc_ever_clumped_true[6,9] <- sc_e_mdd_true_presso_mainres[1, "P-value"]
result_sc_ever_clumped_true[6,10] <- sc_e_mdd_true_presso_mainres[1, "lo_ci"]
result_sc_ever_clumped_true[6,11] <- sc_e_mdd_true_presso_mainres[1, "up_ci"]
result_sc_ever_clumped_true[6,12] <- sc_e_mdd_true_presso_mainres[1, "or"]
result_sc_ever_clumped_true[6,13] <- sc_e_mdd_true_presso_mainres[1, "or_lci95"]
result_sc_ever_clumped_true[6,14] <- sc_e_mdd_true_presso_mainres[1, "or_uci95"]
#Remove MR-Egger
result_sc_ever_clumped_true <- subset(result_sc_ever_clumped_true, method != "MR Egger")

#Save tables
setwd(" ")
write.csv(result_sc_ever_clumped_true, "results_sc_ever_steiger.csv", row.names=FALSE, quote=FALSE)

################################
######### Plot results #########
################################
#https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#plots 

setwd("//rdsfcifs.acrc.bris.ac.uk/UK_Biobank_Project_9142/Working folder Chloe/Results/2sMR Tobacco Cannabis")
#Renaming exposure and outcome for plots
outcome_mdd_dat_sc_never_clumped$exposure <- gsub("exposure", "Smoking Cessation", outcome_mdd_dat_sc_never_clumped$exposure)
outcome_mdd_dat_sc_never_clumped$outcome<- gsub("outcome", "MDD", outcome_mdd_dat_sc_never_clumped$outcome)
outcome_mdd_dat_sc_ever_clumped$exposure <- gsub("exposure", "Smoking Cessation", outcome_mdd_dat_sc_ever_clumped$exposure)
outcome_mdd_dat_sc_ever_clumped$outcome<- gsub("outcome", "MDD", outcome_mdd_dat_sc_ever_clumped$outcome)

# single SNP analyses
sc_never_clmp_single <- mr_singlesnp(outcome_mdd_dat_sc_never_clumped, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
sc_never_clmp_single <- mr_forest_plot(sc_never_clmp_single)
ggsave(sc_never_clmp_single[[1]], file="single_forest_sc_never_clumped.png", width=7, height=10)

# single SNP analyses
sc_ever_clmp_single <- mr_singlesnp(outcome_mdd_dat_sc_ever_clumped, all_method=c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
sc_ever_clmp_single <- mr_forest_plot(sc_ever_clmp_single)
ggsave(sc_ever_clmp_single[[1]], file="single_forest_sc_ever_clumped.png", width=7, height=10)

# leave one out analyses
sc_never_clmp_loo <- mr_leaveoneout(outcome_mdd_dat_sc_never_clumped)
sc_never_clmp_loo <- mr_leaveoneout_plot(sc_never_clmp_loo)
ggsave(sc_never_clmp_loo[[1]], file="loo_sc_never_clumped.png", width=7, height=8)

# leave one out analyses
sc_ever_clmp_loo <- mr_leaveoneout(outcome_mdd_dat_sc_ever_clumped)
sc_ever_clmp_loo <- mr_leaveoneout_plot(sc_ever_clmp_loo)
ggsave(sc_ever_clmp_loo[[1]], file="loo_sc_ever_clumped.png", width=7, height=8)

# scatter plot
plot_result_sc_never_clmp <- subset(result_sc_never_clumped, method != "MR PRESSO")
scatter_sc_never_clmp <- mr_scatter_plot(plot_result_sc_never_clmp, outcome_mdd_dat_sc_never_clumped)
ggsave(scatter_sc_never_clmp[[1]], file="scatter_sc_never_clumped.png", width=7, height=7)

# scatter plot
scatter_sc_ever_clmp <- mr_scatter_plot(result_sc_ever_clumped, outcome_mdd_dat_sc_ever_clumped)
ggsave(scatter_sc_ever_clmp[[1]], file="scatter_sc_ever_clumped.png", width=7, height=7)

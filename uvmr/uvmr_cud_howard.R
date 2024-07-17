# Script was created by Chloe Burke April 2024.
# The script conducts univariable MR exploring the effects of cannabis use disorder on lifetime depression.
# The SNPs used in this script were selected based on the findings from a multi-ancestry GWAS of CUD (Levey et al., 2023)
# The outcome used in this script is based on a GWAS of major depression by Howard et al., (2018), with UK Biobank included

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

#################################
##### Set SNP list for CUD  ##### 
#################################
setwd(" ")
CUD_SNPlist<- read.table("CUD_SNPlist.txt", header=FALSE)
CUD_SNPlist<-rename(CUD_SNPlist, c("SNP" = "V1"))

##############################################################
##### Extract exposure data for MR of CUD and depression #####
##############################################################
# EAF missing from summary statistics
# Contacted lead author via e-mail who confirmed these are not available
# Stated that "MAF values for the vast majority of the samples (those we do have MAF for) are virtually identical to 1KG"

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

##Adding information on EAF from 1000G ref (EUR; e.g., https://www.ensembl.org/Multi/Search/Results?q=rs10986600;site=ensembl_all)
setwd(" ")
cud_eaf <- read.csv("CUD_eaf_1000G.csv", header = T)
cud_dat_mr_clumped_eaf <- merge(cud_dat_mr_clumped, cud_eaf, by = "SNP", all = TRUE)
cud_dat_mr_clumped_eaf <- subset (cud_dat_mr_clumped_eaf, select = -eaf.exposure.x)
cud_dat_mr_clumped_eaf <- rename(cud_dat_mr_clumped_eaf, eaf.exposure = eaf.exposure.y)
cud_dat_mr_clumped_eaf <- select(cud_dat_mr_clumped_eaf, SNP, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure, eaf.exposure, exposure, mr_keep.exposure, pval_origin.exposure, id.exposure)

#Checking formatting of file
unique(cud_dat_mr_clumped_eaf$eaf.exposure)

#Removing trailing spaces
cud_dat_mr_clumped_eaf$eaf.exposure <- str_trim(cud_dat_mr_clumped_eaf$eaf.exposure)

#Convert to numeric
cud_dat_mr_clumped_eaf$eaf.exposure <- as.numeric(cud_dat_mr_clumped_eaf$eaf.exposure)

#######################################
##### Extract outcome data for MR #####
#######################################
setwd(" ")

howard_dat_cud_clumped <- read_outcome_data(
  "howard_gwas_fixed.txt",
  snps = cud_dat_mr_clumped_eaf$SNP,
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
)

##################
## Saving files ##
##################
setwd(" ")
write.csv(cud_dat_mr_clumped_eaf,"cud_dat_uvmr.csv", row.names = FALSE)
write.csv(howard_dat_cud_clumped,"cud_howard_uvmr.csv", row.names = FALSE)

#################################
#####   Convert OR to log   #####
#################################
#Note. SE already for log odds
howard_dat_cud_clumped$beta.outcome<-log(howard_dat_cud_clumped$beta.outcome)

#######################
##### Harmonising #####
#######################
#Removing the following SNPs for being palindromic with intermediate allele frequencies:  
#rs56070621
outcome_howard_dat_cud_clumped <-harmonise_data(cud_dat_mr_clumped_eaf, howard_dat_cud_clumped, action = 2)

##################
## Saving files ##
##################
setwd(" ")
write.csv(outcome_howard_dat_cud_clumped,"outcome_howard_cud_uvmr.csv", row.names = FALSE)

###################
## Reading files ##
####################
#setwd(" ")
#outcome_howard_dat_cud_clumped <- read.csv("outcome_howard_cud_uvmr.csv", header = T)

################################################################################
################################## MR ##########################################
################################################################################

#IVW, Egger, Weighted Median + Weighted Mode
result_howard_cud_clumped <- mr(
  outcome_howard_dat_cud_clumped,
  method_list = c("mr_ivw", "mr_weighted_median", "mr_weighted_mode", "mr_egger_regression")
)
result_howard_cud_clumped <- generate_odds_ratios(result_howard_cud_clumped)
outcome_howard_dat_cud_clumped <- subset(outcome_howard_dat_cud_clumped, mr_keep)
mr_how_cud_clumped <- TwoSampleMR::mr_ivw(
  outcome_howard_dat_cud_clumped$beta.exposure,
  outcome_howard_dat_cud_clumped$beta.outcome,
  outcome_howard_dat_cud_clumped$se.exposure,
  outcome_howard_dat_cud_clumped$se.outcome,
  parameters = default_parameters()
)
ptr_how_cud_clumped <- data.frame(mr_how_cud_clumped["Q"])
egger_how_cud_clumped <- mr_egger_regression(
  outcome_howard_dat_cud_clumped$beta.exposure,
  outcome_howard_dat_cud_clumped$beta.outcome,
  outcome_howard_dat_cud_clumped$se.exposure,
  outcome_howard_dat_cud_clumped$se.outcome,
  parameters
)
ptr_how_cud_clumped[2, 1] <- egger_how_cud_clumped ["Q"]
F <- outcome_howard_dat_cud_clumped$beta.exposure^2 / outcome_howard_dat_cud_clumped$se.exposure^2
mF <- mean(F)
ptr_how_cud_clumped [3, 1] <- mF

#MR-PRESSO
cud_how_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                           SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                           OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = outcome_howard_dat_cud_clumped, NbDistribution = 1000,  SignifThreshold = 0.05)

cud_how_presso_mainres <- cud_how_presso$`Main MR results`
cud_how_presso_mainres <- cud_how_presso_mainres %>%
  rename(b = `Causal Estimate`)
cud_how_presso_mainres <- cud_how_presso_mainres %>%
  rename(se = `Sd`)
cud_how_presso_mainres <- generate_odds_ratios(cud_how_presso_mainres)
cud_how_presso_mainres$Distortion_P <- cud_how_presso$`MR-PRESSO results`$`Distortion Test`$Pvalue
cud_how_presso$`MR-PRESSO results`$`Outlier Test`$nOUT <- ifelse(cud_how_presso$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05, TRUE, FALSE)
summary(cud_how_presso$`MR-PRESSO results`$`Outlier Test`$nOUT) #No outliers identified
ptr_how_cud_clumped [4, 1] <- NA

#Save results
setwd(" ")
write.csv(cud_how_presso_mainres, "presso_cud_howard_uvmr.csv", row.names=FALSE, quote=FALSE)

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
BetaXG   = outcome_howard_dat_cud_clumped$beta.exposure
seBetaXG = outcome_howard_dat_cud_clumped$se.exposure 
seBetaYG <- outcome_howard_dat_cud_clumped$se.outcome

BXG  = abs(BetaXG)         # gene--exposure estimates are positive  

Isq_unweighted_cud_how <- ISQ(BXG,seBetaXG) #unweighted = 0.52
Isq_weighted_cud_how <- ISQ((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted = 0

#I2GX <0.6, not approrpiate to conduct SIMEX correction or report MR-Egger
#################################################
#########   Binding results tables   ############
#################################################

#Add other results into main results data frame (i.e., PRESSO)
#PRESSO:
#No outliers identified so this is the "raw" estimate (i.e., not outlier corrected)
result_howard_cud_clumped[5, 1] <- result_howard_cud_clumped[1, "id.exposure"]
result_howard_cud_clumped[5,2] <- result_howard_cud_clumped[1, "id.outcome"]
result_howard_cud_clumped[5,3] <- result_howard_cud_clumped[1, "outcome"]
result_howard_cud_clumped[5,4] <- result_howard_cud_clumped[1, "exposure"]
result_howard_cud_clumped[5,5] <- "MR PRESSO"
result_howard_cud_clumped[5,6] <- result_howard_cud_clumped[1, "nsnp"] #No outliers removed so same as other analyses
result_howard_cud_clumped[5,7] <- cud_how_presso_mainres[1, "b"]
result_howard_cud_clumped[5,8] <- cud_how_presso_mainres[1, "se"]
result_howard_cud_clumped[5,9] <- cud_how_presso_mainres[1, "P-value"]
result_howard_cud_clumped[5,10] <- cud_how_presso_mainres[1, "lo_ci"]
result_howard_cud_clumped[5,11] <- cud_how_presso_mainres[1, "up_ci"]
result_howard_cud_clumped[5,12] <- cud_how_presso_mainres[1, "or"]
result_howard_cud_clumped[5,13] <- cud_how_presso_mainres[1, "or_lci95"]
result_howard_cud_clumped[5,14] <- cud_how_presso_mainres[1, "or_uci95"]
#Remove MR-Egger
result_howard_cud_clumped <- subset(result_howard_cud_clumped, method != "MR Egger")

#Save results
setwd(" ")
write.csv(result_howard_cud_clumped, "results_cud_howard_uvmr.csv", row.names=FALSE, quote=FALSE)

#####################################################
#########   Test for reverse causation   ############
#####################################################
### Steiger filtering ###
# Binary exposures/outcomes must contain: samplesize, ncase, ncontrol, prevalence and units="log odds"
# For cannabis use disorder there were 886,025; nCASE = 42,281, nCTRL = 843,744 taken from Table 1 in https://www.nature.com/articles/s41588-023-01563-z#Sec2 
# With iPSYCH2 removed this is 785635, nCASE = 37548, nCTRL = 748087
# For MDD (-23&Me) there were nCASE = 127,552 and nCTRL = 233,763 i.e., total sample = 361315; Table 1 in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6522363/
outcome_howard_dat_cud_clumped$samplesize.exposure <- 785635
outcome_howard_dat_cud_clumped$ncase.exposure <-37548
outcome_howard_dat_cud_clumped$ncontrol.exposure <- 748087
outcome_howard_dat_cud_clumped$samplesize.outcome <- 361315
outcome_howard_dat_cud_clumped$ncase.outcome <-127552
outcome_howard_dat_cud_clumped$ncontrol.outcome <- 233763

# Adding a unit column labelled "log odds" and add a prevalence column
outcome_howard_dat_cud_clumped$prevalence.exposure <-0.05
outcome_howard_dat_cud_clumped$prevalence.outcome <-0.35
outcome_howard_dat_cud_clumped$units.exposure <-"log odds"
outcome_howard_dat_cud_clumped$units.outcome <- "log odds"

#Run Steiger filtering
steiger <- steiger_filtering(outcome_howard_dat_cud_clumped)

#How many SNPs explain more variance in the exposure than in the outcome?
#Generates steiger_dir (which is `TRUE` if the rsq.exposure is larger than rsq.outcome)
table(steiger$steiger_dir) #0 are FALSE
true_cud_howard <-subset(steiger, steiger$steiger_dir==TRUE)
str(true_cud_howard)

################################
######### Plot results #########
################################
#https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html#plots 

setwd(" ")
#Renaming exposure and outcome for plots
outcome_howard_dat_cud_clumped$exposure <- gsub("exposure", "CUD", outcome_howard_dat_cud_clumped$exposure)
outcome_howard_dat_cud_clumped$outcome<- gsub("outcome", "MDD (Howard)", outcome_howard_dat_cud_clumped$outcome)

# single SNP analyses
cud_howard_single <- mr_singlesnp(outcome_howard_dat_cud_clumped, all_method=c("mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
cud_howard_single <- mr_forest_plot(cud_howard_single)
ggsave(cud_howard_single[[1]], file="single_forest_cud_howard.png", width=7, height=10)

# leave one out analyses
cud_howard_loo <- mr_leaveoneout(outcome_howard_dat_cud_clumped)
cud_howard_loo <- mr_leaveoneout_plot(cud_howard_loo)
ggsave(cud_howard_loo[[1]], file="loo_cud_howard.png", width=7, height=10)

# scatter plot
scatter_cud_howard <- mr_scatter_plot(result_howard_cud_clumped, outcome_howard_dat_cud_clumped)
ggsave(scatter_cud_howard[[1]], file="scatter_cud_howard.png", width=7, height=7)

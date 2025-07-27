# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# We have the individual and cluster (facility) level network measures (see 
# attached). But Raj had mentioned group level (i.e., the groups within the 
# clusters/facilities) network measures.

# Here is the summary we want to do (to check if I understand it correctly):
#   Do (pairwise) mediation analysis with mediation packages
# UC vs GMV
# UC vs MF
# UC vs GMVMF
# GMV vs GMVMF
# Do (pairwise) regression (and include 95% CI) where M=Mediator (at 12 month) and Y=SBP change
# M ~ treatment + confounders
# Y ~ treatment + confounders
# Y ~ M + treatment + confounders

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Load R packages
# library(MCMCpack)
# library(miscTools)
# library(mvnfast)
# library(sn)
library(MASS)
library(Rcpp)
library(RcppArmadillo)

library(dplyr)
library(tidyr)
library(mgcv)
library(mediation)
library(lme4)

# library(devtools)
# devtools::install_github("RcppCore/RcppArmadillo")

# Sys.setenv(MAKEFLAGS = paste0("-j",parallel::detectCores()))
# install.packages(c("StanHeaders","rstan"),type="source")
# # https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows
# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(USE_CXX14 = 1)

# library(dbarts)
library(rBARTmediation)

# Extract ID for simulated dataset (specific to LSF computing cluster)
# Note: The LSB_JOBINDEX is specified in the bsub command using the -J
# option
# run_ID = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

set.seed(0)

esttype = "mean" # "median"
saveall = TRUE # FALSE

# the number of MC integration per iteration for a MCMC chain
# gibbs_thin = 5e1 # *gibbs_GLM
# gibbs_iter = 1e4 # *gibbs_GLM
# gibbs_burnin = 5e4 # *gibbs_GLM
gibbs_thin = 1e1
gibbs_iter = 1e2
gibbs_burnin = 1e3

Mdist = "Continuous"
# Mdist = "Binary"
# Mdist = "Ordinal";K = 7 # the number of ordinal mediator when Mdist == "Ordinal"
# Mdist = "Count"

Ydist = "Continuous"
# Ydist = "Binary"
# Ydist = "Ordinal";K = 7 # the number of ordinal mediator when Ydist == "Ordinal"
# Ydist = "Count"

# sparse = TRUE
sparse = FALSE
# ntree = 20
# ntree = 50
# ntree = 100
# ntree = 200
ntree = 500

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# This code must be in your current directory or you can change the path.
setwd("/Users/WooJung/Documents/Rproject/TSBBmediation/source")

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Load R code
source("TSBBmediation_source_r.R")
# source("TSBBmediation_source_rBIGPIC.R")
# Load cpp code
sourceCpp("TSBBmediation_source_cpp.cpp")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# This code must be in your current directory or you can change the path.
setwd("/Users/WooJung/Documents/Rproject/TSBBmediation/RealDataAnalysis/DATA/BIGPICdata")
# Load Data
load("BIGPICdata.Rdata")

# length(BIGPICdata$bigpic_group_id)
# sum(is.na(BIGPICdata$bigpic_group_id))
# sum(Z==0)
# sum(Z==1)
# sum(Z==2)
# sum(Z==3)

# colnames(BIGPICdata)
# 
# AAA = BIGPICdata %>%
#   select(site, study_arm) %>% 
#   distinct(site, study_arm, .keep_all = TRUE)
# 
# BBB = BIGPICdata %>%
#   select(bigpic_group_id, site, study_arm) %>% 
#   distinct(bigpic_group_id, site, study_arm, .keep_all = TRUE)
# 
# CCC = BIGPICdata %>%
#   select(participant_id, bigpic_group_id, site, study_arm) %>% 
#   distinct(participant_id, bigpic_group_id, site, study_arm, .keep_all = TRUE)
# View(CCC)

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Clusters: j=1,...,J=24
Uindex = BIGPICdata$site
sitetable = table(Uindex);sitetable
sitenames = names(sitetable)
J = length(sitenames)

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# treatment: 
# Z=0 (UC: Usual Care)
# Z=1 (UCMF: Usual Care + Microfinance)
# Z=2 (GMV: Group Medical Visit)
# Z=3 (GMVMF: Group Medical Visit + Microfinance)
# summariessss(BIGPICdata$arm.refUC2)
Z = BIGPICdata$arm.refUC2
ZZ = numeric(length(Z))
ZZ[which(Z=="UC")] = 0
ZZ[which(Z=="UCMF")] = 1
ZZ[which(Z=="GMV")] = 2
ZZ[which(Z=="GMVMF")] = 3
Z = ZZ

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# outcome: the 1-year change in systolic blood pressure (SBP).
Y = BIGPICdata$sbp.change.12mos

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# individual-level
# BIGPICdata$bl_sbp_ave # baseline SBP
# ifelse(BIGPICdata$gender=="Male",1,0) # sex
# BIGPICdata$age_est # age
# BIGPICdata$pathway # pathway

# cluster-level
# BIGPICdata$HF # HF 
# BIGPICdata$Gise_Penetration # GISE Penetration
# sapply(1:J, function(j) mean(dataYZCV$C7[which(dataYZCV$Uindex==j)],na.rm=T)) # cluster-specific means of SBP at baseline

# will add up later
# baseline values of the mediator (M). 

# potential individual-level confounder
# age
C1 = BIGPICdata$age_est

# gender
C2 = ifelse(BIGPICdata$gender=="Male",1,0)

# Livestock (Cows, Sheep, Goats, Donkeys, etc) : None or Some 
C3 = ifelse(BIGPICdata$livestock.cat=="Some",1,0)

# binary: Unemployed or Employed
C4 = ifelse(BIGPICdata$earnings_status=="Employed",1,0)

# binary: Derived variable - baseline blood sugars
# "Yes": (dat$bl_bsg_type=="Random"  & dat$bl_bsg>=11.1) or
# (dat$bl_bsg_type=="Fasting"  & dat$bl_bsg >= 7) or 
# (dat$bl_ever_diagnosis_1=="Yes") or
# (dat$bl_ever_diagnosis_2=="Yes")
# otherwise "NO"
C5 = ifelse(BIGPICdata$dm=="Yes",1,0)

C6 = BIGPICdata$bl_sbp_ave

# binary: elevated BP (SBP ≥140 or diastolic BP (DBP) ≥ 90)
C7 = ifelse(BIGPICdata$elavated_bp=="Sbp/Dbp >= 140/90",1,0)

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# potential cluster-level confounder
# BIGPIC County: 
# "Busia", "Kisumu", "Uasin Gishu", "Trans-Nzoia"
factorV1 = factor(BIGPICdata$county)
Kv1 = length(levels(factorV1))
V1 = sapply(2:Kv1, function(l) ifelse(factorV1==levels(factorV1)[l],1,0))

# HF: Prerandomization variable that tell the type/leve of health facility near the community
# Dis: Dispensary
# HC: health center 
# SCH: Sub county hospital
factorV2 = factor(BIGPICdata$HF)
Kv2 = length(levels(factorV2))
V2 = sapply(2:Kv2, function(l) ifelse(factorV2==levels(factorV2)[l],1,0))

# GISE (groups intergrateed saving for empowerment)
# Prerandomization variable that tell the near the level of GISE penetrtaion around the community. 
# Light & Heavy
V3 = ifelse(BIGPICdata$Gise_Penetration=="Heavy",1,0)

V4_temp = sapply(1:J, function(l) mean(C6[which(Uindex==sitenames[l])],na.rm=TRUE))
V4 = numeric(length(Y))
for (l in 1:J) {
  whichsite = which(Uindex==sitenames[l])
  V4[whichsite] = V4_temp[l]
  C6[whichsite] = C6[whichsite] - V4_temp[l]
}

# -----------------------------------------------------------------------
dataYZCV = unique(data.frame(participant_id = BIGPICdata$participant_id, 
                             Y=Y, 
                             Z=Z, 
                             C1=C1, 
                             C2=C2, 
                             C3=C3, 
                             C4=C4, 
                             C5=C5, 
                             C6=C6, 
                             # C7=C7,
                             V1=V1, 
                             V2=V2, 
                             V3=V3, 
                             V4=V4,
                             Uindex = as.numeric(Uindex)))

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# This code must be in your current directory or you can change the path.
setwd("/Users/WooJung/Documents/Rproject/TSBBmediation/RealDataAnalysis/DATA")

# mediator (intermediate factors):
# Individual Level Measures
indBIGPICdata = read.csv("ind_meas_2021_031721.csv")
names(indBIGPICdata)

# # pick the data with Wave = 2 (3 month)
# wave1which = which(indBIGPICdata$wave == 1) # month 0
# wave2which = which(indBIGPICdata$wave == 2) # month 3
# wave3which = which(indBIGPICdata$wave == 3) # month 12
# indBIGPICdataWAVE1 = indBIGPICdata[wave1which,]
# indBIGPICdataWAVE2 = indBIGPICdata[wave2which,]
# indBIGPICdataWAVE3 = indBIGPICdata[wave3which,]
# 
# # Given R = 3 (relation = CGC Network)
# R3wave1which = which(indBIGPICdataWAVE1$relation==3)
# R3wave2which = which(indBIGPICdataWAVE2$relation==3)
# R3wave3which = which(indBIGPICdataWAVE3$relation==3)
# 
# indBIGPICdataWAVE1R3 = indBIGPICdataWAVE1[R3wave1which,]
# indBIGPICdataWAVE2R3 = indBIGPICdataWAVE2[R3wave2which,]
# indBIGPICdataWAVE3R3 = indBIGPICdataWAVE3[R3wave3which,]

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# This code must be in your current directory or you can change the path.
setwd("/Users/WooJung/Documents/Rproject/TSBBmediation/RealDataAnalysis/DATA")

# mediator (intermediate factors): 
# Network Level Measures
netBIGPICdata = unique(read.csv("net_level_2021.csv"))
names(netBIGPICdata)

# # pick the data with Wave = 2 (month 3) or Wave = 3 (month 12)
# wave1which = which(netBIGPICdata$wave == 1) # month 0
# wave2which = which(netBIGPICdata$wave == 2) # month 3
# wave3which = which(netBIGPICdata$wave == 3) # month 12
# netBIGPICdataWAVE1 = netBIGPICdata[wave1which,]
# netBIGPICdataWAVE2 = netBIGPICdata[wave2which,]
# netBIGPICdataWAVE3 = netBIGPICdata[wave3which,]

# -----------------------------------------------------------------------
indkeep = c("participant_id","facility_id","relation")
netkeep = c("facility_id", "relation", "density", "transitivity", "cohesion", "APL", "centralization")

# -----------------------------------------------------------------------
netindBIGPICdata <- indBIGPICdata %>%
  select(participant_id, facility_id) %>%
  inner_join(netBIGPICdata %>%
               select(facility_id, wave,relation,density,transitivity,cohesion,APL,centralization), 
             by = "facility_id", 
             relationship = "many-to-many") %>%
  distinct(participant_id,facility_id,wave,relation,density,transitivity,cohesion,APL,centralization, .keep_all = TRUE) %>%
  arrange(participant_id,facility_id,wave)

dataYZCVnetindBIGPICdata <- netindBIGPICdata %>%
  inner_join(dataYZCV, by = "participant_id") %>%
  select(participant_id,facility_id,wave,relation,density,transitivity,cohesion,APL,centralization,Uindex,
         Y,Z,C1,C2,C3 ,C4,C5,C6,V1.1,V1.2,V1.3,V2.1,V2.2,V3,V4) %>% # C7,
  arrange(participant_id,wave,Uindex) %>%
  drop_na()

NEWdataYZCVnetindBIGPICdata = dataYZCVnetindBIGPICdata %>%
  filter(Z %in% c(2, 3)) %>%
  filter(wave %in% c(3)) %>% # %in% c(1,  3) ==> using c(3)
  filter(relation %in% c(2)) %>%
  arrange(Uindex,participant_id,wave,facility_id,Z,relation) %>%
  drop_na()
NEWdataYZCVnetindBIGPICdata$Z = ifelse(NEWdataYZCVnetindBIGPICdata$Z==2, 0, 1)
NEWdataYZCVnetindBIGPICdata$Uindex = as.factor(NEWdataYZCVnetindBIGPICdata$Uindex)
# dim(NEWdataYZCVnetindBIGPICdata)
# length(unique(NEWdataYZCVnetindBIGPICdata$participant_id))
# table(NEWdataYZCVnetindBIGPICdata$Uindex, NEWdataYZCVnetindBIGPICdata$Z)

Uindex = NEWdataYZCVnetindBIGPICdata$Uindex
uniqueUindex = sort(unique(Uindex))
J = length(uniqueUindex)
NEWdataYZCVnetindBIGPICdata$Uindex = apply(sapply(1:J, function(l) ifelse(Uindex==uniqueUindex[l],l,0)),1,sum)

NEWdataYZCVnetindBIGPICdata = NEWdataYZCVnetindBIGPICdata %>%
  # Join the data with a summarized version of itself
  left_join(
    # This inner part creates the summary table on the fly
    NEWdataYZCVnetindBIGPICdata %>%
      group_by(Uindex) %>%
      summarise(across(C1:C6, ~mean(., na.rm = TRUE), .names = "bar{.col}")),
    # Specify the column to join by
    by = "Uindex"
  )

# View(NEWdataYZCVnetindBIGPICdata)
# View(NEWdataYZCVnetindBIGPICdata[which(NEWdataYZCVnetindBIGPICdata$Uindex==21),])

Y = NEWdataYZCVnetindBIGPICdata$Y
Z = NEWdataYZCVnetindBIGPICdata$Z
Uindex = NEWdataYZCVnetindBIGPICdata$Uindex
uniqueUindex = sort(unique(Uindex))
J = length(uniqueUindex)
Uindex = apply(sapply(1:J, function(l) ifelse(Uindex==uniqueUindex[l],l,0)),1,sum)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
M_NEWdataYZCVnetindBIGPICdata = NEWdataYZCVnetindBIGPICdata %>%
  select(
    wave, relation, density, transitivity, cohesion, APL, centralization, Uindex, Z,
    starts_with("V"),
    starts_with("barC")
  )  %>%
  arrange(Uindex) %>%
  distinct()
names(M_NEWdataYZCVnetindBIGPICdata)
names(NEWdataYZCVnetindBIGPICdata)

run_ID = 1
Y = NEWdataYZCVnetindBIGPICdata$Y
if (run_ID == 1) {
  M = M_NEWdataYZCVnetindBIGPICdata$density
} else if (run_ID == 2) {
  M = M_NEWdataYZCVnetindBIGPICdata$transitivity
} else if (run_ID == 3) {
  M = M_NEWdataYZCVnetindBIGPICdata$cohesion
} else if (run_ID == 4) {
  M = M_NEWdataYZCVnetindBIGPICdata$APL
}
Z = NEWdataYZCVnetindBIGPICdata$Z
Uindex = NEWdataYZCVnetindBIGPICdata$Uindex
uniqueUindex = sort(unique(Uindex))
J = length(uniqueUindex)
Uindex = apply(sapply(1:J, function(l) ifelse(Uindex==uniqueUindex[l],l,0)),1,sum)

matX = cbind(M_NEWdataYZCVnetindBIGPICdata$Z,
             M_NEWdataYZCVnetindBIGPICdata$V1.1,
             M_NEWdataYZCVnetindBIGPICdata$V1.2,
             M_NEWdataYZCVnetindBIGPICdata$V1.3,
             M_NEWdataYZCVnetindBIGPICdata$V2.1,
             M_NEWdataYZCVnetindBIGPICdata$V2.2,
             M_NEWdataYZCVnetindBIGPICdata$V3,
             M_NEWdataYZCVnetindBIGPICdata$V4,
             M_NEWdataYZCVnetindBIGPICdata$barC1,
             M_NEWdataYZCVnetindBIGPICdata$barC2,
             M_NEWdataYZCVnetindBIGPICdata$barC3,
             M_NEWdataYZCVnetindBIGPICdata$barC4,
             M_NEWdataYZCVnetindBIGPICdata$barC5,
             M_NEWdataYZCVnetindBIGPICdata$barC6)
BARTfitM = BART(M, matX, typeY = 'continuous', 
                ndpost=gibbs_iter, nskip=gibbs_burnin, keepevery=gibbs_thin, 
                ntree=ntree, sparse = sparse)
# predict(BARTfitM, matX)

matM = cbind(NEWdataYZCVnetindBIGPICdata$Z,
             NEWdataYZCVnetindBIGPICdata$transitivity,
             NEWdataYZCVnetindBIGPICdata$V1.1,
             NEWdataYZCVnetindBIGPICdata$V1.2,
             NEWdataYZCVnetindBIGPICdata$V1.3,
             NEWdataYZCVnetindBIGPICdata$V2.1,
             NEWdataYZCVnetindBIGPICdata$V2.2,
             NEWdataYZCVnetindBIGPICdata$V3,
             NEWdataYZCVnetindBIGPICdata$V4,
             NEWdataYZCVnetindBIGPICdata$C1,
             NEWdataYZCVnetindBIGPICdata$C2,
             NEWdataYZCVnetindBIGPICdata$C3,
             NEWdataYZCVnetindBIGPICdata$C4,
             NEWdataYZCVnetindBIGPICdata$C5,
             NEWdataYZCVnetindBIGPICdata$C6)
BARTfitY = rBART(NEWdataYZCVnetindBIGPICdata$Y, matM, NEWdataYZCVnetindBIGPICdata$Uindex, typeY = "continuous", 
                 ndpost=gibbs_iter, nskip=gibbs_burnin, keepevery=gibbs_thin,
                 ntree=ntree, sparse = sparse)
# predict(BARTfitY, matM, Uindex)

Vjmat= cbind(M_NEWdataYZCVnetindBIGPICdata$V1.1,
                M_NEWdataYZCVnetindBIGPICdata$V1.2,
                M_NEWdataYZCVnetindBIGPICdata$V1.3,
                M_NEWdataYZCVnetindBIGPICdata$V2.1,
                M_NEWdataYZCVnetindBIGPICdata$V2.2,
                M_NEWdataYZCVnetindBIGPICdata$V3,
                M_NEWdataYZCVnetindBIGPICdata$V4)

Cjmat = cbind(M_NEWdataYZCVnetindBIGPICdata$barC1,
                M_NEWdataYZCVnetindBIGPICdata$barC2,
                M_NEWdataYZCVnetindBIGPICdata$barC3,
                M_NEWdataYZCVnetindBIGPICdata$barC4,
                M_NEWdataYZCVnetindBIGPICdata$barC5,
                M_NEWdataYZCVnetindBIGPICdata$barC6)

Vmat = cbind(NEWdataYZCVnetindBIGPICdata$V1.1,
               NEWdataYZCVnetindBIGPICdata$V1.2,
               NEWdataYZCVnetindBIGPICdata$V1.3,
               NEWdataYZCVnetindBIGPICdata$V2.1,
               NEWdataYZCVnetindBIGPICdata$V2.2,
               NEWdataYZCVnetindBIGPICdata$V3,
               NEWdataYZCVnetindBIGPICdata$V4)

Cmat = cbind(NEWdataYZCVnetindBIGPICdata$C1,
               NEWdataYZCVnetindBIGPICdata$C2,
               NEWdataYZCVnetindBIGPICdata$C3,
               NEWdataYZCVnetindBIGPICdata$C4,
               NEWdataYZCVnetindBIGPICdata$C5,
               NEWdataYZCVnetindBIGPICdata$C6)

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Average Effects
{
  resultBBmediationPOST = BBmediationPOST(BARTfitY, BARTfitM, Cmat, Vmat, Cjmat, Vjmat, Uindex)
  print("Average Effects - BB")
  print(resultBBmediationPOST)
  # allinfo = c(allinfo,c(resultBBmediationPOST$NIE_result_mc,
  #                       resultBBmediationPOST$NDE_result_mc,
  #                       resultBBmediationPOST$ATE_result_mc))
  
  # BB-HBB
  resultHBBmediationPOST = HBBmediationPOST(BARTfitY, BARTfitM, Cmat, Vmat, Cjmat, Vjmat, Uindex)
  print("Average Effects - HBB")
  print(resultHBBmediationPOST)
  # allinfo = c(allinfo,c(resultHBBmediationPOST$NIE_result_mc,
  #                       resultHBBmediationPOST$NDE_result_mc,
  #                       resultHBBmediationPOST$ATE_result_mc))
  
  # TSBB
  resultTSBBmediationPOST = TSBBmediationPOST(BARTfitY, BARTfitM, Cmat, Vmat, Cjmat, Vjmat, Uindex)
  print("Average Effects - TSBB")
  print(resultTSBBmediationPOST)
  # allinfo = c(allinfo,c(resultTSBBmediationPOST$NIE_result_mc,
  #                       resultTSBBmediationPOST$NDE_result_mc,
  #                       resultTSBBmediationPOST$ATE_result_mc))
}

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Conditional Effects
{
  # Conditional Effects - cluster
  # GISE (groups intergrateed saving for empowerment)
  # Prerandomization variable that tell the near the level of GISE penetrtaion around the community.
  # Light & Heavy
  # V3 = ifelse(BIGPICdata$Gise_Penetration=="Heavy",1,0)
  S_V3_0 = (NEWdataYZCVnetindBIGPICdata$V3 == 0) # Light GISE
  S_V3_1 = (NEWdataYZCVnetindBIGPICdata$V3 == 1) # HeavyGISE
  
  TSBBgivenGISElight = VcondTSBBmediationPOST(BARTfitY, BARTfitM, Cmat, Vmat, Cjmat, Vjmat, S_V3_0, Uindex)
  print("Conditional Effects - cluster (GISE (groups intergrateed saving for empowerment) == 0)")
  print(TSBBgivenGISElight)
  # TSBBgivenGISElight$NIE_result_mc
  # TSBBgivenGISElight$NDE_result_mc
  # TSBBgivenGISElight$ATE_result_mc
  # allinfoCOND = c(allinfoCOND,c(TSBBgivenGISElight$NIE_result_mc,
  #                               TSBBgivenGISElight$NDE_result_mc,
  #                               TSBBgivenGISElight$ATE_result_mc))
  
  TSBBgivenGISEheavy = VcondTSBBmediationPOST(BARTfitY, BARTfitM, Cmat, Vmat, Cjmat, Vjmat, S_V3_1, Uindex)
  print("Conditional Effects - cluster (GISE (groups intergrateed saving for empowerment) == 1)")
  print(TSBBgivenGISEheavy)
  # TSBBgivenGISEheavy$NIE_result_mc
  # TSBBgivenGISEheavy$NDE_result_mc
  # TSBBgivenGISEheavy$ATE_result_mc
  # allinfoCOND = c(allinfoCOND,c(TSBBgivenGISEheavy$NIE_result_mc,
  #                               TSBBgivenGISEheavy$NDE_result_mc,
  #                               TSBBgivenGISEheavy$ATE_result_mc))
  
  # Conditional Effects - individual
  # binary: Derived variable - baseline blood sugars
  # "Yes": (dat$bl_bsg_type=="Random"  & dat$bl_bsg>=11.1) or
  # (dat$bl_bsg_type=="Fasting"  & dat$bl_bsg >= 7) or 
  # (dat$bl_ever_diagnosis_1=="Yes") or
  # (dat$bl_ever_diagnosis_2=="Yes")
  # otherwise "NO"
  # C5 = ifelse(BIGPICdata$dm=="Yes",1,0)
  S_C5_0 = (NEWdataYZCVnetindBIGPICdata$C5 == 0)
  S_C5_1 = (NEWdataYZCVnetindBIGPICdata$C5 == 1)
  
  TSBBgivenDM0 = CcondTSBBmediationPOST(BARTfitY, BARTfitM, Cmat, Vmat, Cjmat, Vjmat, S_C5_0, Uindex)
  print("Conditional Effects - individual (baseline blood sugars == 0)")
  print(TSBBgivenDM0)
  # allinfoCOND = c(allinfoCOND,c(TSBBgivenDM0$NIE_result_mc,
  #                               TSBBgivenDM0$NDE_result_mc,
  #                               TSBBgivenDM0$ATE_result_mc))
  
  TSBBgivenDM1 = CcondTSBBmediationPOST(BARTfitY, BARTfitM, Cmat, Vmat, Cjmat, Vjmat, S_C5_1, Uindex)
  print("Conditional Effects - individual (baseline blood sugars == 1)")
  print(TSBBgivenDM1)
  # allinfoCOND = c(allinfoCOND,c(TSBBgivenDM1$NIE_result_mc,
  #                               TSBBgivenDM1$NDE_result_mc,
  #                               TSBBgivenDM1$ATE_result_mc))
  # 
  # allinfoCOND = data.frame(t(allinfoCOND))
  # write.table(allinfoCOND, file = txt.titleCOND, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
}




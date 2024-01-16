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
library(mgcv)
library(mediation)
library(openxlsx) 
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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
summariessss = function(x){
  out = list()
  out$table = table(x)
  out$num = length(x)
  out$na.num = sum(is.na(x))
  
  return(out)
}

# Set Directory
setwd("C:/Users/WooJung/Documents/Rproject/BBmediation/RealDataAnalysis/BIGPIC/BIGPICdata")

# Load Data
load("BIGPICdata.Rdata")

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# LARK data?
# ========== Yij (outcome) 
# => SBP at 12 months (measured at baseline, month 3, and month 12)
#    using a multivariate normal to help address the missingness.
# ========== Mij (mediator)
# => an indicator of linkage to care or a continuous measure of engagement in care
# ========== Cij (individual-level covariates)
# => age, baseline BP, gender, pathway recruitment, ...
# ========== Vj (cluster-level covariates)
# => usage of community savings plans and type of health facility
# ========== Zj (cluster-level treatments)
# => Each cluster is randomized to treatment

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# Clusters: j=1,...,J=24
V0 = BIGPICdata$site
summariessss(V0)
sitetable = table(V0);sitetable
sitenames = names(sitetable)
J = length(sitenames)

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# treatment: 
# UC: Usual Care
# UCMF: Usual Care + Microfinance
# GMV: Group Medical Visit
# GMVMF: Group Medical Visit + Microfinance
summariessss(BIGPICdata$arm.refUC2)
Z = BIGPICdata$arm.refUC2
ZZ = numeric(length(Z))
ZZ[which(Z=="UC")] = 0
ZZ[which(Z=="UCMF")] = 1
ZZ[which(Z=="GMV")] = 2
ZZ[which(Z=="GMVMF")] = 3
Z = ZZ

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# outcome: the 1-year change in systolic blood pressure (SBP).
# BIGPICdata$bl_sbp_ave
# BIGPICdata$m3_sbp_ave
# BIGPICdata$m12_sbp_ave
# BIGPICdata$sbp.change.3mos       
# BIGPICdata$sbp.change.12mos      
summariessss(BIGPICdata$bl_sbp_ave)
summariessss(BIGPICdata$m3_sbp_ave)
summariessss(BIGPICdata$m12_sbp_ave)
summariessss(BIGPICdata$sbp.change.3mos)
summariessss(BIGPICdata$sbp.change.12mos)
Y = BIGPICdata$sbp.change.12mos

# BIGPICdata$bl_dbp_ave
# BIGPICdata$m3_dbp_ave
# BIGPICdata$m12_dbp_ave
# BIGPICdata$dbp.change.3mos
# BIGPICdata$dbp.change.12mos

# cbind(sitenames,as.character(sapply(1:J, function(l) unlist(
#   BIGPICdata$arm.refUC2)[which(sitenames[l] == BIGPICdata$site)[1]])))

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# potential cluster-level confounder
# BIGPIC County: 
# "Busia", "Kisumu", "Uasin Gishu", "Trans-Nzoia"
factorV1 = factor(BIGPICdata$county)
summariessss(factorV1)
Kv1 = length(levels(factorV1))
V1 = sapply(2:Kv1, function(l) ifelse(factorV1==levels(factorV1)[l],1,0))

# HF: Prerandomization variable that tell the type/leve of health facility near the community
# Dis: Dispensary
# HC: health center 
# SCH: Sub county hospital
factorV2 = factor(BIGPICdata$HF)
summariessss(factorV2)
Kv2 = length(levels(factorV2))
V2 = sapply(2:Kv2, function(l) ifelse(factorV2==levels(factorV2)[l],1,0))

# GISE (groups intergrateed saving for empowerment)
# Prerandomization variable that tell the near the leve of GISE penetrtaion around the community. 
# Light & Heavy
V3 = ifelse(BIGPICdata$Gise_Penetration=="Heavy",1,0)

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# potential individual-level confounder
# age
C1 = BIGPICdata$age_est
summariessss(C1)

# gender
C2 = ifelse(BIGPICdata$gender=="Male",1,0)
summariessss(C2)

# International Wealth Index
summariessss(BIGPICdata$iwi.change)
summariessss(BIGPICdata$bl_iwi)
summariessss(BIGPICdata$m12_iwi)

# Derived variable- Variables used for QRISK calculation are; sbp, age, family history of angina, 
# smoking, cholesterol ratio, Atrial fibrillation, type II diabetes, type I diabetes, hypertension 
# treatment/medication, Chronic renal disease and rheumatoid arthritis. Note that QRISK score is 
# not measured at 3 months.
# QRISK3
summariessss(BIGPICdata$qrisk3)
summariessss(BIGPICdata$bl_QRISK3)
summariessss(BIGPICdata$m12_QRISK3)

# Do you have a father, mother, brother, sister, or own child with type 1 or type 2 diabetes?
summariessss(BIGPICdata$fhx_diabetes)

# QF-1: In the past 3 months, have you TAKEN any medicines that were prescribed by your health
# care provider (excluding herbalists or spiritual healer), for either hypertension, diabetes,
# stroke, or heart disease?
summariessss(BIGPICdata$pres_med)

# BMI
summariessss(BIGPICdata$waist)
summariessss(BIGPICdata$height)
summariessss(BIGPICdata$weight)
summariessss(BIGPICdata$bmi)

# Livestock (Cows, Sheep, Goats, Donkeys, etc) : None or Some 
C5 = ifelse(BIGPICdata$livestock.cat=="Some",1,0)
summariessss(C5)

# binary: Unemployed or Employed
C6 = ifelse(BIGPICdata$earnings_status=="Employed",1,0)
summariessss(C6)

# ordial: Not formally employed, less that 
# Ksh 1000, Ksh 1000-2999, Ksh 3000-4999, Ksh >= 5000, Don't Know, Refused
summariessss(BIGPICdata$monthly_earn)

# binary: elevated BP (SBP ≥140 or diastolic BP (DBP) ≥ 90)
C7 = ifelse(BIGPICdata$elavated_bp=="Sbp/Dbp >= 140/90",1,0)
summariessss(C7)

# Derived variable - Leicester score
# The Leicester Risk Assessment score (Appendix B) was adapted from the prospectively
# validated FINDRISC score and has been prospectively validated in multi-ethnic 
# populations to predict development of diabetes
summariessss(BIGPICdata$Lscore)

# binary: Derived variable - baseline blood sugars
# "Yes": (dat$bl_bsg_type=="Random"  & dat$bl_bsg>=11.1) or
# (dat$bl_bsg_type=="Fasting"  & dat$bl_bsg >= 7) or 
# (dat$bl_ever_diagnosis_1=="Yes") or
# (dat$bl_ever_diagnosis_2=="Yes")
# otherwise "NO"
C8 = ifelse(BIGPICdata$dm=="Yes",1,0)
summariessss(C8)

V = cbind(V1, V2, V3)
C = cbind(C1,C2,C5,C6,C7,C8)

dataYZCV = data.frame(ID = BIGPICdata$participant_id, Y=Y, Z=Z, C=C, V=V, V0 = as.numeric(V0))

# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# -----------------------------------------------------------------------
# mediator (intermediate factors): 
# (h) physical activity data => {the number of days/week} or {minutes/week}
setwd("C:/Users/WooJung/Documents/Rproject/BBmediation/RealDataAnalysis/BIGPIC/BIGPICdata/DataExcel")
a=read.csv("cost_if_singlevisits.csv")
aa = a$participant_id
unique_aa = unique(aa)
M = a$moderate_days[a$event_name == "3_months_follow-up"]

dataM = data.frame(ID=a$participant_id[a$event_name == "3_months_follow-up"], M=M)
tempdata = merge(dataYZCV, dataM, by = c("ID"))
nrow(tempdata)

tempdata[which(apply(ifelse(is.na(tempdata),0,1),1,prod)==0),]
apply(ifelse(is.na(tempdata),0,1),2,prod)

tempdata = tempdata[-which(apply(ifelse(is.na(tempdata),0,1),1,prod)==0),]
nrow(tempdata)

# # (j) a lot of variables social network characteristics; mostly ordinal or categorical level larger than 2.
# setwd("C:/Users/WooJung/Documents/Rproject/BBmediation/RealDataAnalysis/BIGPIC/BIGPICdata/DataExcel")
# a=read.csv("sns.csv")
# 
# aa = a$participant_id
# unique_aa = unique(aa)
# aaa = a$ppl_talkto
# M = a$ppl_talkto[a$event_name == "3_months_follow-up"]
# summariessss(M)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# esttype = "mean"
esttype = "median"
save_cluster = TRUE
# save_cluster = FALSE

# Extract ID for simulated dataset (specific to LSF computing cluster)
# Note: The LSB_JOBINDEX is specified in the bsub command using the -J
# option
run_ID = 1

# Specify the seed so can repeat simulation later if necessary
set.seed(run_ID)

# Y = tempdata$Y
# M = tempdata$M
# Z = tempdata$Z
# C = cbind(tempdata$C.C1,tempdata$C.C2,tempdata$C.C5,tempdata$C.C6,tempdata$C.C7,tempdata$C.C8)
# V = cbind(tempdata$V.V1,tempdata$V.V2,tempdata$V.V3,tempdata$V.V4,tempdata$V.V5,tempdata$V.V3.1)
# V0 = tempdata$V0
tempdata01 = tempdata[-which(tempdata$Z==2),]
tempdata01 = tempdata01[-which(tempdata01$Z==3),]

Y_01 = tempdata01$Y
M_01 = tempdata01$M
Z_01 = tempdata01$Z
C_01 = cbind(tempdata01$C.C1,tempdata01$C.C2,tempdata01$C.C5,tempdata01$C.C6,tempdata01$C.C7,tempdata01$C.C8)
V_01 = cbind(tempdata01$V.V1,tempdata01$V.V2,tempdata01$V.V3,tempdata01$V.V4,tempdata01$V.V5,tempdata01$V.V3.1)
V0_01 = tempdata01$V0

M_01 = M_01 + 1
BBmediation_result = BBmediation(Y_01, M_01, Z_01, C_01, V_01, V0_01,
                                 gibbs_iter, gibbs_burnin, gibbs_mc,
                                 esttype, save_cluster)

BBmediation_result$NIE_result_mc
BBmediation_result$NDE_result_mc
BBmediation_result$ATE_result_mc

# mean(tempdata$M[tempdata$Z==0])
# mean(tempdata$M[tempdata$Z==1])
# mean(tempdata$M[tempdata$Z==2])
# mean(tempdata$M[tempdata$Z==3])

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# # Set Directory
# setwd("C:/Users/WooJung/Documents/Rproject/BBmediation/RealDataAnalysis/BIGPICdata")
# 
# # Load Data
# load("BIGPICdata.Rdata")
# names(BIGPICdata)
# BIGPICdata$participant_id        
# BIGPICdata$encounter_date        
# BIGPICdata$county                
# BIGPICdata$site                  
# BIGPICdata$study_arm             
# BIGPICdata$arm.refUC2            
# BIGPICdata$pathway               
# BIGPICdata$HF                    
# BIGPICdata$Gise_Penetration      
# BIGPICdata$age_est               
# BIGPICdata$gender                
# BIGPICdata$age.c                 
# BIGPICdata$fhx_diabetes          
# BIGPICdata$pres_med              
# BIGPICdata$waist                 
# BIGPICdata$height                
# BIGPICdata$bmi                   
# BIGPICdata$weight                
# BIGPICdata$Lscore                
# BIGPICdata$Lscore.cat            
# BIGPICdata$dm                    
# BIGPICdata$elavated_bp           
# BIGPICdata$qrisk3                
# BIGPICdata$iwi.cat               
# BIGPICdata$iwi.change            
# BIGPICdata$land_cat              
# BIGPICdata$land_changes          
# BIGPICdata$livestock.cat         
# BIGPICdata$earnings_status       
# BIGPICdata$monthly_earn          
# BIGPICdata$sbp.change.3mos       
# BIGPICdata$dbp.change.3mos       
# BIGPICdata$sbp.change.12mos      
# BIGPICdata$dbp.change.12mos      
# BIGPICdata$GroupsScheduled       
# BIGPICdata$GroupsAttended        
# BIGPICdata$PercentAttended       
# BIGPICdata$PropAttended          
# BIGPICdata$first_group_date      
# BIGPICdata$bigpic_group_id       
# BIGPICdata$disposition_3m        
# BIGPICdata$disposition_12m       
# BIGPICdata$fup_status            
# BIGPICdata$bl_clinicform_date    
# BIGPICdata$bl_sbp2               
# BIGPICdata$bl_dbp2               
# BIGPICdata$bl_sbp3               
# BIGPICdata$bl_dbp3               
# BIGPICdata$bl_sbp_ave            
# BIGPICdata$bl_dbp_ave            
# BIGPICdata$bl_sbp_ave_high       
# BIGPICdata$bl_dbp_ave_high       
# BIGPICdata$bl_bsg                
# BIGPICdata$bl_bsg_type           
# BIGPICdata$bl_total_chol         
# BIGPICdata$bl_hdl_chol           
# BIGPICdata$bl_ldl_chol           
# BIGPICdata$bl_trigly             
# BIGPICdata$bl_total_chol.cat     
# BIGPICdata$bl_ldl_chol.cat       
# BIGPICdata$bl_hdl_chol.cat       
# BIGPICdata$bl_trigly.cat         
# BIGPICdata$bl_bp_control         
# BIGPICdata$bl_encounter_date_sns 
# BIGPICdata$bl_primary_job        
# BIGPICdata$bl_other_prijobs      
# BIGPICdata$bl_hours_worked       
# BIGPICdata$bl_earn               
# BIGPICdata$bl_reason_nowork      
# BIGPICdata$bl_land_acres         
# BIGPICdata$bl_land               
# BIGPICdata$bl_livestock          
# BIGPICdata$bl_cows               
# BIGPICdata$bl_sheep              
# BIGPICdata$bl_goats              
# BIGPICdata$bl_donkey             
# BIGPICdata$bl_QRISK3             
# BIGPICdata$bl_iwi                
# BIGPICdata$bl_highest_level      
# BIGPICdata$bl_primary_level      
# BIGPICdata$bl_secondary_level    
# BIGPICdata$bl_university_level   
# BIGPICdata$bl_postgrad_level     
# BIGPICdata$bl_nhif               
# BIGPICdata$bl_none_everdiag      
# BIGPICdata$bl_ever_diagnosis_1   
# BIGPICdata$bl_ever_diagnosis_2   
# BIGPICdata$bl_earnings           
# BIGPICdata$bl_land_acres.cat     
# BIGPICdata$m3_clinicform_date    
# BIGPICdata$m3_sbp_ave            
# BIGPICdata$m3_dbp_ave            
# BIGPICdata$m3_bp_control         
# BIGPICdata$m12_clinicform_date   
# BIGPICdata$m12_sbp2              
# BIGPICdata$m12_dbp2              
# BIGPICdata$m12_sbp3              
# BIGPICdata$m12_dbp3              
# BIGPICdata$m12_sbp_ave           
# BIGPICdata$m12_dbp_ave           
# BIGPICdata$m12_sbp_ave_high      
# BIGPICdata$m12_dbp_ave_high      
# BIGPICdata$m12_bsg               
# BIGPICdata$m12_bsg_type          
# BIGPICdata$m12_total_chol        
# BIGPICdata$m12_hdl_chol          
# BIGPICdata$m12_ldl_chol          
# BIGPICdata$m12_trigly            
# BIGPICdata$m12_encounter_date_sns
# BIGPICdata$m12_primary_job       
# BIGPICdata$m12_other_prijobs     
# BIGPICdata$m12_hours_worked      
# BIGPICdata$m12_earn              
# BIGPICdata$m12_reason_nowork     
# BIGPICdata$m12_land_acres        
# BIGPICdata$m12_land              
# BIGPICdata$m12_livestock         
# BIGPICdata$m12_cows              
# BIGPICdata$m12_sheep             
# BIGPICdata$m12_goats             
# BIGPICdata$m12_donkey            
# BIGPICdata$m12_QRISK3            
# BIGPICdata$m12_iwi               
# BIGPICdata$m12_highest_level     
# BIGPICdata$m12_primary_level     
# BIGPICdata$m12_secondary_level   
# BIGPICdata$m12_university_level  
# BIGPICdata$m12_postgrad_level    
# BIGPICdata$m12_nhif              
# BIGPICdata$m12_none_everdiag     
# BIGPICdata$m12_ever_diagnosis_1  
# BIGPICdata$m12_ever_diagnosis_2  
# BIGPICdata$m12_bp_control

# # -----------------------------------------------------------------------------
# # -----------------------------------------------------------------------------
# # -----------------------------------------------------------------------------
# # names(larkdata)
# # N = nrow(larkdata)
# # ncol(larkdata)
# # apply(is.na(larkdata), 2, sum)
# # head(larkdata)
# 
# # larkdata$pstudy_id = as.numeric(larkdata$pstudy_id)
# # larkdata$encounter_date.d
# # larkdata$followup_date.d
# # larkdata$dob
# # larkdata$age
# # Male = 1 / Female = 0
# larkdata$gender = ifelse(larkdata$gender == "Male", 1, 0)
# # larkdata$comm_unit
# # Mosoriot = 1 / Turbo = 0
# larkdata$Division = ifelse(larkdata$Division == "Mosoriot", 1, 0)
# # UsualCare = 0 / PaperBased = 1 / SmartPhone = 2
# larkdata$Arm[which(larkdata$Arm=="UsualCare")] = 0
# larkdata$Arm[which(larkdata$Arm=="PaperBased")] = 1
# larkdata$Arm[which(larkdata$Arm=="SmartPhone")] = 2
# larkdata$Arm = as.numeric(larkdata$Arm)
# # aBaraza = 1 / bNewly identified PHCT = 2 / cPreviously identified not linked = 3
# # dPreviously identified not retained = 4 / eHealth Facility = 5 / fIdentified by CHW = 6
# larkdata$enrollment_strategy[which(larkdata$enrollment_strategy=="aBaraza")] = 1
# larkdata$enrollment_strategy[which(larkdata$enrollment_strategy=="bNewly identified PHCT")] = 2
# larkdata$enrollment_strategy[which(larkdata$enrollment_strategy=="cPreviously identified not linked")] = 3
# larkdata$enrollment_strategy[which(larkdata$enrollment_strategy=="dPreviously identified not retained")] = 4
# larkdata$enrollment_strategy[which(larkdata$enrollment_strategy=="eHealth Facility")] = 5
# larkdata$enrollment_strategy[which(larkdata$enrollment_strategy=="fIdentified by CHW")] = 6
# larkdata$enrollment_strategy = as.numeric(larkdata$enrollment_strategy)
# # Follow-up entered = 1 / Relocated = 2 / Denied = 3 / Deceased = 4
# # Out of Catchment = 5 / Unreachable = 6 / Entry Missing = 7 / Unknown
# larkdata$Disposition[which(larkdata$Disposition=="Follow-up entered")] = 1
# larkdata$Disposition[which(larkdata$Disposition=="Relocated")] = 2
# larkdata$Disposition[which(larkdata$Disposition=="Denied")] = 3
# larkdata$Disposition[which(larkdata$Disposition=="Deceased")] = 4
# larkdata$Disposition[which(larkdata$Disposition=="Out of Catchment")] = 5
# larkdata$Disposition[which(larkdata$Disposition=="Unreachable")] = 6
# larkdata$Disposition[which(larkdata$Disposition=="Entry Missing")] = 7
# larkdata$Disposition[which(larkdata$Disposition=="Unknown")] = "Unknown"
# # larkdata$systolic_bp
# # larkdata$fu_systolic_bp
# # larkdata$diastolic_bp
# # larkdata$fu_diastolic_bp
# # No = 0 / Yes = 1 / Unknown
# # Link to the chronic disease management program (CDM) per 
# # AMPATH Medical Record System (AMRS) within a year of enrollment
# larkdata$Link_CDM_1Year.unk[which(larkdata$Link_CDM_1Year.unk=="No")] = 0
# larkdata$Link_CDM_1Year.unk[which(larkdata$Link_CDM_1Year.unk=="Yes")] = 1
# larkdata$Link_CDM_1Year.unk[which(larkdata$Link_CDM_1Year.unk=="Unknown")] = "Unknown"
# # No = 0 / Yes = 1 / Unknown
# # Link to CDM per AMRS at any time after enrollment
# larkdata$Link_CDM.unk[which(larkdata$Link_CDM.unk=="No")] = 0
# larkdata$Link_CDM.unk[which(larkdata$Link_CDM.unk=="Yes")] = 1
# larkdata$Link_CDM.unk[which(larkdata$Link_CDM.unk=="Unknown")] = "Unknown"
# # No = 0 / Yes = 1 / Unknown
# # Link to care based on checking "Yes" on tracker
# larkdata$Link_Track.unk[which(larkdata$Link_Track.unk=="No")] = 0
# larkdata$Link_Track.unk[which(larkdata$Link_Track.unk=="Yes")] = 1
# larkdata$Link_Track.unk[which(larkdata$Link_Track.unk=="Unknown")] = "Unknown"
# # No = 0 / Yes = 1 / Unknown
# # Link per tracker based on any information in the tracker form
# larkdata$Link_All.unk[which(larkdata$Link_All.unk=="No")] = 0
# larkdata$Link_All.unk[which(larkdata$Link_All.unk=="Yes")] = 1
# larkdata$Link_All.unk[which(larkdata$Link_All.unk=="Unknown")] = "Unknown"
# # No = 0 / Yes = 1 / z.Missing
# larkdata$mod.equiv.gt150.txt[which(larkdata$mod.equiv.gt150.txt=="No")] = 0
# larkdata$mod.equiv.gt150.txt[which(larkdata$mod.equiv.gt150.txt=="Yes")] = 1
# larkdata$mod.equiv.gt150.txt[which(larkdata$mod.equiv.gt150.txt=="z.Missing")] = "Missing"
# # aNo_Job = 1 / bFarmer = 2 / cBusiness_person = 3 / dPublic_sector = 4
# # eStudent = 5 / fOther = 6 / z.Missing = (-1e5)
# larkdata$job[which(larkdata$job=="aNo_Job")] = 1
# larkdata$job[which(larkdata$job=="bFarmer")] = 2
# larkdata$job[which(larkdata$job=="cBusiness_person")] = 3
# larkdata$job[which(larkdata$job=="dPublic_sector")] = 4
# larkdata$job[which(larkdata$job=="eStudent")] = 5
# larkdata$job[which(larkdata$job=="fOther")] = 6
# larkdata$job[which(larkdata$job=="z.Missing")] = "Missing"
# # a.nojob = 1 / b.lt5000 = 2 / c.ge5000.lt10000 = 3 / d.ge10000.lt20000 = 4
# # e.ge20000.lt30000 = 5 / f.ge30000 = 6 / z.Missing = (-1e5)
# larkdata$ord.earn[which(larkdata$ord.earn=="a.nojob")] = 1
# larkdata$ord.earn[which(larkdata$ord.earn=="b.lt5000")] = 2
# larkdata$ord.earn[which(larkdata$ord.earn=="c.ge5000.lt10000")] = 3
# larkdata$ord.earn[which(larkdata$ord.earn=="d.ge10000.lt20000")] = 4
# larkdata$ord.earn[which(larkdata$ord.earn=="e.ge20000.lt30000")] = 5
# larkdata$ord.earn[which(larkdata$ord.earn=="f.ge30000")] = 6
# larkdata$ord.earn[which(larkdata$ord.earn=="z.Missing")] = "Missing"
# # No = 0 / Yes = 1
# larkdata$Employed = ifelse(larkdata$Employed == "Yes", 1, 0)
# # No = 0 / Yes = 1 / Missing = (-1e5)
# larkdata$high_bp_12mos[which(larkdata$high_bp_12mos=="No")] = 0
# larkdata$high_bp_12mos[which(larkdata$high_bp_12mos=="Yes")] = 1
# larkdata$high_bp_12mos[which(larkdata$high_bp_12mos=="z.Missing")] = "Missing"
# # No = 0 / Yes = 1 / z.Missing = (-1e5)
# larkdata$nhif.yn[which(larkdata$nhif.yn=="No")] = 0
# larkdata$nhif.yn[which(larkdata$nhif.yn=="Yes")] = 1
# larkdata$nhif.yn[which(larkdata$nhif.yn=="z.Missing")] = "Missing"
# # No = 0 / Yes = 1 / z.Missing = (-1e5)
# larkdata$tobacco.txt[which(larkdata$tobacco.txt=="No")] = 0
# larkdata$tobacco.txt[which(larkdata$tobacco.txt=="Yes")] = 1
# larkdata$tobacco.txt[which(larkdata$tobacco.txt=="z.Missing")] = "Missing"
# # Chew = 1 / Cigarettes = 2 / Pipe = 3 / Sniffed = 4 / NA
# larkdata$tobacco.type.txt[which(larkdata$tobacco.type.txt=="Chew")] = 1
# larkdata$tobacco.type.txt[which(larkdata$tobacco.type.txt=="Cigarettes")] = 2
# larkdata$tobacco.type.txt[which(larkdata$tobacco.type.txt=="Pipe")] = 3
# larkdata$tobacco.type.txt[which(larkdata$tobacco.type.txt=="Sniffed")] = 4
# # No = 0 / Yes = 1 / z.Missing = (-1e5)
# larkdata$alcohol.txt[which(larkdata$alcohol.txt=="No")] = 0
# larkdata$alcohol.txt[which(larkdata$alcohol.txt=="Yes")] = 1
# larkdata$alcohol.txt[which(larkdata$alcohol.txt=="z.Missing")] = "Missing"
# # a<5 = 1 / b6-10 = 2 / c11-15 = 3 / z.Missing = (-1e5) / NA
# larkdata$alcohol.servingwk.txt[which(larkdata$alcohol.servingwk.txt=="a<5")] = 1
# larkdata$alcohol.servingwk.txt[which(larkdata$alcohol.servingwk.txt=="b6-10")] = 2
# larkdata$alcohol.servingwk.txt[which(larkdata$alcohol.servingwk.txt=="c11-15")] = 3
# larkdata$alcohol.servingwk.txt[which(larkdata$alcohol.servingwk.txt=="z.Missing")] = "Missing"


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# 1. Compute the NIE, NDE and ATE for particular combinations of the 
# cluster level covariates (the cluster level causal mediation effects) 
# - maybe just based on GISE and type of facility combo (‘integrating’ over county)







# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# 2. Compute the causal effects for a few combinations of the individual 
# level covariates (individual level causal mediation effects) - maybe just 
# pick 3 or 4 different combinations of the individual level covariates.  





















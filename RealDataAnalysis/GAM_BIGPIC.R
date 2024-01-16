# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Extract ID for simulated dataset (specific to LSF computing cluster)
# Note: The LSB_JOBINDEX is specified in the bsub command using the -J
# option
run_ID = 1

# Specify the seed so can repeat simulation later if necessary
set.seed(run_ID)

# Y = Y_23
# M = M_23
# Z = Z_23
# C = C_23
# V = V_23
# V0 = V0_23

# # This code must be in your current directory or you can change the path.
# setwd("C:/Users/WooJung/Documents/Rproject/BBmediation/RealDataAnalysis/BIGPICattendace(continuous)")
# 
# # Load R code
# source("BIGPIC_BBsource_r.R")

gibbs_iter = 2e3;gibbs_burnin = 2e3;gibbs_thin = 10;esttype = "median";save_cluster = TRUE
MCMCresult = BBmediationMCMCconcon(Y_23, M_23, Z_23, C_23, V_23, V0_23,
                                   gibbs_iter, gibbs_burnin, gibbs_thin)
POSTresult = BBmediationPOSTconcon(MCMCresult, esttype, save_cluster)

POSTresult$NIE_result_mc
POSTresult$NDE_result_mc
POSTresult$ATE_result_mc

# 
a0 = unique(V0[Z==0]);aa0 = length(a0)
a1 = unique(V0[Z==1]);aa1 = length(a1)
mean(sapply(1:aa0, function(l) mean(Y[which(V0==a0[l])],na.rm=TRUE)))
mean(sapply(1:aa1, function(l) mean(Y[which(V0==a1[l])],na.rm=TRUE)))
mean(Y[Z==0],na.rm=T)
mean(Y[Z==1],na.rm=T)

# -----------------------------------------------------------------------
med.gam.reg = gam(M ~ Z + C.C1 + C.C2 + C.C3 + C.C4 + C.C5 + C.C6 + V.V1 + V.V2 + V.V3 + V.V4 + V.V5 + V.V3.1, data = tempdata23)
out.gam.reg = gam(Y ~ M + Z + C.C1 + C.C2 + C.C3 + C.C4 + C.C5 + C.C6 + V.V1 + V.V2 + V.V3 + V.V4 + V.V5 + V.V3.1, data = tempdata23)
mediation.result = mediate(med.gam.reg, out.gam.reg, treat = "Z", mediator = "M", boot = TRUE)

# Causal Effects = c("NIE","NDE","TE")
CE = c(mediation.result$d0,mediation.result$z0,mediation.result$tau.coef);CE
c(mediation.result$d0,mediation.result$d0.ci,as.numeric(diff(mediation.result$d0.ci))) # indirect effect
c(mediation.result$z0,mediation.result$z0.ci,as.numeric(diff(mediation.result$z0.ci))) # direct effect
c(mediation.result$tau.coef,mediation.result$tau.ci,as.numeric(diff(mediation.result$tau.ci))) # total effect

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




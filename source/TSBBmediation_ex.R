# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Load R packages
library(MASS)
library(Rcpp)
library(RcppArmadillo)

library(TSBBmediation)
# library(dbarts)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
# Sys.setenv(USE_CXX14 = 1)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# This code must be in your current directory or you can change the path.
setwd("/Users/WooJung/Documents/Rproject/TSBBmediation/source")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Load R code
source("TSBBmediation_source_r.R")

# Load cpp code
sourceCpp("TSBBmediation_source_cpp.cpp")

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Scenario
Scenario = 3
run_ID = 1
set.seed(run_ID)

# Distribution of M
# Mdist = "continuous"
Mdist = "binary"
# Mdist = "ordinal" # ; K = 7 # the number of ordinal mediator when Mdist == "Ordinal"
# Mdist = "count"

# Distribution of Y
Ydist = "continuous"
# Ydist = "binary"

# simple = TRUE
simple = FALSE

{
  esttype = "mean"
  # esttype = "median"
  saveall = TRUE # FALSE
  
  # Extract ID for simulated dataset (specific to LSF computing cluster)
  # Note: The LSB_JOBINDEX is specified in the bsub command using the -J
  # option
  # run_ID = 1
  
  # Specify the seed so can repeat simulation later if necessary
  # set.seed(run_ID)
  
  # txt.title = paste0("Results/D",Scenario,"_result.txt")
  # if (run_ID == 1) {
  #   df = data.frame(matrix(ncol = 16, nrow = 0))
  #   df_col_names = c("run_ID",
  #                    "simEY00", "simEY01", "simEY11",
  #                    "EY0", "EY01", "EY1",
  #                    "NIE_table", "NDE_table", "E_Yz1m1_mc", 
  #                    "E_Mz0_mc", "E_Mz1_mc",
  #                    "EM0", "EM1", "simEM00", "simEM11")
  #   colnames(df)=df_col_names
  #   write.table(df, file = txt.title, sep="\t", row.names = FALSE, col.names = TRUE)
  # }
  
  # Define of constants (adjust to fit the data generating scenario)---------
  
  # Define number of clusters
  # Define number of observations for each dataset
  # Define mean of mediators (count)
  J = 12
  n_j = 50
  # n_j = 100
  # n_j = 200
  N = J * n_j
  
  # ------------------------------------------------------------------------------
  # Define of constants (adjust to fit the data generating scenario) -------------
  # ------------------------------------------------------------------------------
  
  # the number of MC integration per iteration for a MCMC chain
  gibbs_thin = 1e0
  
  # Define number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  # Define number of MCMC draws to 'burn' in Gibbs Sampler (check convergence)
  # gibbs_iter = 1e3 * gibbs_thin
  # gibbs_burnin = 5e4
  if (N > 1500) {
    # for big N (e.g. N = 1500), can use less samples (check posterior)
    gibbs_iter = 1e3
    gibbs_burnin = 1e3
  } else if (N > 1000) {
    # for mid N (e.g. N = 1000), need more samples
    gibbs_iter = 1e3
    gibbs_burnin = 2e3
  } else {
    # for small N (e.g. N = 500), need more samples
    gibbs_iter = 2e3
    gibbs_burnin = 2e3
  }
  gibbs_GLM = 1
  
  # ------------------------------------------------------------------------------
  # Data Generation --------------------------------------------------------------
  # ------------------------------------------------------------------------------
  
  # Data Generation --------------------------------------------------------------
  temp_data = generate_data(J,N,Scenario,Mdist,simple,Ydist)
  
  # Load data for the specific dataset by run_ID
  
  # Load outcome
  Y = temp_data$Y
  
  # Load mediator
  M = temp_data$M
  
  # Load cluster-level treatment
  Z = temp_data$Z
  
  # Load individual-level confounder
  C = temp_data$C
  
  # Load cluster-level confounder
  V = temp_data$V
  
  # Load cluster
  V0 = temp_data$V0
  
  # Causal Effects
  E_true = temp_data$E_true
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Define the number of clusters
  uniqueV0 = sort(unique(V0))
  J = length(uniqueV0)
  V0 = apply(sapply(1:J, function(l) ifelse(V0==uniqueV0[l],l,0)),1,sum)
  
  # Define the number of observations
  N = length(Y)
  
  # Define the number of observations for each cluster
  n_j = as.numeric(table(V0))
  
  # Data reordering
  order_V0 = order(V0)
  Y = Y[order_V0]
  M = M[order_V0]
  Z = Z[order_V0]
  C = C[order_V0,]
  V = V[order_V0,]
  V0 = V0[order_V0]
  
  colnames(C) = paste0("C",1:ncol(C))
  colnames(V) = paste0("V",1:ncol(V))
  
  matX = cbind(Z, C, V)
  matM = cbind(M, Z, C, V)
  
  df = data.frame(Y = Y, M = M, Z = Z, C = C, V = V, V0 = V0, matX = matX, matM = matM)
}

gibbs_iter = 2e3
gibbs_burnin = 2e3

# sparse = TRUE
sparse = FALSE
# ntree = 20
# ntree = 50
# ntree = 100
ntree = 200
# ntree = 2000

c(E_true[1],E_true[2],E_true[3])
c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])

length(unique(BARTfit$uMdraw[,12]))
length(unique(BARTfit$uYdraw[,12]))
length(unique(BARTfit$RHOdraw[,12]))

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# rBARTmediation
BARTfit = rBARTmediation(Y, M, Z, C, V, V0,
                         typeY = Ydist, typeM = Mdist, 
                         ndpost=gibbs_iter, nskip=gibbs_burnin, keepevery=gibbs_thin, 
                         ntree=ntree, sparse = sparse)

# BB-BB
rBARTmediationresultBBmediationPOST = BBmediationPOST(BARTfit, C, V, V0, esttype, saveall, T)
(rBARTmediationBBtableNIE = rBARTmediationresultBBmediationPOST$NIE_result_mc)
(rBARTmediationBBtableNDE = rBARTmediationresultBBmediationPOST$NDE_result_mc)
(rBARTmediationBBtableATE = rBARTmediationresultBBmediationPOST$ATE_result_mc)

plot(rBARTmediationresultBBmediationPOST$NIE)
plot(rBARTmediationresultBBmediationPOST$NDE)
plot(rBARTmediationresultBBmediationPOST$ATE)

# BB-HBB
rBARTmediationresultHBBmediationPOST = HBBmediationPOST(BARTfit, C, V, V0, esttype, saveall, T)
(rBARTmediationHBBtableNIE = rBARTmediationresultHBBmediationPOST$NIE_result_mc)
(rBARTmediationHBBtableNDE = rBARTmediationresultHBBmediationPOST$NDE_result_mc)
(rBARTmediationHBBtableATE = rBARTmediationresultHBBmediationPOST$ATE_result_mc)

# TSBB
# rBARTmediationresultTSBBmediationPOST = TSBBmediationPOST(BARTfit, C, V, V0, esttype, saveall, chi = 1, zeta = 0.5, T)
# (rBARTmediationTSBBtableNIE = rBARTmediationresultTSBBmediationPOST$NIE_result_mc)
# (rBARTmediationTSBBtableNDE = rBARTmediationresultTSBBmediationPOST$NDE_result_mc)
# (rBARTmediationTSBBtableATE = rBARTmediationresultTSBBmediationPOST$ATE_result_mc)

rBARTmediationresultTSBBmediationPOSTsim = TSBBmediationPOSTsim(BARTfit, C, V, V0, esttype, saveall,
                                                                list_chi = c(1e-2, 1e-1, 1, 1e1, 1e2),
                                                                list_zeta = c(0, 0.5, 1), T)
# cbind(rBARTmediationresultTSBBmediationPOSTsim$NIE_result_mc,
#       rBARTmediationresultTSBBmediationPOSTsim$NDE_result_mc,
#       rBARTmedia# c(E_true[1],E_true[2],E_true[3])tionresultTSBBmediationPOSTsim$ATE_result_mc)

round(rbind(c(rBARTmediationresultBBmediationPOST$NIE_result_mc,
              rBARTmediationresultBBmediationPOST$NDE_result_mc,
              rBARTmediationresultBBmediationPOST$ATE_result_mc),
            c(rBARTmediationresultHBBmediationPOST$NIE_result_mc,
              rBARTmediationresultHBBmediationPOST$NDE_result_mc,
              rBARTmediationresultHBBmediationPOST$ATE_result_mc),
            cbind(rBARTmediationresultTSBBmediationPOSTsim$NIE_result_mc,
                  rBARTmediationresultTSBBmediationPOSTsim$NDE_result_mc,
                  rBARTmediationresultTSBBmediationPOSTsim$ATE_result_mc)),3)
c(E_true[1],E_true[2],E_true[3])
c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])

# allinfo = c()
# list_chi = c(1e-2, 1e-1, 1, 1e1, 1e2)
# list_zeta = c(0, 0.5, 1)
# for (chi in list_chi){
#   for (zeta in list_zeta){
#     # ------------------------------------------------------------------------------
#     # Compputation steps -----------------------------------------------------------
#     # ------------------------------------------------------------------------------
#     # TSBB
#     resultTSBBmediationPOST = TSBBmediationPOST(BARTfit , C, V, V0, esttype, saveall, chi, zeta, T)
#     allinfo = c(allinfo,c(resultTSBBmediationPOST$NIE_result_mc,
#                           resultTSBBmediationPOST$NDE_result_mc,
#                           resultTSBBmediationPOST$ATE_result_mc))
#     rm(resultTSBBmediationPOST);gc()
#   }
# }
# allinfo = matrix(allinfo, nrow = 15)
# allinfo = t(cbind(c(rBARTmediationBBtableNIE,
#                     rBARTmediationBBtableNDE,
#                     rBARTmediationBBtableATE),
#                   c(rBARTmediationHBBtableNIE,
#                     rBARTmediationHBBtableNDE,
#                     rBARTmediationHBBtableATE),
#                   allinfo))

# # -----------------------------------------------------------------------------
# # -----------------------------------------------------------------------------
# # dbarts
# Mpost = rbart_vi(M ~ matX, df, group.by = V0,
#                  n.samples = gibbs_iter, n.burn = gibbs_burnin, n.thin = gibbs_thin,
#                  n.chains = 1, n.threads = 1)
# Ypost = rbart_vi(Y ~ matM, df, group.by = V0,
#                  n.samples = gibbs_iter, n.burn = gibbs_burnin, n.thin = gibbs_thin,
#                  n.chains = 1, n.threads = 1)
# rbartBARTfit = list(Yobject = Ypost, Mobject = Mpost)
# class(rbartBARTfit) = class(Ypost)
# 
# # BB-BB
# dbartresultBBmediationPOSTrbart = BBmediationPOST(rbartBARTfit, C, V, V0, esttype, saveall)
# (dbartBBtableNIErbart = dbartresultBBmediationPOSTrbart$NIE_result_mc)
# (dbartBBtableNDErbart = dbartresultBBmediationPOSTrbart$NDE_result_mc)
# (dbartBBtableATErbart = dbartresultBBmediationPOSTrbart$ATE_result_mc)
# 
# # BB-HBB
# dbartresultHBBmediationPOSTrbart = HBBmediationPOST(rbartBARTfit, C, V, V0, esttype, saveall)
# (dbartHBBtableNIErbart = dbartresultHBBmediationPOSTrbart$NIE_result_mc)
# (dbartHBBtableNDErbart = dbartresultHBBmediationPOSTrbart$NDE_result_mc)
# (dbartHBBtableATErbart = dbartresultHBBmediationPOSTrbart$ATE_result_mc)
# 
# # TSBB
# dbartresultTSBBmediationPOSTrbart = TSBBmediationPOST(rbartBARTfit, C, V, V0, esttype, saveall, chi = 1, zeta = 0.5)
# (dbartTSBBtableNIErbart = dbartresultTSBBmediationPOSTrbart$NIE_result_mc)
# (dbartTSBBtableNDErbart = dbartresultTSBBmediationPOSTrbart$NDE_result_mc)
# (dbartTSBBtableATErbart = dbartresultTSBBmediationPOSTrbart$ATE_result_mc)
# c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])

source("TSBBmediation_source_r.R")
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# GLM
GLMfit = GLMmediation(Y, M, Z, C, V, V0,
                      typeY = Ydist, typeM = Mdist,
                      gibbs_iter = gibbs_iter, gibbs_burnin = gibbs_burnin, gibbs_thin = gibbs_thin)

# BB-BB
GLMresultBBmediationPOST = BBmediationPOST(GLMfit, C, V, V0, esttype, saveall, F)
(GLMBBtableNIE = GLMresultBBmediationPOST$NIE_result_mc)
(GLMBBtableNDE = GLMresultBBmediationPOST$NDE_result_mc)
(GLMBBtableATE = GLMresultBBmediationPOST$ATE_result_mc)

# BB-HBB
GLMresultHBBmediationPOST = HBBmediationPOST(GLMfit, C, V, V0, esttype, saveall, F)
(GLMHBBtableNIE = GLMresultHBBmediationPOST$NIE_result_mc)
(GLMHBBtableNDE = GLMresultHBBmediationPOST$NDE_result_mc)
(GLMHBBtableATE = GLMresultHBBmediationPOST$ATE_result_mc)

# # TSBB
# GLMresultTSBBmediationPOST = TSBBmediationPOST(GLMfit, C, V, V0, esttype, saveall, chi = 1, zeta = 0.5, F)
# (GLMTSBBtableNIE = GLMresultTSBBmediationPOST$NIE_result_mc)
# (GLMTSBBtableNDE = GLMresultTSBBmediationPOST$NDE_result_mc)
# (GLMTSBBtableATE = GLMresultTSBBmediationPOST$ATE_result_mc)
# c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])

GLMresultTSBBmediationPOSTsim = TSBBmediationPOSTsim(GLMfit, C, V, V0, esttype, saveall,
                                                     list_chi = c(1e-2, 1e-1, 1, 1e1, 1e2),
                                                     list_zeta = c(0, 0.5, 1), F)
# cbind(GLMresultTSBBmediationPOSTsim$NIE_result_mc,
#       GLMresultTSBBmediationPOSTsim$NDE_result_mc,
#       GLMresultTSBBmediationPOSTsim$ATE_result_mc)

round(rbind(c(GLMresultBBmediationPOST$NIE_result_mc,
              GLMresultBBmediationPOST$NDE_result_mc,
              GLMresultBBmediationPOST$ATE_result_mc),
            c(GLMresultHBBmediationPOST$NIE_result_mc,
              GLMresultHBBmediationPOST$NDE_result_mc,
              GLMresultHBBmediationPOST$ATE_result_mc),
            cbind(GLMresultTSBBmediationPOSTsim$NIE_result_mc,
                  GLMresultTSBBmediationPOSTsim$NDE_result_mc,
                  GLMresultTSBBmediationPOSTsim$ATE_result_mc)),3)
c(E_true[1],E_true[2],E_true[3])
c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])

# allinfo = c()
# list_chi = c(1e-2, 1e-1, 1, 1e1, 1e2)
# list_zeta = c(0, 0.5, 1)
# for (chi in list_chi){
#   for (zeta in list_zeta){
#     # ------------------------------------------------------------------------------
#     # Compputation steps -----------------------------------------------------------
#     # ------------------------------------------------------------------------------
#     # TSBB
#     resultTSBBmediationPOST = TSBBmediationPOST(GLMfit , C, V, V0, esttype, saveall, chi, zeta, T)
#     allinfo = c(allinfo,c(resultTSBBmediationPOST$NIE_result_mc,
#                           resultTSBBmediationPOST$NDE_result_mc,
#                           resultTSBBmediationPOST$ATE_result_mc))
#     rm(resultTSBBmediationPOST);gc()
#   }
# }
# allinfo = matrix(allinfo, nrow = 15)
# allinfo = t(cbind(c(rBARTmediationBBtableNIE,
#                     rBARTmediationBBtableNDE,
#                     rBARTmediationBBtableATE),
#                   c(rBARTmediationHBBtableNIE,
#                     rBARTmediationHBBtableNDE,
#                     rBARTmediationHBBtableATE),
#                   allinfo));allinfo

round(rbind(c(rBARTmediationresultBBmediationPOST$NIE_result_mc,
              rBARTmediationresultBBmediationPOST$NDE_result_mc,
              rBARTmediationresultBBmediationPOST$ATE_result_mc),
            c(rBARTmediationresultHBBmediationPOST$NIE_result_mc,
              rBARTmediationresultHBBmediationPOST$NDE_result_mc,
              rBARTmediationresultHBBmediationPOST$ATE_result_mc),
            cbind(rBARTmediationresultTSBBmediationPOSTsim$NIE_result_mc,
                  rBARTmediationresultTSBBmediationPOSTsim$NDE_result_mc,
                  rBARTmediationresultTSBBmediationPOSTsim$ATE_result_mc)),3)

round(rbind(c(GLMresultBBmediationPOST$NIE_result_mc,
              GLMresultBBmediationPOST$NDE_result_mc,
              GLMresultBBmediationPOST$ATE_result_mc),
            c(GLMresultHBBmediationPOST$NIE_result_mc,
              GLMresultHBBmediationPOST$NDE_result_mc,
              GLMresultHBBmediationPOST$ATE_result_mc),
            cbind(GLMresultTSBBmediationPOSTsim$NIE_result_mc,
                  GLMresultTSBBmediationPOSTsim$NDE_result_mc,
                  GLMresultTSBBmediationPOSTsim$ATE_result_mc)),3)

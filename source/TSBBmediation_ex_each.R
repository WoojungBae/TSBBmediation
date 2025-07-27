# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Load R packages
library(MASS)
library(Rcpp)
library(RcppArmadillo)

library(rBARTmediation)
# library(dbarts)
library(rstan)
# options(mc.cores = parallel::detectCores())
# rstan_options(auto_write = TRUE)
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
Scenario = 7
run_ID = 1
set.seed(run_ID)

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
  # n_j = 50
  n_j = 100
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
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Scenario
  #  1 -  6: Ydist - "continuous" & Mdist - "binary"
  #  7 - 12: Ydist - "continuous" & Mdist - "continuous"
  #  1 -  3 &   7 -  9: simple mean for Y
  #  4 -  6 &  10 - 12: complex mean for Y
  #  1, 4, 7, 10: CVdependence = 1 # "indep"
  #  2, 5, 8, 11: CVdependence = 2 # "weak"
  #  3, 6, 9, 12: CVdependence = 3 # "strong"
  
  if (Scenario==1) {
    # ------------------------------------------------------------------------------
    Ydist = "continuous"
    Mdist = "binary"
    simple = TRUE
    CVdependence = 1 # "indep"
  } else if (Scenario==2) {
    Ydist = "continuous"
    Mdist = "binary"
    simple = TRUE
    CVdependence = 2 # "weak"
  } else if (Scenario==3) {
    Ydist = "continuous"
    Mdist = "binary"
    simple = TRUE
    CVdependence = 3 # "strong"
  } else if (Scenario==4) {
    # ------------------------------------------------------------------------------
    Ydist = "continuous"
    Mdist = "binary"
    simple = FALSE
    CVdependence = 1 # "indep"
  } else if (Scenario==5) {
    Ydist = "continuous"
    Mdist = "binary"
    simple = FALSE
    CVdependence = 2 # "weak"
  } else if (Scenario==6) {
    Ydist = "continuous"
    Mdist = "binary"
    simple = FALSE
    CVdependence = 3 # "strong"
  } else if (Scenario==7) {
    # ------------------------------------------------------------------------------
    Ydist = "continuous"
    Mdist = "continuous"
    simple = TRUE
    CVdependence = 1 # "indep"
  } else if (Scenario==8) {
    Ydist = "continuous"
    Mdist = "continuous"
    simple = TRUE
    CVdependence = 2 # "weak"
  } else if (Scenario==9) {
    Ydist = "continuous"
    Mdist = "continuous"
    simple = TRUE
    CVdependence = 3 # "strong"
  } else if (Scenario==10) {
    # ------------------------------------------------------------------------------
    Ydist = "continuous"
    Mdist = "continuous"
    simple = FALSE
    CVdependence = 1 # "indep"
  } else if (Scenario==11) {
    Ydist = "continuous"
    Mdist = "continuous"
    simple = FALSE
    CVdependence = 2 # "weak"
  } else if (Scenario==12) {
    Ydist = "continuous"
    Mdist = "continuous"
    simple = FALSE
    CVdependence = 3 # "strong"
  } else if (Scenario==13) {
    # ------------------------------------------------------------------------------
    Ydist = "binary"
    Mdist = "binary"
    simple = TRUE
    CVdependence = 1 # "indep"
  } else if (Scenario==14) {
    Ydist = "binary"
    Mdist = "binary"
    simple = TRUE
    CVdependence = 2 # "weak"
  } else if (Scenario==15) {
    Ydist = "binary"
    Mdist = "binary"
    simple = TRUE
    CVdependence = 3 # "strong"
  } else if (Scenario==16) {
    # ------------------------------------------------------------------------------
    Ydist = "binary"
    Mdist = "binary"
    simple = FALSE
    CVdependence = 1 # "indep"
  } else if (Scenario==17) {
    Ydist = "binary"
    Mdist = "binary"
    simple = FALSE
    CVdependence = 2 # "weak"
  } else if (Scenario==18) {
    Ydist = "binary"
    Mdist = "binary"
    simple = FALSE
    CVdependence = 3 # "strong"
  }
  
  # ------------------------------------------------------------------------------
  # Data Generation --------------------------------------------------------------
  # ------------------------------------------------------------------------------
  
  # Data Generation --------------------------------------------------------------
  temp_data = generate_data(J,N,Scenario)
  
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
  Uindex = temp_data$Uindex
  
  # Causal Effects
  E_true = temp_data$E_true
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  uniqueUindex0 = sort(unique(Uindex[Z==0]))
  J0 = length(uniqueUindex0)
  uniqueUindex1 = sort(unique(Uindex[Z==1]))
  J1 = length(uniqueUindex1)
  UindexNEW = numeric(N)
  for (j in 1:J0) {
    whichj = which(Uindex == uniqueUindex0[j])
    UindexNEW[whichj] = j + J
  }
  for (j in 1:J1) {
    whichj = which(Uindex == uniqueUindex1[j])
    UindexNEW[whichj] = j + J0 + J
  }
  UindexNEW = UindexNEW - J
  Uindex = UindexNEW
  
  # Define the number of clusters
  uniqueUindex = sort(unique(Uindex))
  J = length(uniqueUindex)
  Uindex = apply(sapply(1:J, function(l) ifelse(Uindex==uniqueUindex[l],l,0)),1,sum)
  
  # Define the number of observations
  N = length(Y)
  
  # Define the number of observations for each cluster
  n_j = as.numeric(table(Uindex))
  
  # Data reordering
  order_Uindex = order(Uindex)
  Y = Y[order_Uindex]
  M = M[order_Uindex]
  Z = Z[order_Uindex]
  C = C[order_Uindex,]
  V = V[order_Uindex,]
  Uindex = Uindex[order_Uindex]
  
  colnames(C) = paste0("C",1:ncol(C))
  colnames(V) = paste0("V",1:ncol(V))
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  Z0 = 0
  which_Z0 = which(Z == Z0)
  Y0 = Y[which_Z0]
  M0 = M[which_Z0]
  C0 = C[which_Z0,]
  V0 = V[which_Z0,]
  Uindex0 = Uindex[which_Z0]
  matX0 = cbind(C0, V0)
  matM0 = cbind(M0, C0, V0)
  N0 = length(Y0)
  uniqueUindex0 = sort(unique(Uindex0))
  J0 = length(uniqueUindex0)
  Uindex0 = apply(sapply(1:J0, function(l) ifelse(Uindex0==uniqueUindex0[l],l,0)),1,sum)
  n_j0 = as.numeric(table(Uindex0))
  order_Uindex0 = order(Uindex0)
  Y0 = Y0[order_Uindex0]
  M0 = M0[order_Uindex0]
  C0 = C0[order_Uindex0,]
  V0 = V0[order_Uindex0,]
  Uindex0 = Uindex0[order_Uindex0]
  # df0 = data.frame(Y = Y0, M = M0, Z = Z0, C = C0, V = V0, Uindex = Uindex0, matX = matX0, matM = matM0)
  
  Z1 = 1
  which_Z1 = which(Z == Z1)
  Y1 = Y[which_Z1]
  M1 = M[which_Z1]
  C1 = C[which_Z1,]
  V1 = V[which_Z1,]
  Uindex1 = Uindex[which_Z1]
  matX1 = cbind(C1, V1)
  matM1 = cbind(M1, C1, V1)
  N1 = length(Y1)
  uniqueUindex1 = sort(unique(Uindex1))
  J1 = length(uniqueUindex1)
  Uindex1 = apply(sapply(1:J1, function(l) ifelse(Uindex1==uniqueUindex1[l],l,0)),1,sum)
  n_j1 = as.numeric(table(Uindex1))
  order_Uindex1 = order(Uindex1)
  Y1 = Y1[order_Uindex1]
  M1 = M1[order_Uindex1]
  C1 = C1[order_Uindex1,]
  V1 = V1[order_Uindex1,]
  Uindex1 = Uindex1[order_Uindex1]
  # df1 = data.frame(Y = Y1, M = M1, Z = Z1, C = C1, V = V1, Uindex = Uindex1, matX = matX1, matM = matM1)
}

gibbs_thin = 1e1
gibbs_iter = 2e3
gibbs_burnin = 2e4

# gibbs_thin = 1e1
# gibbs_iter = 2e2
# gibbs_burnin = 2e3

# sparse = TRUE
sparse = FALSE
# ntree = 20
# ntree = 50
# ntree = 100
ntree = 200
# ntree = 300
# ntree = 500
# ntree = 2000

c(E_true[1],E_true[2],E_true[3])
c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])

# meanY = mean(Y)
# meanYj = sapply(1:J, function(l) mean(Y[Uindex==l]))
# aY = sapply(1:N, function(i) Y[i] - meanYj[Uindex[i]])

# gibbs_thin = 1e1
# gibbs_iter = 1e5
# gibbs_burnin = 1e5
# gibbs_iter = gibbs_iter/gibbs_thin

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
source("TSBBmediation_source_r.R")
sourceCpp("TSBBmediation_source_cpp.cpp")
# rBARTmediation
BARTfit0 = rBARTmediation(Y0, M0, C0, V0, Uindex0,
                          typeY = Ydist, typeM = Mdist, 
                          ndpost=gibbs_iter, nskip=gibbs_burnin, keepevery=gibbs_thin,
                          ntree=ntree, sparse = sparse)

# rBARTmediation
BARTfit1 = rBARTmediation(Y1, M1, C1, V1, Uindex1,
                          typeY = Ydist, typeM = Mdist, 
                          ndpost=gibbs_iter, nskip=gibbs_burnin, keepevery=gibbs_thin, 
                          ntree=ntree, sparse = sparse)

BARTfit = list(object0 = BARTfit0, object1 = BARTfit1)
class(BARTfit) = "rBARTmediation"

c(BARTfit$object0$B_uM,BARTfit$object0$B_uY,
  BARTfit$object1$B_uM,BARTfit$object1$B_uY)

source("TSBBmediation_source_r.R")
sourceCpp("TSBBmediation_source_cpp.cpp")
# BB-BB
rBARTmediationresultBBmediationPOST = BBmediationPOST(BARTfit, C, V, Uindex, esttype, saveall, T)
(rBARTmediationBBtableNIE = rBARTmediationresultBBmediationPOST$NIE_result_mc)
(rBARTmediationBBtableNDE = rBARTmediationresultBBmediationPOST$NDE_result_mc)
(rBARTmediationBBtableATE = rBARTmediationresultBBmediationPOST$ATE_result_mc)

# BB-HBB
rBARTmediationresultHBBmediationPOST = HBBmediationPOST(BARTfit, C, V, Uindex, esttype, saveall, T)
(rBARTmediationHBBtableNIE = rBARTmediationresultHBBmediationPOST$NIE_result_mc)
(rBARTmediationHBBtableNDE = rBARTmediationresultHBBmediationPOST$NDE_result_mc)
(rBARTmediationHBBtableATE = rBARTmediationresultHBBmediationPOST$ATE_result_mc)

# # TSBB
# rBARTmediationresultTSBBmediationPOST = TSBBmediationPOST(BARTfit, C, V, Uindex, esttype, saveall, chi = 1, zeta = 0.5, T)
# (rBARTmediationTSBBtableNIE = rBARTmediationresultTSBBmediationPOST$NIE_result_mc)
# (rBARTmediationTSBBtableNDE = rBARTmediationresultTSBBmediationPOST$NDE_result_mc)
# (rBARTmediationTSBBtableATE = rBARTmediationresultTSBBmediationPOST$ATE_result_mc)

source("TSBBmediation_source_r.R")
sourceCpp("TSBBmediation_source_cpp.cpp")
rBARTmediationresultTSBBmediationPOSTsim = TSBBmediationPOSTsim(BARTfit, C, V, Uindex, esttype, saveall,
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

c(BARTfit$object0$B_uM,BARTfit$object0$B_uY,
  BARTfit$object1$B_uM,BARTfit$object1$B_uY)

# allinfo = c()
# list_chi = c(1e-2, 1e-1, 1, 1e1, 1e2)
# list_zeta = c(0, 0.5, 1)
# for (chi in list_chi){
#   for (zeta in list_zeta){
#     # ------------------------------------------------------------------------------
#     # Compputation steps -----------------------------------------------------------
#     # ------------------------------------------------------------------------------
#     # TSBB
#     resultTSBBmediationPOST = TSBBmediationPOST(BARTfit , C, V, Uindex, esttype, saveall, chi, zeta, T)
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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
source("TSBBmediation_source_r.R")
# source("TSBBmediation_source_r_includeZ.R")
# source("TSBBmediation_source_r_excludeZ.R")

# GLM
GLMfit = GLMmediation(Y, M, Z, C, V, Uindex,
                      typeY = Ydist, typeM = Mdist,
                      gibbs_iter = gibbs_iter * gibbs_thin, 
                      gibbs_burnin = gibbs_burnin, gibbs_thin = gibbs_thin)

# BB-BB
GLMresultBBmediationPOST = BBmediationPOST(GLMfit, C, V, Uindex, esttype, saveall, T)
(GLMBBtableNIE = GLMresultBBmediationPOST$NIE_result_mc)
(GLMBBtableNDE = GLMresultBBmediationPOST$NDE_result_mc)
(GLMBBtableATE = GLMresultBBmediationPOST$ATE_result_mc)

# BB-HBB
GLMresultHBBmediationPOST = HBBmediationPOST(GLMfit, C, V, Uindex, esttype, saveall, T)
(GLMHBBtableNIE = GLMresultHBBmediationPOST$NIE_result_mc)
(GLMHBBtableNDE = GLMresultHBBmediationPOST$NDE_result_mc)
(GLMHBBtableATE = GLMresultHBBmediationPOST$ATE_result_mc)

# # TSBB
# GLMresultTSBBmediationPOST = TSBBmediationPOST(GLMfit, C, V, Uindex, esttype, saveall, chi = 1, zeta = 0.5, T)
# (GLMTSBBtableNIE = GLMresultTSBBmediationPOST$NIE_result_mc)
# (GLMTSBBtableNDE = GLMresultTSBBmediationPOST$NDE_result_mc)
# (GLMTSBBtableATE = GLMresultTSBBmediationPOST$ATE_result_mc)
# c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])

GLMresultTSBBmediationPOSTsim = TSBBmediationPOSTsim(GLMfit, C, V, Uindex, esttype, saveall,
                                                     list_chi = c(1e-2, 1e-1, 1, 1e1, 1e2),
                                                     list_zeta = c(0, 0.5, 1), T)
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
#     resultTSBBmediationPOST = TSBBmediationPOST(GLMfit , C, V, Uindex, esttype, saveall, chi, zeta, T)
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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# # -----------------------------------------------------------------------------
# # -----------------------------------------------------------------------------
# # dbarts
# Mpost = rbart_vi(M ~ matX, df, group.by = Uindex,
#                  n.samples = gibbs_iter, n.burn = gibbs_burnin, n.thin = gibbs_thin,
#                  n.chains = 1, n.threads = 1)
# Ypost = rbart_vi(Y ~ matM, df, group.by = Uindex,
#                  n.samples = gibbs_iter, n.burn = gibbs_burnin, n.thin = gibbs_thin,
#                  n.chains = 1, n.threads = 1)
# rbartBARTfit = list(Yobject = Ypost, Mobject = Mpost)
# class(rbartBARTfit) = class(Ypost)
# 
# # BB-BB
# dbartresultBBmediationPOSTrbart = BBmediationPOST(rbartBARTfit, C, V, Uindex, esttype, saveall)
# (dbartBBtableNIErbart = dbartresultBBmediationPOSTrbart$NIE_result_mc)
# (dbartBBtableNDErbart = dbartresultBBmediationPOSTrbart$NDE_result_mc)
# (dbartBBtableATErbart = dbartresultBBmediationPOSTrbart$ATE_result_mc)
# 
# # BB-HBB
# dbartresultHBBmediationPOSTrbart = HBBmediationPOST(rbartBARTfit, C, V, Uindex, esttype, saveall)
# (dbartHBBtableNIErbart = dbartresultHBBmediationPOSTrbart$NIE_result_mc)
# (dbartHBBtableNDErbart = dbartresultHBBmediationPOSTrbart$NDE_result_mc)
# (dbartHBBtableATErbart = dbartresultHBBmediationPOSTrbart$ATE_result_mc)
# 
# # TSBB
# dbartresultTSBBmediationPOSTrbart = TSBBmediationPOST(rbartBARTfit, C, V, Uindex, esttype, saveall, chi = 1, zeta = 0.5)
# (dbartTSBBtableNIErbart = dbartresultTSBBmediationPOSTrbart$NIE_result_mc)
# (dbartTSBBtableNDErbart = dbartresultTSBBmediationPOSTrbart$NDE_result_mc)
# (dbartTSBBtableATErbart = dbartresultTSBBmediationPOSTrbart$ATE_result_mc)
# c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])

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

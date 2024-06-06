# -----------------------------------------------------------------------------
# Define functions ------------------------------------------------------------
# -----------------------------------------------------------------------------
# logistic function
logit = function(x_temp) { 
  return(log(x_temp/(1-x_temp)))
}

# inverse-logistic function
invlogit = function(x_temp) { 
  return((1/(1+exp(-x_temp))))
}

# Define function of distance for matrix
matrix_dist = function(X,V0){
  J = max(V0)
  n_j = as.numeric(table(V0))
  X_V0 = mapply(function(l) matrix(X[which(V0==l),], nrow=sum(V0==l)) ,1:J, SIMPLIFY = FALSE)
  
  X_dist = matrix(NA,nrow=J,ncol=J)
  for (l in 1:J) {
    n_k_temp = n_j[l]
    for (j in 1:J) {
      X_dist[l,j] = mean(sapply(1:n_k_temp, function(i) mean((X_V0[[l]][i,] - t(X_V0[[j]]))^{2})))
    }
  }
  
  return(X_dist)
}

# -----------------------------------------------------------------------------
# Define function to generate datasets
generate_data = function(J,       # the number of cluster
                         N,       # the number of total subjects
                         Scenario,
                         Mdist = c("Continuous", "Binary", "Ordinal", "Count"),
                         K = 7){  # the number of ordinal mediator when Mdist == "Ordinal"
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Generate cluster, V0cat
  pi_V0 = rep(1,J)
  V0cat = t(rmultinom(N,1,pi_V0))
  V0 = apply(V0cat,1,function(l) which(l == 1))
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Generate cluster-level treatment, Z
  Z = ifelse(V0 <= (J+1)/2,0,1)
  
  z0=0
  z1=1
  ind_z0=(Z==z0)
  ind_z1=(Z==z1)
  n0=sum(ind_z0)
  n1=sum(ind_z1)
  
  prob_Y1=0.6
  I_Y=numeric(N)
  I_Y[ind_z0]=rbinom(n0,1,prob_Y1)
  I_Y[ind_z1]=rbinom(n1,1,prob_Y1)
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Generate cluster-level and individual-level confounders (V, C)
  
  # ------------------------------------------------------------------------------
  # Generate cluster-level confounders, V
  VboolFalse = F
  while (VboolFalse==F) {
    # p_V = 3
    pi_V1 = rep(0.5, J)
    V1 = rbinom(J,1,pi_V1)
    pi_V2 = 0.4+0.2*V1
    V2 = rbinom(J,1,pi_V2)
    VboolFalse = (rank_cpp(cbind(Z, V1, V2))==3)
  }
  mu_V3 = -1 + V1 + 2*V2; sig2_V3 = rep(4, J)
  V3 = rnorm(J,mu_V3,sqrt(sig2_V3))
  V = cbind(V1 , V2, V3)
  V = t(sapply(1:N, function(i) V[V0[i],]))
  
  matV = cbind(1,V)
  
  # ------------------------------------------------------------------------------
  # Generate individual-level confounders, C
  p1_C = 3
  p2_C = 3
  p_C = p1_C + p2_C
  if (Scenario == 1) {
    beta_CV1 = rbind(rep(-1,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
    beta_CV2 = rbind(rep( 1,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
    beta_CV3 = rbind(rep( 0,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
    beta_CV4 = rbind(rep(-2,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
    beta_CV5 = rbind(rep( 2,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
    beta_CV6 = rbind(rep( 0,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
  } else if (Scenario == 2) {
    timesCV = 5
    beta_CV1 = rbind(-0.5*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV2 = rbind( 0.2*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV3 = rbind( 0.3*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV4 = rbind(-1.0*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV5 = rbind( 0.4*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV6 = rbind( 0.6*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
  } else if (Scenario == 3) {
    timesCV = 10
    beta_CV1 = rbind(-0.5*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV2 = rbind( 0.2*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV3 = rbind( 0.3*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV4 = rbind(-1.0*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV5 = rbind( 0.4*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV6 = rbind( 0.6*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
  }
  
  SIG_p2_C  = 0.7 * diag(p2_C) + 0.3 * rep(1,p2_C) %*% t(rep(1,p2_C))
  C1 = numeric(N)
  C2 = numeric(N)
  C3 = numeric(N)
  C_p2 = matrix(NA, nrow=N, ncol=p2_C)
  for (j in 1:J) {
    ind_V0 = which(V0 == j)
    num_V0 = length(ind_V0)
    
    matV_temp = matV[ind_V0,]
    
    pi_C1 = invlogit(matV_temp %*% beta_CV1[,j])
    C1[ind_V0] = rbinom(num_V0,1,pi_C1)
    
    pi_C2 = invlogit(matV_temp %*% beta_CV2[,j])
    C2[ind_V0] = rbinom(num_V0,1,pi_C2)
    
    pi_C3 = invlogit(matV_temp %*% beta_CV3[,j])
    C3[ind_V0] = rbinom(num_V0,1,pi_C3)
    
    MU_p2C = cbind(matV_temp %*% beta_CV4[,j],
                   matV_temp %*% beta_CV5[,j],
                   matV_temp %*% beta_CV6[,j])
    C_p2[ind_V0,] = t(sapply(1:num_V0, function(i) rmvn_cpp(1,MU_p2C[i,], SIG_p2_C)))
  }
  C = cbind(C1,C2,C3,C_p2)
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Generate random effects for M and Y
  # psi_M =  0.5*(-1)^{1:J}
  # psi_Y = -0.5*(-1)^{1:J}
  psi_M = rep(0,J)
  psi_Y = rep(0,J)
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  
  if (Mdist == "Continuous") {
    # Generate continuous mediator M
    beta_M0 =  0.5
    beta_MZ =  0.5
    beta_MV = c(-0.1, -0.1,  0.2)
    beta_MC = c( 0.1, -0.1,  0.1, -0.1,  0.1, -0.1)
    
    Mhat_z0 = numeric(N)
    Mhat_z1 = numeric(N)
    for (j in 1:J) {
      ind_V0 = which(V0 == j)
      num_V0 = length(ind_V0)
      C_temp = C[ind_V0,]
      V_temp = V[ind_V0,]
      Mhat_z0[ind_V0] = beta_M0 + z0 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j]
      Mhat_z1[ind_V0] = beta_M0 + z1 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j]
    }
    M_z0 = rnorm(N,Mhat_z0)
    M_z1 = rnorm(N,Mhat_z1)
    
  } else if (Mdist == "Binary") {
    # ------------------------------------------------------------------------------
    # Generate binary mediator M
    
  } else if (Mdist == "Ordinal") {
    # ------------------------------------------------------------------------------
    # Generate ordinal mediator M
    gamma_M0 = c(-2, 2, 4, 6, 8)
    beta_M0 =  1.0
    beta_MZ =  2.0
    beta_MV = c(-0.5, -0.5,  1.0)
    beta_MC = c( 0.5, -0.5,  0.5, -0.5,  0.5, -0.5)
    
    MboolFalse = F
    while (MboolFalse==F) {
      M_z0 = numeric(N)
      M_z1 = numeric(N)
      for (j in 1:J) {
        ind_V0 = which(V0 == j)
        num_V0 = length(ind_V0)
        
        C_temp = C[ind_V0,]
        V_temp = V[ind_V0,]
        
        M_z0_temp = sapply(1:(K-1), function(l) gamma_M0[l] - (beta_M0 + z0 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j]))
        M_z0_temp = t(apply(cbind(0,invlogit(M_z0_temp),1), 1, diff))
        M_z1_temp = sapply(1:(K-1), function(l) gamma_M0[l] - (beta_M0 + z1 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j]))
        M_z1_temp = t(apply(cbind(0,invlogit(M_z1_temp),1), 1, diff))
        
        M_z0[ind_V0] = sapply(1:num_V0, function(l) which(rmultinom(1, 1, M_z0_temp[l,])==1))
        M_z1[ind_V0] = sapply(1:num_V0, function(l) which(rmultinom(1, 1, M_z1_temp[l,])==1))
      }
      
      M = ifelse(Z==z0, M_z0, M_z1)
      uniqueM = sort(unique(M))
      if (length(uniqueM) < K) {
        M = sapply(1:N, function(l) which(uniqueM==M[l]))
      }
      
      matX = cbind(Z, C, V) # no intercept for proportional odds model
      factorM = ordered(M, levels = paste0(1:K))
      glm_M = try(polr(factorM ~ matX, Hess=TRUE), silent = TRUE)
      MboolFalse = !inherits(glm_M, "try-error")
    }
  } else if (Mdist == "Count") {
    # ------------------------------------------------------------------------------
    # Generate count mediator M
    beta_M0 =  0.5
    beta_MZ =  0.5
    beta_MV = c(-0.1, -0.1,  0.2)
    beta_MC = c( 0.1, -0.1,  0.1, -0.1,  0.1, -0.1)
    
    Mhat_z0 = numeric(N)
    Mhat_z1 = numeric(N)
    for (j in 1:J) {
      ind_V0 = which(V0 == j)
      num_V0 = length(ind_V0)
      C_temp = C[ind_V0,]
      V_temp = V[ind_V0,]
      Mhat_z0[ind_V0] = exp(beta_M0 + z0 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j])
      Mhat_z1[ind_V0] = exp(beta_M0 + z1 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j])
    }
    M_z0 = rpois(N,Mhat_z0)
    M_z1 = rpois(N,Mhat_z1)
  }
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Generate continuous outcome Y
  sig2_Y = 10^{2}
  beta_Y0 = -20
  beta_YM =  0.5
  beta_YZ =  2.0
  beta_YV = c(-0.5, -0.5,  1.0)
  beta_YC = c( 0.5, -0.5,  0.5, -0.5,  0.5, -0.5)
  
  # cbind(1,M,C,Z,V,1)
  # intercept, mediator,individual-level, treatment, culster-level, random-effect
  Y_z0m0 = numeric(N)
  Y_z1m0 = numeric(N)
  Y_z1m1 = numeric(N)
  for (j in 1:J) {
    ind_V0 = which(V0 == j)
    num_V0 = length(ind_V0)
    
    C_temp = C[ind_V0,]
    V_temp = V[ind_V0,]
    M_z0_temp = M_z0[ind_V0]
    M_z1_temp = M_z1[ind_V0]
    
    Y_z0m0_temp = beta_Y0 + M_z0_temp * beta_YM + z0 * beta_YZ + C_temp %*% beta_YC + V_temp %*% beta_YV + psi_Y[j]
    Y_z1m0_temp = beta_Y0 + M_z0_temp * beta_YM + z1 * beta_YZ + C_temp %*% beta_YC + V_temp %*% beta_YV + psi_Y[j]
    Y_z1m1_temp = beta_Y0 + M_z1_temp * beta_YM + z1 * beta_YZ + C_temp %*% beta_YC + V_temp %*% beta_YV + psi_Y[j]
    
    Y_z0m0[ind_V0] = rnorm(num_V0, Y_z0m0_temp, sqrt(sig2_Y))
    Y_z1m0[ind_V0] = rnorm(num_V0, Y_z1m0_temp, sqrt(sig2_Y))
    Y_z1m1[ind_V0] = rnorm(num_V0, Y_z1m1_temp, sqrt(sig2_Y))
  }
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  M = ifelse(Z==z0, M_z0, M_z1)
  Y = ifelse(Z==z0, Y_z0m0, Y_z1m1)
  
  E_true = c(mean(Y_z1m1), mean(Y_z1m0), mean(Y_z0m0))
  # E_true = c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])
  
  return(list(Y=Y, M=M, C=C, Z=Z, V=V, V0=V0, E_true=E_true))
}

# -----------------------------------------------------------------------------
BBmediationPOST_each = function(object0, # object from rBARTmediation
                                object1, # object from rBARTmediation
                                X.test,	# matrix X to predict at
                                Uindex0,
                                Uindex1,
                                esttype = "mean",
                                saveall = FALSE){
  
  # predict(rBARTmediation)
  BARTfitPRED = prBARTmediation_each(object0, object1, X.test, Uindex0, Uindex1)
  
  # Constants
  N = nrow(X.test)
  n_MCMC = nrow(BARTfitPRED$Yz0m0.test)
  
  # Significance Level alpha
  level = 0.05
  
  # Parameters for individual-level confounder
  a_pi = rep(1,N)
  
  # Define interation check
  iter_check = floor(n_MCMC/10)
  
  # Causal Effects
  NIE = numeric(n_MCMC)
  NDE = numeric(n_MCMC)
  ATE = numeric(n_MCMC)
  
  for (post_reps in 1:n_MCMC) {
    # Update individual&culster-level confounders parameters
    piPars = c(rdirichlet_cpp(1,a_pi))
    
    # 
    Yhat_z0m0_mc = BARTfitPRED$Yz0m0.test[post_reps,]
    Yhat_z1m0_mc = BARTfitPRED$Yz1m0.test[post_reps,]
    Yhat_z1m1_mc = BARTfitPRED$Yz1m1.test[post_reps,]
    E_Y_z0m0_mc  = sum(piPars * Yhat_z0m0_mc)
    E_Y_z1m0_mc  = sum(piPars * Yhat_z1m0_mc)
    E_Y_z1m1_mc  = sum(piPars * Yhat_z1m1_mc)
    
    # 
    NIE_mc = E_Y_z1m1_mc - E_Y_z1m0_mc
    NDE_mc = E_Y_z1m0_mc - E_Y_z0m0_mc
    ATE_mc = E_Y_z1m1_mc - E_Y_z0m0_mc
    NIE[post_reps] = mean(NIE_mc)
    NDE[post_reps] = mean(NDE_mc)
    ATE[post_reps] = mean(ATE_mc)
    
    if (post_reps %% iter_check == 0){
      cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  
  # Estimates: mean or median
  if (esttype == "median"){
    NIE_est_mc = median(NIE)
    NDE_est_mc = median(NDE)
    ATE_est_mc = median(ATE)
  } else {
    NIE_est_mc = mean(NIE)
    NDE_est_mc = mean(NDE)
    ATE_est_mc = mean(ATE)
  }
  
  quantile_alpha = c(level/2,1-level/2)
  
  # Calculate median of posterior for NIE
  NIE_sd_mc = sd(NIE)
  NIE_quantile_mc = quantile(NIE, quantile_alpha)
  NIE_CIlength_mc = abs(diff(NIE_quantile_mc))
  NIE_quantile025_mc = min(NIE_quantile_mc)
  NIE_quantile975_mc = max(NIE_quantile_mc)
  
  # Calculate median of posterior for NDE
  NDE_sd_mc = sd(NDE)
  NDE_quantile_mc = quantile(NDE, quantile_alpha)
  NDE_CIlength_mc = abs(diff(NDE_quantile_mc))
  NDE_quantile025_mc = min(NDE_quantile_mc)
  NDE_quantile975_mc = max(NDE_quantile_mc)
  
  # Calculate median of posterior for ATE
  ATE_sd_mc = sd(ATE)
  ATE_quantile_mc = quantile(ATE, quantile_alpha)
  ATE_CIlength_mc = abs(diff(ATE_quantile_mc))
  ATE_quantile025_mc = min(ATE_quantile_mc)
  ATE_quantile975_mc = max(ATE_quantile_mc)
  
  NIE_result_mc = cbind(NIE_est_mc, NIE_sd_mc, NIE_quantile025_mc, NIE_quantile975_mc, NIE_CIlength_mc)
  colnames(NIE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(NIE_result_mc) = c("")
  
  NDE_result_mc = cbind(NDE_est_mc, NDE_sd_mc, NDE_quantile025_mc, NDE_quantile975_mc, NDE_CIlength_mc)
  colnames(NDE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(NDE_result_mc) = c("")
  
  ATE_result_mc = cbind(ATE_est_mc, ATE_sd_mc, ATE_quantile025_mc, ATE_quantile975_mc, ATE_CIlength_mc)
  colnames(ATE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(ATE_result_mc) = c("")
  
  if (saveall) {
    POSTresult = list(NIE_result_mc = NIE_result_mc,
                      NDE_result_mc = NDE_result_mc,
                      ATE_result_mc = ATE_result_mc,
                      NIE = NIE, NDE = NDE, ATE = ATE)
  } else {
    POSTresult = list(NIE_result_mc = NIE_result_mc,
                      NDE_result_mc = NDE_result_mc,
                      ATE_result_mc = ATE_result_mc)
  }
  
  return(POSTresult)
}

# -----------------------------------------------------------------------------
HBBmediationPOST_each = function(object0, # object from rBARTmediation
                                 object1, # object from rBARTmediation
                                 X.test,	# matrix X to predict at
                                 Uindex0,
                                 Uindex1,
                                 esttype = "mean",
                                 saveall = FALSE){
  
  # predict(rBARTmediation)
  BARTfitPRED = prBARTmediation_each(object0, object1, X.test, Uindex0, Uindex1)
  
  # Constants
  N = nrow(X.test)
  n_j = as.numeric(table(Uindex))
  n_MCMC = nrow(BARTfitPRED$Yz0m0.test)
  
  # Significance Level alpha
  level = 0.05
  
  # Parameters for cluster-level confounder
  a_rho = rep(1,J)
  
  # Parameters for individual-level confounder
  a_pi = rep(1,N)
  
  # Parameters for alpha (concentration parameter)
  tau_ome = 100
  alpha_ome = (N/n_j)*tau_ome
  
  # Define interation check
  iter_check = floor(n_MCMC/10)
  
  # Causal Effects
  NIE = numeric(n_MCMC)
  NDE = numeric(n_MCMC)
  ATE = numeric(n_MCMC)
  
  for (post_reps in 1:n_MCMC) {
    # Update culster-level confounders parameters
    VrhoPars = c(rdirichlet_cpp(1,a_rho))
    
    # Update individual-level confounders parameters
    CpiPars = c(rdirichlet_cpp(1,a_pi))
    
    # Update the parameters for individual-level confounders given culster-level confounders
    CpiJPars = sapply(1:J, function(j) rdirichlet_cpp(1,(alpha_ome[j] * CpiPars + (Uindex==j))))
    
    # 
    Yhat_z0m0_mc = BARTfitPRED$Yz0m0.test[post_reps,]
    Yhat_z1m0_mc = BARTfitPRED$Yz1m0.test[post_reps,]
    Yhat_z1m1_mc = BARTfitPRED$Yz1m1.test[post_reps,]
    E_Y_z0m0_mc  = sum(VrhoPars * apply(CpiJPars * Yhat_z0m0_mc,2,sum))
    E_Y_z1m0_mc  = sum(VrhoPars * apply(CpiJPars * Yhat_z1m0_mc,2,sum))
    E_Y_z1m1_mc  = sum(VrhoPars * apply(CpiJPars * Yhat_z1m1_mc,2,sum))
    
    # 
    NIE_mc = E_Y_z1m1_mc - E_Y_z1m0_mc
    NDE_mc = E_Y_z1m0_mc - E_Y_z0m0_mc
    ATE_mc = E_Y_z1m1_mc - E_Y_z0m0_mc
    NIE[post_reps] = mean(NIE_mc)
    NDE[post_reps] = mean(NDE_mc)
    ATE[post_reps] = mean(ATE_mc)
    
    if (post_reps %% iter_check == 0){
      cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  
  # Estimates: mean or median
  if (esttype == "median"){
    NIE_est_mc = median(NIE)
    NDE_est_mc = median(NDE)
    ATE_est_mc = median(ATE)
  } else {
    NIE_est_mc = mean(NIE)
    NDE_est_mc = mean(NDE)
    ATE_est_mc = mean(ATE)
  }
  
  quantile_alpha = c(level/2,1-level/2)
  
  # Calculate median of posterior for NIE
  NIE_sd_mc = sd(NIE)
  NIE_quantile_mc = quantile(NIE, quantile_alpha)
  NIE_CIlength_mc = abs(diff(NIE_quantile_mc))
  NIE_quantile025_mc = min(NIE_quantile_mc)
  NIE_quantile975_mc = max(NIE_quantile_mc)
  
  # Calculate median of posterior for NDE
  NDE_sd_mc = sd(NDE)
  NDE_quantile_mc = quantile(NDE, quantile_alpha)
  NDE_CIlength_mc = abs(diff(NDE_quantile_mc))
  NDE_quantile025_mc = min(NDE_quantile_mc)
  NDE_quantile975_mc = max(NDE_quantile_mc)
  
  # Calculate median of posterior for ATE
  ATE_sd_mc = sd(ATE)
  ATE_quantile_mc = quantile(ATE, quantile_alpha)
  ATE_CIlength_mc = abs(diff(ATE_quantile_mc))
  ATE_quantile025_mc = min(ATE_quantile_mc)
  ATE_quantile975_mc = max(ATE_quantile_mc)
  
  NIE_result_mc = cbind(NIE_est_mc, NIE_sd_mc, NIE_quantile025_mc, NIE_quantile975_mc, NIE_CIlength_mc)
  colnames(NIE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(NIE_result_mc) = c("")
  
  NDE_result_mc = cbind(NDE_est_mc, NDE_sd_mc, NDE_quantile025_mc, NDE_quantile975_mc, NDE_CIlength_mc)
  colnames(NDE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(NDE_result_mc) = c("")
  
  ATE_result_mc = cbind(ATE_est_mc, ATE_sd_mc, ATE_quantile025_mc, ATE_quantile975_mc, ATE_CIlength_mc)
  colnames(ATE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(ATE_result_mc) = c("")
  
  if (saveall) {
    POSTresult = list(NIE_result_mc = NIE_result_mc,
                      NDE_result_mc = NDE_result_mc,
                      ATE_result_mc = ATE_result_mc,
                      NIE = NIE, NDE = NDE, ATE = ATE)
  } else {
    POSTresult = list(NIE_result_mc = NIE_result_mc,
                      NDE_result_mc = NDE_result_mc,
                      ATE_result_mc = ATE_result_mc)
  }
  
  return(POSTresult)
}

# -----------------------------------------------------------------------------
TSBBmediationPOST_each = function(object0, # object from rBARTmediation
                                  object1, # object from rBARTmediation
                                  Z.test,	# matrix X to predict at
                                  C.test,	# matrix X to predict at
                                  V.test,	# matrix X to predict at
                                  Uindex,
                                  esttype = "mean",
                                  saveall = FALSE, 
                                  chi = 1, zeta = 0.5){
  # a constant to scale the desired rate of decay
  # chi = 1
  # a constant ratio of cluster-level confounders and individual-level confounders
  # zeta = 0.5
  
  N = nrow(C.test)
  n_j = as.numeric(table(Uindex))
  J = length(unique(Uindex))
  p.test = ncol(C.test) + ncol(V.test)
  
  UindexTSBB = numeric(N * J)
  Z.testTSBB = numeric(N * J)
  X.testTSBB = matrix(nrow = N * J, ncol = p.test)
  for (l in 1:J) {
    ind_V0_l = which(Uindex == l)
    V_lunique = V.test[ind_V0_l[1],]
    Xtemp = matrix(nrow = N, ncol = p.test)
    for (j in 1:J) {
      ind_V0_j = which(Uindex == j)
      num_V0_j = length(ind_V0_j)
      C_j = matrix(C.test[ind_V0_j,], nrow = num_V0_j)
      V_l = matrix(rep(V_lunique,num_V0_j), nrow = num_V0_j, byrow = T)
      Xtemp[ind_V0_j,] = cbind(C_j, V_l)
    }
    Z.testTSBB[((N*(l-1)+1):(N*l))] = Z.test[ind_V0_l[1]]
    X.testTSBB[((N*(l-1)+1):(N*l)),] = Xtemp
  }
  
  # predict(rBARTmediation)
  BARTfitPRED = prBARTmediation_each(object0, object1, X.testTSBB, UindexTSBB[Z.test==0], UindexTSBB[Z.test==1])
  
  # Constants
  n_MCMC = nrow(BARTfitPRED$Yz0m0.test)
  
  # Significance Level alpha
  level = 0.05
  
  # Parameters for cluster-level confounder
  a_rho = rep(1,J)
  
  # Parameters for individual-level confounder
  a_pi = rep(1,N)
  
  # Parameters for alpha (concentration parameter)
  tau_ome = 100
  alpha_ome = (N/n_j)*tau_ome
  
  # Distances
  # \lambda_{lj}^{\oemga} = \alpha_{l}^{\oemga} \dist_{lj}: adding additional "pseudosubjects
  # tau_ome = 1
  # tau_ome = 100
  # alpha_ome = rep(tau_ome,J)
  # alpha_ome = n_j/N
  alpha_ome = N/n_j
  
  # j (row) times l (column)
  Vxi_dist = matrix_dist(scale(V.test),Uindex)
  Cxi_dist = matrix_dist(scale(C.test),Uindex)
  xi_dist = exp(-((zeta)*Vxi_dist + (1-zeta)*Cxi_dist)/chi)
  xi_distPOST = sapply(1:J, function(l) alpha_ome[l]*xi_dist[,l]) + diag(J)
  
  # Define interation check
  iter_check = floor(n_MCMC/10)
  
  # Causal Effects
  NIE = numeric(n_MCMC)
  NDE = numeric(n_MCMC)
  ATE = numeric(n_MCMC)
  
  for (post_reps in 1:n_MCMC) {
    # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
    piPars = sapply(1:J, function(l) rdirichlet_cpp(1, rep(1,n_j[l])))
    
    # Update parameters for link cluster-individual (omega: J X J matrix)
    omePars = sapply(1:J, function(l) rdirichlet_cpp(1, xi_distPOST[,l]))
    
    # Update parameters for cluster-level confounder (rho: J X 1 vector)
    rhoPars = c(rdirichlet_cpp(1,rep(1,J)))
    
    RhoOmegaPi = numeric(N * J)
    for (l in 1:J) {
      RhoOmegaPitemp = numeric(N)
      for (j in 1:J) {
        ind_V0_j = which(Uindex == j)
        RhoOmegaPitemp[ind_V0_j] = omePars[j,l]*c(piPars[[j]])
      }
      RhoOmegaPi[((N*(l-1)+1):(N*l))] = rhoPars[l]*RhoOmegaPitemp
    }
    
    # 
    Yhat_z0m0_mc = BARTfitPRED$Yz0m0.test[post_reps,]
    Yhat_z1m0_mc = BARTfitPRED$Yz1m0.test[post_reps,]
    Yhat_z1m1_mc = BARTfitPRED$Yz1m1.test[post_reps,]
    E_Y_z0m0_mc  = sum(RhoOmegaPi * Yhat_z0m0_mc)
    E_Y_z1m0_mc  = sum(RhoOmegaPi * Yhat_z1m0_mc)
    E_Y_z1m1_mc  = sum(RhoOmegaPi * Yhat_z1m1_mc)
    
    # calculate causal effects
    # NIE_mc = E_Y_z1m1_mc
    # NDE_mc = E_Y_z1m0_mc
    # ATE_mc = E_Y_z1m1_mc
    NIE_mc = E_Y_z1m1_mc - E_Y_z1m0_mc
    NDE_mc = E_Y_z1m0_mc - E_Y_z0m0_mc
    ATE_mc = E_Y_z1m1_mc - E_Y_z0m0_mc
    NIE[post_reps] = mean(NIE_mc)
    NDE[post_reps] = mean(NDE_mc)
    ATE[post_reps] = mean(ATE_mc)
    
    if (post_reps %% iter_check == 0){
      cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
    }
  }
  
  # Estimates: mean or median
  if (esttype == "median"){
    NIE_est_mc = median(NIE)
    NDE_est_mc = median(NDE)
    ATE_est_mc = median(ATE)
  } else {
    NIE_est_mc = mean(NIE)
    NDE_est_mc = mean(NDE)
    ATE_est_mc = mean(ATE)
  }
  
  quantile_alpha = c(level/2,1-level/2)
  
  # Calculate median of posterior for NIE
  NIE_sd_mc = sd(NIE)
  NIE_quantile_mc = quantile(NIE, quantile_alpha)
  NIE_CIlength_mc = abs(diff(NIE_quantile_mc))
  NIE_quantile025_mc = min(NIE_quantile_mc)
  NIE_quantile975_mc = max(NIE_quantile_mc)
  
  # Calculate median of posterior for NDE
  NDE_sd_mc = sd(NDE)
  NDE_quantile_mc = quantile(NDE, quantile_alpha)
  NDE_CIlength_mc = abs(diff(NDE_quantile_mc))
  NDE_quantile025_mc = min(NDE_quantile_mc)
  NDE_quantile975_mc = max(NDE_quantile_mc)
  
  # Calculate median of posterior for ATE
  ATE_sd_mc = sd(ATE)
  ATE_quantile_mc = quantile(ATE, quantile_alpha)
  ATE_CIlength_mc = abs(diff(ATE_quantile_mc))
  ATE_quantile025_mc = min(ATE_quantile_mc)
  ATE_quantile975_mc = max(ATE_quantile_mc)
  
  # Save the results
  NIE_result_mc = cbind(NIE_est_mc, NIE_sd_mc, NIE_quantile025_mc, NIE_quantile975_mc, NIE_CIlength_mc)
  colnames(NIE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(NIE_result_mc) = c("")
  
  NDE_result_mc = cbind(NDE_est_mc, NDE_sd_mc, NDE_quantile025_mc, NDE_quantile975_mc, NDE_CIlength_mc)
  colnames(NDE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(NDE_result_mc) = c("")
  
  ATE_result_mc = cbind(ATE_est_mc, ATE_sd_mc, ATE_quantile025_mc, ATE_quantile975_mc, ATE_CIlength_mc)
  colnames(ATE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(ATE_result_mc) = c("")
  
  if (saveall) {
    POSTresult = list(NIE_result_mc = NIE_result_mc,
                      NDE_result_mc = NDE_result_mc,
                      ATE_result_mc = ATE_result_mc,
                      NIE = NIE, NDE = NDE, ATE = ATE)
  } else {
    POSTresult = list(NIE_result_mc = NIE_result_mc,
                      NDE_result_mc = NDE_result_mc,
                      ATE_result_mc = ATE_result_mc)
  }
  
  return(POSTresult)
}

# -----------------------------------------------------------------------------
# End function definitions ----------------------------------------------------
# -----------------------------------------------------------------------------
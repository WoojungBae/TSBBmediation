# -----------------------------------------------------------------------------
# Define functions ------------------------------------------------------------
# -----------------------------------------------------------------------------
# Define function of distance for matrix
matrix_dist = function(X,Uindex){
  J = max(Uindex)
  n_j = as.numeric(table(Uindex))
  X_Uindex = mapply(function(l) matrix(X[which(Uindex==l),], nrow=sum(Uindex==l)) ,1:J, SIMPLIFY = FALSE)
  
  X_dist = matrix(NA,nrow=J,ncol=J)
  for (l in 1:J) {
    n_k_temp = n_j[l]
    for (j in 1:J) {
      X_dist[l,j] = mean(sapply(1:n_k_temp, function(i) mean((X_Uindex[[l]][i,] - t(X_Uindex[[j]]))^{2})))
    }
  }
  
  return(X_dist)
}

# -----------------------------------------------------------------------------
# Define function to generate datasets
generate_data = function(J,       # the number of cluster
                         N,       # the number of total subjects
                         Scenario){
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
  }
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  if (Scenario==13) {
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
  # ------------------------------------------------------------------------------
  # Generate cluster, Uindexcat
  Uindexcat = rmultinom(N, 1, rep(1,J))
  Uindex = apply(Uindexcat, 2, function(l) which(l == 1))
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Generate cluster-level treatment, Z
  Z_Uindex = sort(sample(J,(J/2)))
  Z = (Uindex %in% Z_Uindex) * 1
  z0 = 0
  z1 = 1
  
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
  V = t(sapply(1:N, function(i) V[Uindex[i],]))
  
  matV = cbind(1,V)
  
  # ------------------------------------------------------------------------------
  # Generate individual-level confounders, C
  p1_C = 3
  p2_C = 3
  p_C = p1_C + p2_C
  if (CVdependence == 1) {
    beta_CV1 = rbind(rep(-1,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
    beta_CV2 = rbind(rep( 1,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
    beta_CV3 = rbind(rep( 0,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
    beta_CV4 = rbind(rep(-2,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
    beta_CV5 = rbind(rep( 2,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
    beta_CV6 = rbind(rep( 0,J),matrix(rep(c( 0, 0, 0), J),ncol=J))
  } else if (CVdependence == 2) {
    timesCV = 10
    beta_CV1 = rbind(-0.5*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV2 = rbind( 0.2*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV3 = rbind( 0.3*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV4 = rbind(-1.0*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV5 = rbind( 0.4*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
    beta_CV6 = rbind( 0.6*(-1)^{1:J},matrix(rep(c( 0.3, -0.2, -0.1)*timesCV, J),ncol=J))
  } else if (CVdependence == 3) {
    timesCV = 20
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
    ind_Uindex = which(Uindex == j)
    num_Uindex = length(ind_Uindex)
    
    matV_temp = matV[ind_Uindex,]
    
    pi_C1 = plogis(matV_temp %*% beta_CV1[,j])
    C1[ind_Uindex] = rbinom(num_Uindex,1,pi_C1)
    
    pi_C2 = plogis(matV_temp %*% beta_CV2[,j])
    C2[ind_Uindex] = rbinom(num_Uindex,1,pi_C2)
    
    pi_C3 = plogis(matV_temp %*% beta_CV3[,j])
    C3[ind_Uindex] = rbinom(num_Uindex,1,pi_C3)
    
    MU_p2C = cbind(matV_temp %*% beta_CV4[,j],
                   matV_temp %*% beta_CV5[,j],
                   matV_temp %*% beta_CV6[,j])
    C_p2[ind_Uindex,] = t(sapply(1:num_Uindex, function(i) rmvn_cpp(1,MU_p2C[i,], SIG_p2_C)))
  }
  C = cbind(C1,C2,C3,C_p2)
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  # Generate random effects for M and Y
  B_psi_M = 1; B_psi_Y = 1; b_psi_YM = 0.3
  bB_psi_YM = b_psi_YM * sqrt(B_psi_M * B_psi_Y)
  psi_YM = rmvn_cpp(J, c(0,0),
                    rbind(c(B_psi_Y, bB_psi_YM),
                          c(bB_psi_YM, B_psi_M)))
  psi_M = psi_YM[1,]
  psi_Y = psi_YM[2,]
  # psi_M =  1*(-1)^{1:J}
  # psi_Y = -1*(-1)^{1:J}
  # psi_M = rep(0,J)
  # psi_Y = rep(0,J)
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  if (Mdist == "continuous") {
    # Generate continuous mediator M
    beta_M0 =  0.5
    beta_MZ =  0.5
    beta_MV = c(-0.1, -0.1,  0.2)
    beta_MC = c( 0.1, -0.1,  0.1, -0.1,  0.1, -0.1)
    
    Mhat_z0 = numeric(N)
    Mhat_z1 = numeric(N)
    for (j in 1:J) {
      ind_Uindex = which(Uindex == j)
      num_Uindex = length(ind_Uindex)
      C_temp = C[ind_Uindex,]
      V_temp = V[ind_Uindex,]
      Mhat_z0[ind_Uindex] = beta_M0 + z0 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j]
      Mhat_z1[ind_Uindex] = beta_M0 + z1 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j]
    }
    M_z0 = rnorm(N, Mhat_z0)
    M_z1 = rnorm(N, Mhat_z1)
    
  } else if (Mdist == "binary") {
    # ------------------------------------------------------------------------------
    # Generate binary mediator M
    beta_M0 = -1.0
    beta_MZ = -beta_M0*2
    beta_MV = c(-0.1, -0.1,  0.2) * 1
    beta_MC = c( 0.1, -0.1,  0.1, -0.1,  0.1, -0.1) * 1
    
    Mhat_z0 = numeric(N)
    Mhat_z1 = numeric(N)
    for (j in 1:J) {
      ind_Uindex = which(Uindex == j)
      num_Uindex = length(ind_Uindex)
      C_temp = C[ind_Uindex,]
      V_temp = V[ind_Uindex,]
      Mhat_z0[ind_Uindex] = beta_M0 + z0 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j]
      Mhat_z1[ind_Uindex] = beta_M0 + z1 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j]
    }
    M_z0 = rbinom(N, 1, pnorm(Mhat_z0))
    M_z1 = rbinom(N, 1, pnorm(Mhat_z1))
    
  } else if (Mdist == "ordinal") {
    # ------------------------------------------------------------------------------
    # Generate ordinal mediator M
    gamma_M0 = c(-2, 2, 4, 6, 8)
    beta_M0 =  1.0
    beta_MZ =  2.0
    beta_MV = c(-0.1, -0.1,  0.2) * 5
    beta_MC = c( 0.1, -0.1,  0.1, -0.1,  0.1, -0.1) * 5
    
    K = length(gamma_M0)
    
    MboolFalse = F
    while (MboolFalse==F) {
      M_z0 = numeric(N)
      M_z1 = numeric(N)
      for (j in 1:J) {
        ind_Uindex = which(Uindex == j)
        num_Uindex = length(ind_Uindex)
        
        C_temp = C[ind_Uindex,]
        V_temp = V[ind_Uindex,]
        
        M_z0_temp = sapply(1:(K-1), function(l) gamma_M0[l] - (beta_M0 + z0 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j]))
        M_z0_temp = t(apply(cbind(0,plogis(M_z0_temp),1), 1, diff))
        M_z1_temp = sapply(1:(K-1), function(l) gamma_M0[l] - (beta_M0 + z1 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j]))
        M_z1_temp = t(apply(cbind(0,plogis(M_z1_temp),1), 1, diff))
        
        M_z0[ind_Uindex] = sapply(1:num_Uindex, function(l) which(rmultinom(1, 1, M_z0_temp[l,])==1))
        M_z1[ind_Uindex] = sapply(1:num_Uindex, function(l) which(rmultinom(1, 1, M_z1_temp[l,])==1))
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
  } else if (Mdist == "count") {
    # ------------------------------------------------------------------------------
    # Generate count mediator M
    beta_M0 =  0.5
    beta_MZ =  1.0
    beta_MV = c(-0.1, -0.1,  0.2)
    beta_MC = c( 0.1, -0.1,  0.1, -0.1,  0.1, -0.1)
    
    Mhat_z0 = numeric(N)
    Mhat_z1 = numeric(N)
    for (j in 1:J) {
      ind_Uindex = which(Uindex == j)
      num_Uindex = length(ind_Uindex)
      C_temp = C[ind_Uindex,]
      V_temp = V[ind_Uindex,]
      Mhat_z0[ind_Uindex] = exp(beta_M0 + z0 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j])
      Mhat_z1[ind_Uindex] = exp(beta_M0 + z1 * beta_MZ + C_temp %*% beta_MC + V_temp %*% beta_MV + psi_M[j])
    }
    M_z0 = rpois(N,Mhat_z0)
    M_z1 = rpois(N,Mhat_z1)
  }
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  if (Ydist == "continuous"){
    # ------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------
    # Generate continuous outcome Y
    sig2_Y = 5^{2}
    beta_Y0 = -20
    beta_YM =  1.0
    beta_YZ =  1.0
    beta_YV = c(-0.5, -0.5,  1.0)
    beta_YC = c( 0.5, -0.5,  0.5, -0.5,  0.5, -0.5)
    
    # cbind(1,M,C,Z,V,1)
    # intercept, mediator,individual-level, treatment, culster-level, random-effect
    Y_z0m0 = numeric(N)
    Y_z1m0 = numeric(N)
    Y_z1m1 = numeric(N)
    for (j in 1:J) {
      ind_Uindex = which(Uindex == j)
      num_Uindex = length(ind_Uindex)
      
      C_temp = C[ind_Uindex,]
      V_temp = V[ind_Uindex,]
      M_z0_temp = M_z0[ind_Uindex]
      M_z1_temp = M_z1[ind_Uindex]
      
      Y_z0m0_temp = beta_Y0 + z0 * beta_YZ + C_temp %*% beta_YC + V_temp %*% beta_YV + psi_Y[j]
      Y_z1m0_temp = beta_Y0 + z1 * beta_YZ + C_temp %*% beta_YC + V_temp %*% beta_YV + psi_Y[j]
      Y_z1m1_temp = beta_Y0 + z1 * beta_YZ + C_temp %*% beta_YC + V_temp %*% beta_YV + psi_Y[j]
      if (simple) {
        Y_z0m0_temp = Y_z0m0_temp + M_z0_temp * beta_YM
        Y_z1m0_temp = Y_z1m0_temp + M_z0_temp * beta_YM
        Y_z1m1_temp = Y_z1m1_temp + M_z1_temp * beta_YM
      } else {
        Y_z0m0_temp = Y_z0m0_temp + 0.3 * (M_z0_temp * z0 + M_z0_temp^{2})
        Y_z1m0_temp = Y_z1m0_temp + 0.3 * (M_z0_temp * z1 + M_z0_temp^{2})
        Y_z1m1_temp = Y_z1m1_temp + 0.3 * (M_z1_temp * z1 + M_z1_temp^{2})
      }
      
      Y_z0m0[ind_Uindex] = rnorm(num_Uindex, Y_z0m0_temp, sqrt(sig2_Y))
      Y_z1m0[ind_Uindex] = rnorm(num_Uindex, Y_z1m0_temp, sqrt(sig2_Y))
      Y_z1m1[ind_Uindex] = rnorm(num_Uindex, Y_z1m1_temp, sqrt(sig2_Y))
    }
  } else if (Ydist == "binary"){
    # Generate binary outcome Y
    beta_Y0 = -1.0
    beta_YM =  1.0
    beta_YZ = -beta_Y0*2
    beta_YV = c(-0.5, -0.5,  1.0)
    beta_YC = c( 0.5, -0.5,  0.5, -0.5,  0.5, -0.5)
    
    Y_z0m0 = numeric(N)
    Y_z1m0 = numeric(N)
    Y_z1m1 = numeric(N)
    for (j in 1:J) {
      ind_Uindex = which(Uindex == j)
      num_Uindex = length(ind_Uindex)
      
      C_temp = C[ind_Uindex,]
      V_temp = V[ind_Uindex,]
      M_z0_temp = M_z0[ind_Uindex]
      M_z1_temp = M_z1[ind_Uindex]
      
      Y_z0m0_temp = beta_Y0 + z0 * beta_YZ + C_temp %*% beta_YC + V_temp %*% beta_YV + psi_Y[j]
      Y_z1m0_temp = beta_Y0 + z1 * beta_YZ + C_temp %*% beta_YC + V_temp %*% beta_YV + psi_Y[j]
      Y_z1m1_temp = beta_Y0 + z1 * beta_YZ + C_temp %*% beta_YC + V_temp %*% beta_YV + psi_Y[j]
      if (simple) {
        Y_z0m0_temp = Y_z0m0_temp + M_z0_temp * beta_YM
        Y_z1m0_temp = Y_z1m0_temp + M_z0_temp * beta_YM
        Y_z1m1_temp = Y_z1m1_temp + M_z1_temp * beta_YM
      } else {
        Y_z0m0_temp = Y_z0m0_temp + 0.3 * (M_z0_temp * z0 + M_z0_temp^{2})
        Y_z1m0_temp = Y_z1m0_temp + 0.3 * (M_z0_temp * z1 + M_z0_temp^{2})
        Y_z1m1_temp = Y_z1m1_temp + 0.3 * (M_z1_temp * z1 + M_z1_temp^{2})
      }
      Y_z0m0[ind_Uindex] = rbinom(num_Uindex, 1, pnorm(Y_z0m0_temp))
      Y_z1m0[ind_Uindex] = rbinom(num_Uindex, 1, pnorm(Y_z1m0_temp))
      Y_z1m1[ind_Uindex] = rbinom(num_Uindex, 1, pnorm(Y_z1m1_temp))
    }
  }
  
  # ------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------
  M = ifelse(Z==z0, M_z0, M_z1)
  Y = ifelse(Z==z0, Y_z0m0, Y_z1m1)
  
  E_true = c(mean(Y_z1m1), mean(Y_z1m0), mean(Y_z0m0))
  # E_true = c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])
  
  return(list(Y=Y, M=M, C=C, Z=Z, V=V, Uindex=Uindex, E_true=E_true))
}

# ------------------------------------------------------------------------
# Define posterior summary function
POSTsummary = function(Y, esttype, quantile_alpha, Ycon = NULL){
  Ytmp = Y
  Yquantile = quantile(Ytmp, c(0.25,0.75), na.rm = T)
  Yiqr = diff(Yquantile)
  if(length(Ycon)==0){
    Ycon = 10
  }
  YiqrLOWER = which(Ytmp < (Yquantile[1] - Ycon*Yiqr))
  if (length(YiqrLOWER) > 0) {
    Ytmp = Ytmp[-YiqrLOWER]
  }
  YiqrUPPER = which(Ytmp > (Yquantile[2] + Ycon*Yiqr))
  if (length(YiqrUPPER) > 0) {
    Ytmp = Ytmp[-YiqrUPPER]
  }
  
  # Estimates: mean or median
  if (esttype == "median"){
    Y_est_mc = median(Ytmp, na.rm = T)
  } else {
    Y_est_mc = mean(Ytmp, na.rm = T)
  }
  
  # Calculate median of posterior for Y
  Y_sd_mc = sd(Ytmp, na.rm = T)
  Y_quantile_mc = quantile(Ytmp, quantile_alpha, na.rm = T)
  Y_CIlength_mc = abs(diff(Y_quantile_mc))
  Y_quantile025_mc = min(Y_quantile_mc, na.rm = T)
  Y_quantile975_mc = max(Y_quantile_mc, na.rm = T)
  
  # Save the results
  Y_result_mc = cbind(Y_est_mc, Y_sd_mc, Y_quantile025_mc, Y_quantile975_mc, Y_CIlength_mc)
  colnames(Y_result_mc) = c("estimYs","sd", "quantile025", "quantile975", "CIlength95")
  rownames(Y_result_mc) = c("")
  
  return(Y_result_mc)
}

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
BBmediationPOST = function(object,  # object from rBARTmediation
                           C.test,	# matrix C to predict at
                           V.test,	# matrix V to predict at
                           Uindex,
                           esttype = "mean",
                           saveall = FALSE,
                           CE = TRUE){
  # object = GLMfit
  # # object = BARTfit
  # C.test = C
  # V.test = V
  # Uindex = Uindex
  
  # matrix X to predict at
  X.test = cbind(C.test, V.test)
  
  # Significance Level alpha
  level = 0.05
  
  if (inherits(object,"rBARTmediation")) {
    object0 = object$object0
    object1 = object$object1
    
    # Define POST function for continuous Y and continuous M using BB
    # predict(rBARTmediation)
    # BARTfitPRED = predict(object, X.test, Uindex)
    BARTfitPRED = predict(object0, object1, X.test, Uindex)
    
    # Constants
    N = nrow(X.test)
    J = length(unique(Uindex))
    n_j = as.numeric(table(Uindex))
    n_MCMC = nrow(BARTfitPRED$Yz0m0.test)
    
    # Parameters for individual-level confounder
    a_pi = n_j
    a_rho = rep(1,J)
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = numeric(n_MCMC)
    E_Y_z1m0 = numeric(n_MCMC)
    E_Y_z1m1 = numeric(n_MCMC)
    for (post_reps in 1:n_MCMC) {
      Yhat_z0m0_mc = BARTfitPRED$Yz0m0.test[post_reps,]
      Yhat_z1m0_mc = BARTfitPRED$Yz1m0.test[post_reps,]
      Yhat_z1m1_mc = BARTfitPRED$Yz1m1.test[post_reps,]
      
      # Update culster-level confounders parameters
      # Update individual&culster-level confounders parameters
      piPars = sapply(1:J, function(l) c(rdirichlet_cpp(1,rep(1,a_pi[l]))))
      rhoPars = c(rdirichlet_cpp(1,a_rho))
      rhopiPars = unlist(sapply(1:J, function(l) rhoPars[l]*piPars[[l]]))
      
      E_Y_z0m0_mc = sum(rhopiPars * Yhat_z0m0_mc)
      E_Y_z1m0_mc = sum(rhopiPars * Yhat_z1m0_mc)
      E_Y_z1m1_mc = sum(rhopiPars * Yhat_z1m1_mc)
      
      # E_Y_z0m0_mc = rnorm(1, E_Y_z0m0_mc, object$iYsigest[post_reps]/sqrt(N))
      # E_Y_z1m0_mc = rnorm(1, E_Y_z1m0_mc, object$iYsigest[post_reps]/sqrt(N))
      # E_Y_z1m1_mc = rnorm(1, E_Y_z1m1_mc, object$iYsigest[post_reps]/sqrt(N))
      # rnorm(N, Yhat_z1m1_mc - Yhat_z1m0_mc, object$iYsigest[post_reps])
      # rnorm(N, Yhat_z1m0_mc - Yhat_z0m0_mc, object$iYsigest[post_reps])
      # rnorm(N, Yhat_z1m1_mc - Yhat_z0m0_mc, object$iYsigest[post_reps])
      
      # Causal Effects
      E_Y_z0m0[post_reps] = mean(E_Y_z0m0_mc)
      E_Y_z1m0[post_reps] = mean(E_Y_z1m0_mc)
      E_Y_z1m1[post_reps] = mean(E_Y_z1m1_mc)
      
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  } else if (inherits(object,"rGLMmediation")) {
    # object0 = object$object0
    # object1 = object$object1
    # n_MCMC = object0$constants$n_MCMC
    
    n_MCMC = object$constants$n_MCMC
    
    # 
    z0 = 0
    z1 = 1
    
    # Define POST function for continuous Y and continuous M using BB
    # Constants
    N = nrow(X.test)
    J = length(unique(Uindex))
    n_j = as.numeric(table(Uindex))
    
    matX.test = cbind(1, X.test)
    
    # Parameters for individual-level confounder
    a_pi = n_j
    a_rho = rep(1,J)
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = numeric(n_MCMC)
    E_Y_z1m0 = numeric(n_MCMC)
    E_Y_z1m1 = numeric(n_MCMC)
    for (post_reps in 1:n_MCMC) {
      # ------------------------------------------------------------------------------
      YbetaPars = object$MCMCposteriors$YbetaLists[[post_reps]]
      Ysig2Pars = object$MCMCposteriors$Ysig2Lists[[post_reps]]
      # YpsiPars  = object$MCMCposteriors$YpsiLists[[post_reps]]
      MbetaPars = object$MCMCposteriors$MbetaLists[[post_reps]]
      Msig2Pars = object$MCMCposteriors$Msig2Lists[[post_reps]]
      # MpsiPars  = object$MCMCposteriors$MpsiLists[[post_reps]]
      YMpsiPars = object$MCMCposteriors$YMpsiLists[[post_reps]]
      meanYMpsi = object$MCMCposteriors$meanYMpsiLists[[post_reps]]
      covYMpsiPars = object$MCMCposteriors$covYMpsiLists[[post_reps]]

      # random-effects
      YMpsi_J = rmvn_cpp(J, meanYMpsi, covYMpsiPars) # *sqrt(diag(YMpsiPars)) # rbind(YpsiPars,MpsiPars)
      Ypsi_N = YMpsi_J[1,Uindex] # rep(0, N)
      Mpsi_N = YMpsi_J[2,Uindex] # rep(0, N)

      # Normal model for continuous M
      # draw M0_mc and M1_mc
      matX_z0_mc = cbind(1, z0, X.test)
      matX_z1_mc = cbind(1, z1, X.test)
      Mhat_z0_mc = matX_z0_mc %*% MbetaPars + Mpsi_N
      Mhat_z1_mc = matX_z1_mc %*% MbetaPars + Mpsi_N
      if (object$typeM == "continuous") {
        MsigPars = sqrt(Msig2Pars)
        M0_mc = rnorm(N, Mhat_z0_mc, MsigPars)
        M1_mc = rnorm(N, Mhat_z1_mc, MsigPars)
      } else if (object$typeM == "binary") {
      if (is.null(object$modelM)) {
        callM = plogis
      } else {
        callM = pnorm
      }
        M0_mc = rbinom(N,1, callM(Mhat_z0_mc))
        M1_mc = rbinom(N,1, callM(Mhat_z1_mc))
      }

      # Normal model for continuous Y
      # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      matM_z0m0_mc = cbind(1, M0_mc, z0, X.test)
      matM_z1m0_mc = cbind(1, M0_mc, z1, X.test)
      matM_z1m1_mc = cbind(1, M1_mc, z1, X.test)
      Yhat_z0m0_mc = c(matM_z0m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m0_mc = c(matM_z1m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m1_mc = c(matM_z1m1_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)

      if (object$typeY == "continuous") {
        YsigPars = sqrt(Ysig2Pars)
        Yhat_z0m0_mc = rnorm(N, Yhat_z0m0_mc, YsigPars)
        Yhat_z1m0_mc = rnorm(N, Yhat_z1m0_mc, YsigPars)
        Yhat_z1m1_mc = rnorm(N, Yhat_z1m1_mc, YsigPars)
      } else if (object$typeY == "binary") {
        Yhat_z0m0_mc = rbinom(N, 1, plogis(Yhat_z0m0_mc))
        Yhat_z1m0_mc = rbinom(N, 1, plogis(Yhat_z1m0_mc))
        Yhat_z1m1_mc = rbinom(N, 1, plogis(Yhat_z1m1_mc))
      }
      
      # # -----------------------------------------------------------------------------
      # YbetaPars0 = object0$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars0 = object0$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars0  = object0$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars0 = object0$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars0 = object0$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiPars0  = object0$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars0 = object0$MCMCposteriors$YMpsiLists[[post_reps]]
      # meanYMpsi0 = object0$MCMCposteriors$meanYMpsiLists[[post_reps]]
      # covYMpsiPars0 = object0$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # YbetaPars1 = object1$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars1 = object1$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars1  = object1$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars1 = object1$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars1 = object1$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiPars1  = object1$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars1 = object1$MCMCposteriors$YMpsiLists[[post_reps]]
      # meanYMpsi1 = object1$MCMCposteriors$meanYMpsiLists[[post_reps]]
      # covYMpsiPars1 = object1$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # # random-effects
      # YMpsi_J0 = rmvn_cpp(J, meanYMpsi0, covYMpsiPars0) # * sqrt(diag(YMpsiPars0)) # rbind(YpsiPars0,MpsiPars0)
      # Ypsi_N0 = YMpsi_J0[1,Uindex]
      # Mpsi_N0 = YMpsi_J0[2,Uindex]
      # 
      # YMpsi_J1 = rmvn_cpp(J, meanYMpsi1, covYMpsiPars1) # * sqrt(diag(YMpsiPars1)) # rbind(YpsiPars1,MpsiPars1)
      # Ypsi_N1 = YMpsi_J1[1,Uindex]
      # Mpsi_N1 = YMpsi_J1[2,Uindex]
      # 
      # # Normal model for continuous M
      # # draw M0_mc and M1_mc
      # Mhat_z0_mc = matX.test %*% MbetaPars0 + Mpsi_N0
      # Mhat_z1_mc = matX.test %*% MbetaPars1 + Mpsi_N1
      # 
      # if (object0$typeM == "continuous") {
      #   MsigPars0  = sqrt(Msig2Pars0)
      #   MsigPars1  = sqrt(Msig2Pars1)
      #   M0_mc = rnorm(N, Mhat_z0_mc, MsigPars0)
      #   M1_mc = rnorm(N, Mhat_z1_mc, MsigPars1)
      # } else if (object0$typeM == "binary") {
      #   if (is.null(object0$modelM)) {
      #     callM = plogis
      #   } else {
      #     callM = pnorm
      #   }
      #   M0_mc = rbinom(N,1, callM(Mhat_z0_mc))
      #   M1_mc = rbinom(N,1, callM(Mhat_z1_mc))
      # }
      # 
      # # Normal model for continuous Y
      # # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      # matM_m0_mc = cbind(1, M0_mc, X.test)
      # matM_m1_mc = cbind(1, M1_mc, X.test)
      # Yhat_z0m0_mc = c(matM_m0_mc %*% YbetaPars0)
      # Yhat_z1m0_mc = c(matM_m0_mc %*% YbetaPars1)
      # Yhat_z1m1_mc = c(matM_m1_mc %*% YbetaPars1)
      # 
      # if (object0$typeY == "continuous") {
      #   YsigPars0 = sqrt(Ysig2Pars0)
      #   YsigPars1 = sqrt(Ysig2Pars1)
      #   Yhat_z0m0_mc = rnorm(N, Yhat_z0m0_mc, YsigPars0)
      #   Yhat_z1m0_mc = rnorm(N, Yhat_z1m0_mc, YsigPars1)
      #   Yhat_z1m1_mc = rnorm(N, Yhat_z1m1_mc, YsigPars1)
      # } else if (object0$typeY == "binary") {
      #   Yhat_z0m0_mc = rbinom(N, 1, plogis(Yhat_z0m0_mc))
      #   Yhat_z1m0_mc = rbinom(N, 1, plogis(Yhat_z1m0_mc))
      #   Yhat_z1m1_mc = rbinom(N, 1, plogis(Yhat_z1m1_mc))
      # }
      
      # -----------------------------------------------------------------------------
      # Update culster-level confounders parameters
      # Update individual&culster-level confounders parameters
      piPars = sapply(1:J, function(l) c(rdirichlet_cpp(1,rep(1,a_pi[l]))))
      rhoPars = c(rdirichlet_cpp(1,a_rho))
      rhopiPars = unlist(sapply(1:J, function(l) rhoPars[l]*piPars[[l]]))
      
      E_Y_z0m0_mc = sum(rhopiPars * Yhat_z0m0_mc)
      E_Y_z1m0_mc = sum(rhopiPars * Yhat_z1m0_mc)
      E_Y_z1m1_mc = sum(rhopiPars * Yhat_z1m1_mc)
      
      # Causal Effects
      E_Y_z0m0[post_reps] = mean(E_Y_z0m0_mc)
      E_Y_z1m0[post_reps] = mean(E_Y_z1m0_mc)
      E_Y_z1m1[post_reps] = mean(E_Y_z1m1_mc)
      
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  }
  
  if (CE) {
    NIE = E_Y_z1m1 - E_Y_z1m0
    NDE = E_Y_z1m0 - E_Y_z0m0
    ATE = E_Y_z1m1 - E_Y_z0m0
  } else {
    NIE = E_Y_z1m1
    NDE = E_Y_z1m0
    ATE = E_Y_z0m0
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
  
  # Calculate median of posterior for NIE, NDE, and ATE
  NIE_result_mc = POSTsummary(NIE, "median", quantile_alpha)
  NDE_result_mc = POSTsummary(NDE, "median", quantile_alpha)
  ATE_result_mc = POSTsummary(ATE, "median", quantile_alpha)
  
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
HBBmediationPOST = function(object, # object from predict(rBARTmediation)
                            C.test,	# matrix C to predict at
                            V.test,	# matrix V to predict at
                            Uindex,
                            esttype = "mean",
                            saveall = FALSE,
                            CE = TRUE){
  # object = GLMfit
  # object = BARTfit
  # C.test = C
  # V.test = V
  # Uindex = Uindex
  
  # matrix X to predict at
  X.test = cbind(C.test, V.test)
  
  # Significance Level alpha
  level = 0.05
  
  # Define POST function for continuous Y and continuous M using HBB
  N = nrow(X.test)
  J = length(unique(Uindex))
  n_j = as.numeric(table(Uindex))
  
  matX.test = cbind(1, X.test)
  
  # Parameters for cluster-level confounder
  a_rho = rep(1,J)
  
  # Parameters for individual-level confounder
  a_pi = rep(1,N)
  
  # Parameters for alpha (concentration parameter)
  tau_ome = 1 # HBB
  alpha_ome = (N/n_j)*tau_ome
  
  if (inherits(object,"rBARTmediation")) {
    object0 = object$object0
    object1 = object$object1
    
    # predict(rBARTmediation)
    # BARTfitPRED = predict(object, X.test, Uindex)
    BARTfitPRED = predict(object0, object1, X.test, Uindex)
    
    n_MCMC = nrow(BARTfitPRED$Yz0m0.test)

    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = numeric(n_MCMC)
    E_Y_z1m0 = numeric(n_MCMC)
    E_Y_z1m1 = numeric(n_MCMC)
    for (post_reps in 1:n_MCMC) {
      Yhat_z0m0_mc = BARTfitPRED$Yz0m0.test[post_reps,]
      Yhat_z1m0_mc = BARTfitPRED$Yz1m0.test[post_reps,]
      Yhat_z1m1_mc = BARTfitPRED$Yz1m1.test[post_reps,]
      
      # Update culster-level confounders parameters
      # Update individual-level confounders parameters
      # Update the parameters for individual-level confounders given culster-level confounders
      VrhoPars = c(rdirichlet_cpp(1,a_rho))
      CpiPars = c(rdirichlet_cpp(1,a_pi))
      CpiJPars = sapply(1:J, function(j) rdirichlet_cpp(1,(alpha_ome[j] * CpiPars + (Uindex==j))))
      
      E_Y_z0m0_mc = sum(VrhoPars * apply(CpiJPars * Yhat_z0m0_mc, 2, sum))
      E_Y_z1m0_mc = sum(VrhoPars * apply(CpiJPars * Yhat_z1m0_mc, 2, sum))
      E_Y_z1m1_mc = sum(VrhoPars * apply(CpiJPars * Yhat_z1m1_mc, 2, sum))
      
      # Causal Effects
      E_Y_z0m0[post_reps] = mean(E_Y_z0m0_mc)
      E_Y_z1m0[post_reps] = mean(E_Y_z1m0_mc)
      E_Y_z1m1[post_reps] = mean(E_Y_z1m1_mc)
      
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  } else if (inherits(object,"rGLMmediation")) {
    # object0 = object$object0
    # object1 = object$object1
    # n_MCMC = object0$constants$n_MCMC
    
    n_MCMC = object$constants$n_MCMC
    
    # 
    z0 = 0
    z1 = 1
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = numeric(n_MCMC)
    E_Y_z1m0 = numeric(n_MCMC)
    E_Y_z1m1 = numeric(n_MCMC)
    for (post_reps in 1:n_MCMC) {
      # ------------------------------------------------------------------------------
      YbetaPars = object$MCMCposteriors$YbetaLists[[post_reps]]
      Ysig2Pars = object$MCMCposteriors$Ysig2Lists[[post_reps]]
      # YpsiPars  = object$MCMCposteriors$YpsiLists[[post_reps]]
      MbetaPars = object$MCMCposteriors$MbetaLists[[post_reps]]
      Msig2Pars = object$MCMCposteriors$Msig2Lists[[post_reps]]
      # MpsiPars  = object$MCMCposteriors$MpsiLists[[post_reps]]
      YMpsiPars = object$MCMCposteriors$YMpsiLists[[post_reps]]
      meanYMpsi = object$MCMCposteriors$meanYMpsiLists[[post_reps]]
      covYMpsiPars = object$MCMCposteriors$covYMpsiLists[[post_reps]]

      # random-effects
      YMpsi_J = rmvn_cpp(J, meanYMpsi, covYMpsiPars) # *sqrt(diag(YMpsiPars)) # rbind(YpsiPars,MpsiPars)
      Ypsi_N = YMpsi_J[1,Uindex] # rep(0, N)
      Mpsi_N = YMpsi_J[2,Uindex] # rep(0, N)

      # Normal model for continuous M
      # draw M0_mc and M1_mc
      matX_z0_mc = cbind(1, z0, X.test)
      matX_z1_mc = cbind(1, z1, X.test)
      Mhat_z0_mc = matX_z0_mc %*% MbetaPars + Mpsi_N
      Mhat_z1_mc = matX_z1_mc %*% MbetaPars + Mpsi_N
      if (object$typeM == "continuous") {
        MsigPars = sqrt(Msig2Pars)
        M0_mc = rnorm(N, Mhat_z0_mc, MsigPars)
        M1_mc = rnorm(N, Mhat_z1_mc, MsigPars)
      } else if (object$typeM == "binary") {
      if (is.null(object$modelM)) {
        callM = plogis
      } else {
        callM = pnorm
      }
        M0_mc = rbinom(N,1, callM(Mhat_z0_mc))
        M1_mc = rbinom(N,1, callM(Mhat_z1_mc))
      }

      # Normal model for continuous Y
      # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      matM_z0m0_mc = cbind(1, M0_mc, z0, X.test)
      matM_z1m0_mc = cbind(1, M0_mc, z1, X.test)
      matM_z1m1_mc = cbind(1, M1_mc, z1, X.test)
      Yhat_z0m0_mc = c(matM_z0m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m0_mc = c(matM_z1m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m1_mc = c(matM_z1m1_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)

      if (object$typeY == "continuous") {
        YsigPars = sqrt(Ysig2Pars)
        Yhat_z0m0_mc = rnorm(N, Yhat_z0m0_mc, YsigPars)
        Yhat_z1m0_mc = rnorm(N, Yhat_z1m0_mc, YsigPars)
        Yhat_z1m1_mc = rnorm(N, Yhat_z1m1_mc, YsigPars)
      } else if (object$typeY == "binary") {
        Yhat_z0m0_mc = rbinom(N, 1, plogis(Yhat_z0m0_mc))
        Yhat_z1m0_mc = rbinom(N, 1, plogis(Yhat_z1m0_mc))
        Yhat_z1m1_mc = rbinom(N, 1, plogis(Yhat_z1m1_mc))
      }
      
      # # -----------------------------------------------------------------------------
      # YbetaPars0 = object0$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars0 = object0$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars0  = object0$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars0 = object0$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars0 = object0$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiPars0  = object0$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars0 = object0$MCMCposteriors$YMpsiLists[[post_reps]]
      # meanYMpsi0 = object0$MCMCposteriors$meanYMpsiLists[[post_reps]]
      # covYMpsiPars0 = object0$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # YbetaPars1 = object1$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars1 = object1$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars1  = object1$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars1 = object1$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars1 = object1$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiPars1  = object1$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars1 = object1$MCMCposteriors$YMpsiLists[[post_reps]]
      # meanYMpsi1 = object1$MCMCposteriors$meanYMpsiLists[[post_reps]]
      # covYMpsiPars1 = object1$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # # random-effects
      # YMpsi_J0 = rmvn_cpp(J, meanYMpsi0, covYMpsiPars0) # * sqrt(diag(YMpsiPars0)) # rbind(YpsiPars0,MpsiPars0)
      # Ypsi_N0 = YMpsi_J0[1,Uindex]
      # Mpsi_N0 = YMpsi_J0[2,Uindex]
      # 
      # YMpsi_J1 = rmvn_cpp(J, meanYMpsi1, covYMpsiPars1) # * sqrt(diag(YMpsiPars1)) # rbind(YpsiPars1,MpsiPars1)
      # Ypsi_N1 = YMpsi_J1[1,Uindex]
      # Mpsi_N1 = YMpsi_J1[2,Uindex]
      # 
      # # Normal model for continuous M
      # # draw M0_mc and M1_mc
      # Mhat_z0_mc = matX.test %*% MbetaPars0 + Mpsi_N0
      # Mhat_z1_mc = matX.test %*% MbetaPars1 + Mpsi_N1
      # 
      # if (object0$typeM == "continuous") {
      #   MsigPars0  = sqrt(Msig2Pars0)
      #   MsigPars1  = sqrt(Msig2Pars1)
      #   M0_mc = rnorm(N, Mhat_z0_mc, MsigPars0)
      #   M1_mc = rnorm(N, Mhat_z1_mc, MsigPars1)
      # } else if (object0$typeM == "binary") {
      #   if (is.null(object0$modelM)) {
      #     callM = plogis
      #   } else {
      #     callM = pnorm
      #   }
      #   M0_mc = rbinom(N,1, callM(Mhat_z0_mc))
      #   M1_mc = rbinom(N,1, callM(Mhat_z1_mc))
      # }
      # 
      # # Normal model for continuous Y
      # # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      # matM_m0_mc = cbind(1, M0_mc, X.test)
      # matM_m1_mc = cbind(1, M1_mc, X.test)
      # Yhat_z0m0_mc = c(matM_m0_mc %*% YbetaPars0)
      # Yhat_z1m0_mc = c(matM_m0_mc %*% YbetaPars1)
      # Yhat_z1m1_mc = c(matM_m1_mc %*% YbetaPars1)
      # 
      # if (object0$typeY == "continuous") {
      #   YsigPars0 = sqrt(Ysig2Pars0)
      #   YsigPars1 = sqrt(Ysig2Pars1)
      #   Yhat_z0m0_mc = rnorm(N, Yhat_z0m0_mc, YsigPars0)
      #   Yhat_z1m0_mc = rnorm(N, Yhat_z1m0_mc, YsigPars1)
      #   Yhat_z1m1_mc = rnorm(N, Yhat_z1m1_mc, YsigPars1)
      # } else if (object0$typeY == "binary") {
      #   Yhat_z0m0_mc = rbinom(N, 1, plogis(Yhat_z0m0_mc))
      #   Yhat_z1m0_mc = rbinom(N, 1, plogis(Yhat_z1m0_mc))
      #   Yhat_z1m1_mc = rbinom(N, 1, plogis(Yhat_z1m1_mc))
      # }
      
      # -----------------------------------------------------------------------------
      # Update culster-level confounders parameters
      # Update individual-level confounders parameters
      # Update the parameters for individual-level confounders given culster-level confounders
      VrhoPars = c(rdirichlet_cpp(1,a_rho))
      CpiPars = c(rdirichlet_cpp(1,a_pi))
      CpiJPars = sapply(1:J, function(j) rdirichlet_cpp(1,(alpha_ome[j] * CpiPars + (Uindex==j))))
      
      E_Y_z0m0_mc = sum(VrhoPars * apply(CpiJPars * Yhat_z0m0_mc,2,sum))
      E_Y_z1m0_mc = sum(VrhoPars * apply(CpiJPars * Yhat_z1m0_mc,2,sum))
      E_Y_z1m1_mc = sum(VrhoPars * apply(CpiJPars * Yhat_z1m1_mc,2,sum))
      
      # Causal Effects
      E_Y_z0m0[post_reps] = mean(E_Y_z0m0_mc)
      E_Y_z1m0[post_reps] = mean(E_Y_z1m0_mc)
      E_Y_z1m1[post_reps] = mean(E_Y_z1m1_mc)
      
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  }
  
  if (CE) {
    NIE = E_Y_z1m1 - E_Y_z1m0
    NDE = E_Y_z1m0 - E_Y_z0m0
    ATE = E_Y_z1m1 - E_Y_z0m0
  } else {
    NIE = E_Y_z1m1
    NDE = E_Y_z1m0
    ATE = E_Y_z0m0
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
  
  # Calculate median of posterior for NIE, NDE, and ATE
  NIE_result_mc = POSTsummary(NIE, "median", quantile_alpha)
  NDE_result_mc = POSTsummary(NDE, "median", quantile_alpha)
  ATE_result_mc = POSTsummary(ATE, "median", quantile_alpha)
  
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
TSBBmediationPOST = function(object, # object from predict(rBARTmediation)
                             C.test,	# matrix C to predict at
                             V.test,	# matrix V to predict at
                             Uindex,
                             esttype = "mean",
                             saveall = FALSE, 
                             chi = 1, zeta = 0.5,
                             CE = TRUE){
  # object = GLMfit
  # object = BARTfit
  # C.test = C
  # V.test = V
  # Uindex = Uindex
  
  # a constant to scale the desired rate of decay
  # chi = 1
  # a constant ratio of cluster-level confounders and individual-level confounders
  # zeta = 0.5
  
  # Significance Level alpha
  level = 0.05
  
  N = nrow(C.test)
  J = length(unique(Uindex))
  n_j = as.numeric(table(Uindex))
  p.test = ncol(C.test) + ncol(V.test)
  NJ = N * J
  
  UindexTSBB = numeric(NJ)
  X.testTSBB = matrix(nrow = NJ, ncol = p.test)
  for (l in 1:J) {
    ind_Uindex_l = which(Uindex == l)
    V_lunique = V.test[ind_Uindex_l[1],]
    Xtemp = matrix(nrow = N, ncol = p.test)
    for (j in 1:J) {
      ind_Uindex_j = which(Uindex == j)
      num_Uindex_j = length(ind_Uindex_j)
      C_j = matrix(C.test[ind_Uindex_j,], nrow = num_Uindex_j)
      V_l = matrix(rep(V_lunique,num_Uindex_j), nrow = num_Uindex_j, byrow = T)
      Xtemp[ind_Uindex_j,] = cbind(C_j, V_l)
    }
    # UindexTSBB[((N*(l-1)+1):(N*l))] = l
    X.testTSBB[((N*(l-1)+1):(N*l)),] = Xtemp
  }
  UindexTSBB = rep(1:J, each=N)
  
  matX.testTSBB = cbind(1, X.testTSBB)
  
  # Parameters for cluster-level confounder
  a_rho = rep(1,J)
  
  # Define POST function for continuous Y and binary M using TSBB
  # Distances
  # \lambda_{lj}^{\oemga} = \alpha_{l}^{\oemga} \dist_{lj}: adding additional "pseudosubjects
  tau_ome = 1 # J # 100 # TSBB
  # alpha_ome = rep(tau_ome,J)
  # alpha_ome = (n_j/N) * tau_ome
  # alpha_ome = n_j * tau_ome
  # alpha_ome = (1/n_j) * tau_ome
  alpha_ome = (N/n_j) * tau_ome
  
  # j (row) times l (column)
  Vxi_dist = matrix_dist(scale(V.test),Uindex)
  Cxi_dist = matrix_dist(scale(C.test),Uindex)
  xi_dist = exp(-((zeta)*Vxi_dist + (1-zeta)*Cxi_dist)/chi)
  xi_distPOST = sapply(1:J, function(l) alpha_ome[l]*xi_dist[,l]) + diag(J)
  
  if (inherits(object,"rBARTmediation")) {
    object0 = object$object0
    object1 = object$object1
    
    # predict(rBARTmediation)
    # BARTfitPRED = predict(object, X.testTSBB, UindexTSBB)
    BARTfitPRED = predict(object0, object1, X.testTSBB, UindexTSBB)
    
    # Constants
    n_MCMC = nrow(BARTfitPRED$Yz0m0.test)
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = numeric(n_MCMC)
    E_Y_z1m0 = numeric(n_MCMC)
    E_Y_z1m1 = numeric(n_MCMC)
    for (post_reps in 1:n_MCMC) {
      Yhat_z0m0_mc = BARTfitPRED$Yz0m0.test[post_reps,]
      Yhat_z1m0_mc = BARTfitPRED$Yz1m0.test[post_reps,]
      Yhat_z1m1_mc = BARTfitPRED$Yz1m1.test[post_reps,]
      
      # ------------------------------------------------------------------------------
      # Update parameters for cluster-level confounder (rho: J X 1 vector)
      # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
      # Update parameters for link cluster-individual (omega: J X J matrix)
      rhoPars = c(rdirichlet_cpp(1,rep(1,J)))
      piPars = sapply(1:J, function(l) rdirichlet_cpp(1, rep(1,n_j[l])))
      omePars = sapply(1:J, function(l) rdirichlet_cpp(1, xi_distPOST[,l]))
      
      RhoOmegaPi = numeric(NJ)
      for (l in 1:J) {
        RhoOmegaPitemp = numeric(N)
        for (j in 1:J) {
          ind_Uindex_j = which(Uindex == j)
          RhoOmegaPitemp[ind_Uindex_j] = omePars[j,l]*c(piPars[[j]])
        }
        RhoOmegaPi[((N*(l-1)+1):(N*l))] = rhoPars[l]*RhoOmegaPitemp
      }
      
      E_Y_z0m0_mc = sum(RhoOmegaPi * Yhat_z0m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m0_mc = sum(RhoOmegaPi * Yhat_z1m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m1_mc = sum(RhoOmegaPi * Yhat_z1m1_mc)/sum(RhoOmegaPi)
      
      # Causal Effects
      E_Y_z0m0[post_reps] = mean(E_Y_z0m0_mc)
      E_Y_z1m0[post_reps] = mean(E_Y_z1m0_mc)
      E_Y_z1m1[post_reps] = mean(E_Y_z1m1_mc)
      
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  } else if (inherits(object,"rGLMmediation")) {
    # object0 = object$object0
    # object1 = object$object1
    # n_MCMC = object0$constants$n_MCMC
    
    n_MCMC = object$constants$n_MCMC
    
    # 
    z0 = 0
    z1 = 1
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = numeric(n_MCMC)
    E_Y_z1m0 = numeric(n_MCMC)
    E_Y_z1m1 = numeric(n_MCMC)
    for (post_reps in 1:n_MCMC) {
      # ------------------------------------------------------------------------------
      YbetaPars = object$MCMCposteriors$YbetaLists[[post_reps]]
      Ysig2Pars = object$MCMCposteriors$Ysig2Lists[[post_reps]]
      # YpsiPars  = object$MCMCposteriors$YpsiLists[[post_reps]]
      MbetaPars = object$MCMCposteriors$MbetaLists[[post_reps]]
      Msig2Pars = object$MCMCposteriors$Msig2Lists[[post_reps]]
      # MpsiPars  = object$MCMCposteriors$MpsiLists[[post_reps]]
      YMpsiPars = object$MCMCposteriors$YMpsiLists[[post_reps]]
      meanYMpsi = object$MCMCposteriors$meanYMpsiLists[[post_reps]]
      covYMpsiPars = object$MCMCposteriors$covYMpsiLists[[post_reps]]

      # random-effects
      YMpsi_J = rmvn_cpp(J, meanYMpsi, covYMpsiPars) # *sqrt(diag(YMpsiPars)) # rbind(YpsiPars,MpsiPars)
      Ypsi_N = YMpsi_J[1,UindexTSBB] # rep(0, N)
      Mpsi_N = YMpsi_J[2,UindexTSBB] # rep(0, N)

      # Normal model for continuous M
      # draw M0_mc and M1_mc
      matX_z0_mc = cbind(1, z0, X.testTSBB)
      matX_z1_mc = cbind(1, z1, X.testTSBB)
      Mhat_z0_mc = matX_z0_mc %*% MbetaPars + Mpsi_N
      Mhat_z1_mc = matX_z1_mc %*% MbetaPars + Mpsi_N
      if (object$typeM == "continuous") {
        MsigPars = sqrt(Msig2Pars)
        M0_mc = rnorm(NJ, Mhat_z0_mc, MsigPars)
        M1_mc = rnorm(NJ, Mhat_z1_mc, MsigPars)
      } else if (object$typeM == "binary") {
      if (is.null(object$modelM)) {
        callM = plogis
      } else {
        callM = pnorm
      }
        M0_mc = rbinom(NJ,1, callM(Mhat_z0_mc))
        M1_mc = rbinom(NJ,1, callM(Mhat_z1_mc))
      }

      # Normal model for continuous Y
      # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      matM_z0m0_mc = cbind(1, M0_mc, z0, X.testTSBB)
      matM_z1m0_mc = cbind(1, M0_mc, z1, X.testTSBB)
      matM_z1m1_mc = cbind(1, M1_mc, z1, X.testTSBB)
      Yhat_z0m0_mc = c(matM_z0m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m0_mc = c(matM_z1m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m1_mc = c(matM_z1m1_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)

      if (object$typeY == "continuous") {
        YsigPars = sqrt(Ysig2Pars)
        Yhat_z0m0_mc = rnorm(NJ, Yhat_z0m0_mc, YsigPars)
        Yhat_z1m0_mc = rnorm(NJ, Yhat_z1m0_mc, YsigPars)
        Yhat_z1m1_mc = rnorm(NJ, Yhat_z1m1_mc, YsigPars)
      } else if (object$typeY == "binary") {
        Yhat_z0m0_mc = rbinom(NJ, 1, plogis(Yhat_z0m0_mc))
        Yhat_z1m0_mc = rbinom(NJ, 1, plogis(Yhat_z1m0_mc))
        Yhat_z1m1_mc = rbinom(NJ, 1, plogis(Yhat_z1m1_mc))
      }
      
      # # -----------------------------------------------------------------------------
      # YbetaPars0 = object0$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars0 = object0$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars0  = object0$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars0 = object0$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars0 = object0$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiParsC0  = object0$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars0 = object0$MCMCposteriors$YMpsiLists[[post_reps]]
      # meanYMpsi0 = object0$MCMCposteriors$meanYMpsiLists[[post_reps]]
      # covYMpsiPars0 = object0$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # YbetaPars1 = object1$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars1 = object1$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars1  = object1$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars1 = object1$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars1 = object1$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiPars1  = object1$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars1 = object1$MCMCposteriors$YMpsiLists[[post_reps]]
      # meanYMpsi1 = object1$MCMCposteriors$meanYMpsiLists[[post_reps]]
      # covYMpsiPars1 = object1$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # # random-effects
      # YMpsi_J0 = rmvn_cpp(J, meanYMpsi0, covYMpsiPars0) # * sqrt(diag(YMpsiPars0)) # rbind(YpsiPars0,MpsiPars0)
      # Ypsi_N0 = YMpsi_J0[1,UindexTSBB]
      # Mpsi_N0 = YMpsi_J0[2,UindexTSBB]
      # 
      # YMpsi_J1 = rmvn_cpp(J, meanYMpsi1, covYMpsiPars1) # * sqrt(diag(YMpsiPars1)) # rbind(YpsiPars1,MpsiPars1)
      # Ypsi_N1 = YMpsi_J1[1,UindexTSBB]
      # Mpsi_N1 = YMpsi_J1[2,UindexTSBB]
      # 
      # # Normal model for continuous M
      # # draw M0_mc and M1_mc
      # Mhat_z0_mc = matX.testTSBB %*% MbetaPars0 + Mpsi_N0
      # Mhat_z1_mc = matX.testTSBB %*% MbetaPars1 + Mpsi_N1
      # 
      # if (object0$typeM == "continuous") {
      #   MsigPars0  = sqrt(Msig2Pars0)
      #   MsigPars1  = sqrt(Msig2Pars1)
      #   M0_mc = rnorm(NJ, Mhat_z0_mc, MsigPars0)
      #   M1_mc = rnorm(NJ, Mhat_z1_mc, MsigPars1)
      # } else if (object0$typeM == "binary") {
      #   if (is.null(object0$modelM)) {
      #     callM = plogis
      #   } else {
      #     callM = pnorm
      #   }
      #   M0_mc = rbinom(NJ,1, callM(Mhat_z0_mc))
      #   M1_mc = rbinom(NJ,1, callM(Mhat_z1_mc))
      # }
      # 
      # # Normal model for continuous Y
      # # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      # matM_m0_mc = cbind(1, M0_mc, X.testTSBB)
      # matM_m1_mc = cbind(1, M1_mc, X.testTSBB)
      # Yhat_z0m0_mc = c(matM_m0_mc %*% YbetaPars0)
      # Yhat_z1m0_mc = c(matM_m0_mc %*% YbetaPars1)
      # Yhat_z1m1_mc = c(matM_m1_mc %*% YbetaPars1)
      # 
      # if (object0$typeY == "continuous") {
      #   YsigPars0 = sqrt(Ysig2Pars0)
      #   YsigPars1 = sqrt(Ysig2Pars1)
      #   Yhat_z0m0_mc = rnorm(NJ, Yhat_z0m0_mc, YsigPars0)
      #   Yhat_z1m0_mc = rnorm(NJ, Yhat_z1m0_mc, YsigPars1)
      #   Yhat_z1m1_mc = rnorm(NJ, Yhat_z1m1_mc, YsigPars1)
      # } else if (object0$typeY == "binary") {
      #   Yhat_z0m0_mc = rbinom(NJ, 1, plogis(Yhat_z0m0_mc))
      #   Yhat_z1m0_mc = rbinom(NJ, 1, plogis(Yhat_z1m0_mc))
      #   Yhat_z1m1_mc = rbinom(NJ, 1, plogis(Yhat_z1m1_mc))
      # }
      
      # ------------------------------------------------------------------------------
      # Update parameters for cluster-level confounder (rho: J X 1 vector)
      # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
      # Update parameters for link cluster-individual (omega: J X J matrix)
      rhoPars = c(rdirichlet_cpp(1,rep(1,J)))
      piPars = sapply(1:J, function(l) rdirichlet_cpp(1, rep(1,n_j[l])))
      omePars = sapply(1:J, function(l) rdirichlet_cpp(1, xi_distPOST[,l]))
      
      RhoOmegaPi = numeric(NJ)
      for (l in 1:J) {
        RhoOmegaPitemp = numeric(N)
        for (j in 1:J) {
          ind_Uindex_j = which(Uindex == j)
          RhoOmegaPitemp[ind_Uindex_j] = omePars[j,l]*c(piPars[[j]])
        }
        RhoOmegaPi[((N*(l-1)+1):(N*l))] = rhoPars[l]*RhoOmegaPitemp
      }
      
      E_Y_z0m0_mc = sum(RhoOmegaPi * Yhat_z0m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m0_mc = sum(RhoOmegaPi * Yhat_z1m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m1_mc = sum(RhoOmegaPi * Yhat_z1m1_mc)/sum(RhoOmegaPi)
      
      # Causal Effects
      E_Y_z0m0[post_reps] = mean(E_Y_z0m0_mc)
      E_Y_z1m0[post_reps] = mean(E_Y_z1m0_mc)
      E_Y_z1m1[post_reps] = mean(E_Y_z1m1_mc)
      
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  }
  
  if (CE) {
    NIE = E_Y_z1m1 - E_Y_z1m0
    NDE = E_Y_z1m0 - E_Y_z0m0
    ATE = E_Y_z1m1 - E_Y_z0m0
  } else {
    NIE = E_Y_z1m1
    NDE = E_Y_z1m0
    ATE = E_Y_z0m0
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
  
  # Calculate median of posterior for NIE, NDE, and ATE
  NIE_result_mc = POSTsummary(NIE, "median", quantile_alpha)
  NDE_result_mc = POSTsummary(NDE, "median", quantile_alpha)
  ATE_result_mc = POSTsummary(ATE, "median", quantile_alpha)
  
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

TSBBmediationPOSTsim = function(object, # object from predict(rBARTmediation)
                                C.test,	# matrix C to predict at
                                V.test,	# matrix V to predict at
                                Uindex,
                                esttype = "mean",
                                saveall = FALSE, 
                                list_chi = c(1e-2, 1e-1, 1, 1e1, 1e2), 
                                list_zeta = c(0, 0.5, 1),
                                CE = TRUE){
  # object = BARTfit
  # object = GLMfit
  # C.test = C
  # V.test = V
  # Uindex = Uindex
  # esttype = "mean"
  # saveall = FALSE
  # list_chi = c(1e-2, 1e-1, 1, 1e1, 1e2)
  # list_zeta = c(0, 0.5, 1)
  # CE = TRUE
  
  # Significance Level alpha
  level = 0.05
  
  N = nrow(C.test)
  J = length(unique(Uindex))
  n_j = as.numeric(table(Uindex))
  p.test = ncol(C.test) + ncol(V.test)
  NJ = N * J
  
  UindexTSBB = numeric(NJ)
  X.testTSBB = matrix(nrow = NJ, ncol = p.test)
  for (l in 1:J) {
    ind_Uindex_l = which(Uindex == l)
    V_lunique = V.test[ind_Uindex_l[1],]
    Xtemp = matrix(nrow = N, ncol = p.test)
    for (j in 1:J) {
      ind_Uindex_j = which(Uindex == j)
      num_Uindex_j = length(ind_Uindex_j)
      C_j = matrix(C.test[ind_Uindex_j,], nrow = num_Uindex_j)
      V_l = matrix(rep(V_lunique,num_Uindex_j), nrow = num_Uindex_j, byrow = T)
      Xtemp[ind_Uindex_j,] = cbind(C_j, V_l)
    }
    # UindexTSBB[((N*(l-1)+1):(N*l))] = l
    X.testTSBB[((N*(l-1)+1):(N*l)),] = Xtemp
  }
  UindexTSBB = rep(1:J, each=N)
  
  matX.testTSBB = cbind(1, X.testTSBB)
  
  # Parameters for cluster-level confounder
  a_rho = rep(1,J)
  
  # Define POST function for continuous Y and binary M using TSBB
  # Distances
  # \lambda_{lj}^{\oemga} = \alpha_{l}^{\oemga} \dist_{lj}: adding additional "pseudosubjects
  tau_ome = 1 # J # 100 # TSBB
  # alpha_ome = rep(tau_ome,J)
  # alpha_ome = (n_j/N) * tau_ome
  # alpha_ome = n_j * tau_ome
  # alpha_ome = (1/n_j) * tau_ome
  alpha_ome = (N/n_j) * tau_ome
  
  # j (row) times l (column)
  l_chi = length(list_chi)
  l_zeta = length(list_zeta)
  l_chizeta = l_chi * l_zeta
  
  # ------------------------------------------------------------------------------
  # TSBB
  Vxi_dist = matrix_dist(scale(V.test),Uindex)
  Cxi_dist = matrix_dist(scale(C.test),Uindex)
  
  list_xi_distPOST = list()
  k = 1
  for (chi in list_chi){
    for (zeta in list_zeta){
      xi_dist = exp(-((zeta)*Vxi_dist + (1-zeta)*Cxi_dist)/chi)
      xi_distPOST = sapply(1:J, function(l) alpha_ome[l]*xi_dist[,l]) + diag(J)
      list_xi_distPOST[[k]] = xi_distPOST
      k = k + 1
    }
  }
  
  if (inherits(object,"rBARTmediation")) {
    object0 = object$object0
    object1 = object$object1
    
    # predict(rBARTmediation)
    # BARTfitPRED = predict(object, X.testTSBB, UindexTSBB)
    BARTfitPRED = predict(object0, object1, X.testTSBB, UindexTSBB)
    
    # Constants
    n_MCMC = nrow(BARTfitPRED$Yz0m0.test)
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = matrix(nrow=n_MCMC,ncol=l_chizeta)
    E_Y_z1m0 = matrix(nrow=n_MCMC,ncol=l_chizeta)
    E_Y_z1m1 = matrix(nrow=n_MCMC,ncol=l_chizeta)
    for (post_reps in 1:n_MCMC) {
      Yhat_z0m0_mc = BARTfitPRED$Yz0m0.test[post_reps,]
      Yhat_z1m0_mc = BARTfitPRED$Yz1m0.test[post_reps,]
      Yhat_z1m1_mc = BARTfitPRED$Yz1m1.test[post_reps,]
      
      # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
      # Update parameters for link cluster-individual (omega: J X J matrix)
      # Update parameters for cluster-level confounder (rho: J X 1 vector)
      k = 1
      for (chi in list_chi){
        for (zeta in list_zeta){
          # ------------------------------------------------------------------------------
          # Compputation steps -----------------------------------------------------------
          # ------------------------------------------------------------------------------
          # Update parameters for cluster-level confounder (rho: J X 1 vector)
          # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
          # Update parameters for link cluster-individual (omega: J X J matrix)
          rhoPars = c(rdirichlet_cpp(1,rep(1,J)))
          piPars = sapply(1:J, function(l) rdirichlet_cpp(1, rep(1,n_j[l])))
          omePars = sapply(1:J, function(l) rdirichlet_cpp(1, list_xi_distPOST[[k]][,l]))
          
          RhoOmegaPi = numeric(NJ)
          for (l in 1:J) {
            RhoOmegaPitemp = numeric(N)
            for (j in 1:J) {
              ind_Uindex_j = which(Uindex == j)
              RhoOmegaPitemp[ind_Uindex_j] = omePars[j,l]*c(piPars[[j]])
            }
            RhoOmegaPi[((N*(l-1)+1):(N*l))] = rhoPars[l]*RhoOmegaPitemp
          }
          
          E_Y_z0m0_mc = sum(RhoOmegaPi * Yhat_z0m0_mc)/sum(RhoOmegaPi)
          E_Y_z1m0_mc = sum(RhoOmegaPi * Yhat_z1m0_mc)/sum(RhoOmegaPi)
          E_Y_z1m1_mc = sum(RhoOmegaPi * Yhat_z1m1_mc)/sum(RhoOmegaPi)
          
          # Causal Effects
          E_Y_z0m0[post_reps,k] = mean(E_Y_z0m0_mc)
          E_Y_z1m0[post_reps,k] = mean(E_Y_z1m0_mc)
          E_Y_z1m1[post_reps,k] = mean(E_Y_z1m1_mc)
          
          k = k + 1
        }
      }
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  } else if (inherits(object,"rGLMmediation")) {
    # object0 = object$object0
    # object1 = object$object1
    # n_MCMC = object0$constants$n_MCMC
    
    n_MCMC = object$constants$n_MCMC
    
    # 
    z0 = 0
    z1 = 1
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = matrix(nrow=n_MCMC,ncol=l_chizeta)
    E_Y_z1m0 = matrix(nrow=n_MCMC,ncol=l_chizeta)
    E_Y_z1m1 = matrix(nrow=n_MCMC,ncol=l_chizeta)
    for (post_reps in 1:n_MCMC) {
      # ------------------------------------------------------------------------------
      YbetaPars = object$MCMCposteriors$YbetaLists[[post_reps]]
      Ysig2Pars = object$MCMCposteriors$Ysig2Lists[[post_reps]]
      # YpsiPars  = object$MCMCposteriors$YpsiLists[[post_reps]]
      MbetaPars = object$MCMCposteriors$MbetaLists[[post_reps]]
      Msig2Pars = object$MCMCposteriors$Msig2Lists[[post_reps]]
      # MpsiPars  = object$MCMCposteriors$MpsiLists[[post_reps]]
      YMpsiPars = object$MCMCposteriors$YMpsiLists[[post_reps]]
      meanYMpsi = object$MCMCposteriors$meanYMpsiLists[[post_reps]]
      covYMpsiPars = object$MCMCposteriors$covYMpsiLists[[post_reps]]

      # random-effects
      YMpsi_J = rmvn_cpp(J, meanYMpsi, covYMpsiPars) # *sqrt(diag(YMpsiPars)) # rbind(YpsiPars,MpsiPars)
      Ypsi_N = YMpsi_J[1,UindexTSBB] # rep(0, N)
      Mpsi_N = YMpsi_J[2,UindexTSBB] # rep(0, N)

      # Normal model for continuous M
      # draw M0_mc and M1_mc
      matX_z0_mc = cbind(1, z0, X.testTSBB)
      matX_z1_mc = cbind(1, z1, X.testTSBB)
      Mhat_z0_mc = matX_z0_mc %*% MbetaPars + Mpsi_N
      Mhat_z1_mc = matX_z1_mc %*% MbetaPars + Mpsi_N
      if (object$typeM == "continuous") {
        MsigPars = sqrt(Msig2Pars)
        M0_mc = rnorm(NJ, Mhat_z0_mc, MsigPars)
        M1_mc = rnorm(NJ, Mhat_z1_mc, MsigPars)
      } else if (object$typeM == "binary") {
      if (is.null(object$modelM)) {
        callM = plogis
      } else {
        callM = pnorm
      }
        M0_mc = rbinom(NJ,1, callM(Mhat_z0_mc))
        M1_mc = rbinom(NJ,1, callM(Mhat_z1_mc))
      }

      # Normal model for continuous Y
      # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      matM_z0m0_mc = cbind(1, M0_mc, z0, X.testTSBB)
      matM_z1m0_mc = cbind(1, M0_mc, z1, X.testTSBB)
      matM_z1m1_mc = cbind(1, M1_mc, z1, X.testTSBB)
      Yhat_z0m0_mc = c(matM_z0m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m0_mc = c(matM_z1m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m1_mc = c(matM_z1m1_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)

      if (object$typeY == "continuous") {
        YsigPars = sqrt(Ysig2Pars)
        Yhat_z0m0_mc = rnorm(NJ, Yhat_z0m0_mc, YsigPars)
        Yhat_z1m0_mc = rnorm(NJ, Yhat_z1m0_mc, YsigPars)
        Yhat_z1m1_mc = rnorm(NJ, Yhat_z1m1_mc, YsigPars)
      } else if (object$typeY == "binary") {
        Yhat_z0m0_mc = rbinom(NJ, 1, plogis(Yhat_z0m0_mc))
        Yhat_z1m0_mc = rbinom(NJ, 1, plogis(Yhat_z1m0_mc))
        Yhat_z1m1_mc = rbinom(NJ, 1, plogis(Yhat_z1m1_mc))
      }
      
      # # -----------------------------------------------------------------------------
      # YbetaPars0 = object0$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars0 = object0$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars0  = object0$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars0 = object0$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars0 = object0$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiParsC0  = object0$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars0 = object0$MCMCposteriors$YMpsiLists[[post_reps]]
      # covYMpsiPars0 = object0$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # YbetaPars1 = object1$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars1 = object1$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars1  = object1$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars1 = object1$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars1 = object1$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiPars1  = object1$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars1 = object1$MCMCposteriors$YMpsiLists[[post_reps]]
      # covYMpsiPars1 = object1$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # # random-effects
      # YMpsi_J0 = rmvn_cpp(J, rep(0,2), covYMpsiPars0) # * sqrt(diag(YMpsiPars0)) # rbind(YpsiPars0,MpsiPars0)
      # Ypsi_N0 = YMpsi_J0[1,UindexTSBB]
      # Mpsi_N0 = YMpsi_J0[2,UindexTSBB]
      # 
      # YMpsi_J1 = rmvn_cpp(J, rep(0,2), covYMpsiPars1) # * sqrt(diag(YMpsiPars1)) # rbind(YpsiPars1,MpsiPars1)
      # Ypsi_N1 = YMpsi_J1[1,UindexTSBB]
      # Mpsi_N1 = YMpsi_J1[2,UindexTSBB]
      # 
      # # Normal model for continuous M
      # # draw M0_mc and M1_mc
      # Mhat_z0_mc = matX.testTSBB %*% MbetaPars0 + Mpsi_N0
      # Mhat_z1_mc = matX.testTSBB %*% MbetaPars1 + Mpsi_N1
      # 
      # if (object0$typeM == "continuous") {
      #   MsigPars0  = sqrt(Msig2Pars0)
      #   MsigPars1  = sqrt(Msig2Pars1)
      #   M0_mc = rnorm(NJ, Mhat_z0_mc, MsigPars0)
      #   M1_mc = rnorm(NJ, Mhat_z1_mc, MsigPars1)
      # } else if (object0$typeM == "binary") {
      #   if (is.null(object0$modelM)) {
      #     callM = plogis
      #   } else {
      #     callM = pnorm
      #   }
      #   M0_mc = rbinom(NJ,1, callM(Mhat_z0_mc))
      #   M1_mc = rbinom(NJ,1, callM(Mhat_z1_mc))
      # }
      # 
      # # Normal model for continuous Y
      # # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      # matM_m0_mc = cbind(1, M0_mc, X.testTSBB)
      # matM_m1_mc = cbind(1, M1_mc, X.testTSBB)
      # Yhat_z0m0_mc = c(matM_m0_mc %*% YbetaPars0)
      # Yhat_z1m0_mc = c(matM_m0_mc %*% YbetaPars1)
      # Yhat_z1m1_mc = c(matM_m1_mc %*% YbetaPars1)
      # 
      # if (object0$typeY == "continuous") {
      #   YsigPars0 = sqrt(Ysig2Pars0)
      #   YsigPars1 = sqrt(Ysig2Pars1)
      #   Yhat_z0m0_mc = rnorm(NJ, Yhat_z0m0_mc, YsigPars0)
      #   Yhat_z1m0_mc = rnorm(NJ, Yhat_z1m0_mc, YsigPars1)
      #   Yhat_z1m1_mc = rnorm(NJ, Yhat_z1m1_mc, YsigPars1)
      # } else if (object0$typeY == "binary") {
      #   Yhat_z0m0_mc = rbinom(NJ, 1, plogis(Yhat_z0m0_mc))
      #   Yhat_z1m0_mc = rbinom(NJ, 1, plogis(Yhat_z1m0_mc))
      #   Yhat_z1m1_mc = rbinom(NJ, 1, plogis(Yhat_z1m1_mc))
      # }
      
      # ------------------------------------------------------------------------------
      # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
      # Update parameters for link cluster-individual (omega: J X J matrix)
      # Update parameters for cluster-level confounder (rho: J X 1 vector)
      k = 1
      for (chi in list_chi){
        for (zeta in list_zeta){
          # ------------------------------------------------------------------------------
          # Compputation steps -----------------------------------------------------------
          # ------------------------------------------------------------------------------
          # Update parameters for cluster-level confounder (rho: J X 1 vector)
          # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
          # Update parameters for link cluster-individual (omega: J X J matrix)
          rhoPars = c(rdirichlet_cpp(1,rep(1,J)))
          piPars = sapply(1:J, function(l) rdirichlet_cpp(1, rep(1,n_j[l])))
          omePars = sapply(1:J, function(l) rdirichlet_cpp(1, list_xi_distPOST[[k]][,l]))
          
          RhoOmegaPi = numeric(NJ)
          for (l in 1:J) {
            RhoOmegaPitemp = numeric(N)
            for (j in 1:J) {
              ind_Uindex_j = which(Uindex == j)
              RhoOmegaPitemp[ind_Uindex_j] = omePars[j,l]*c(piPars[[j]])
            }
            RhoOmegaPi[((N*(l-1)+1):(N*l))] = rhoPars[l]*RhoOmegaPitemp
          }
          
          E_Y_z0m0_mc = sum(RhoOmegaPi * Yhat_z0m0_mc)/sum(RhoOmegaPi)
          E_Y_z1m0_mc = sum(RhoOmegaPi * Yhat_z1m0_mc)/sum(RhoOmegaPi)
          E_Y_z1m1_mc = sum(RhoOmegaPi * Yhat_z1m1_mc)/sum(RhoOmegaPi)
          
          # Causal Effects
          E_Y_z0m0[post_reps,k] = mean(E_Y_z0m0_mc)
          E_Y_z1m0[post_reps,k] = mean(E_Y_z1m0_mc)
          E_Y_z1m1[post_reps,k] = mean(E_Y_z1m1_mc)
          
          k = k + 1
        }
      }
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  }
  
  if (CE) {
    NIE = E_Y_z1m1 - E_Y_z1m0
    NDE = E_Y_z1m0 - E_Y_z0m0
    ATE = E_Y_z1m1 - E_Y_z0m0
  } else {
    NIE = E_Y_z1m1
    NDE = E_Y_z1m0
    ATE = E_Y_z0m0
  }
  
  # Estimates: mean or median
  NIE_est_mc = apply(NIE, 2, esttype)
  NDE_est_mc = apply(NDE, 2, esttype)
  ATE_est_mc = apply(ATE, 2, esttype)
  
  quantile_alpha = c(level/2,1-level/2)
  
  # Calculate median of posterior for NIE
  NIE_sd_mc = apply(NIE, 2, sd)
  NIE_quantile_mc = apply(NIE, 2, quantile, quantile_alpha)
  NIE_CIlength_mc = apply(NIE_quantile_mc, 2, diff)
  NIE_quantile025_mc = apply(NIE_quantile_mc, 2, min)
  NIE_quantile975_mc = apply(NIE_quantile_mc, 2, max)
  
  # Calculate median of posterior for NDE
  NDE_sd_mc = apply(NDE, 2, sd)
  NDE_quantile_mc = apply(NDE, 2, quantile, quantile_alpha)
  NDE_CIlength_mc = apply(NDE_quantile_mc, 2, diff)
  NDE_quantile025_mc = apply(NDE_quantile_mc, 2, min)
  NDE_quantile975_mc = apply(NDE_quantile_mc, 2, max)
  
  # Calculate median of posterior for ATE
  ATE_sd_mc = apply(ATE, 2, sd)
  ATE_quantile_mc = apply(ATE, 2, quantile, quantile_alpha)
  ATE_CIlength_mc = apply(ATE_quantile_mc, 2, diff)
  ATE_quantile025_mc = apply(ATE_quantile_mc, 2, min)
  ATE_quantile975_mc = apply(ATE_quantile_mc, 2, max)
  
  # Save the results
  NIE_result_mc = cbind(NIE_est_mc, NIE_sd_mc, NIE_quantile025_mc, NIE_quantile975_mc, NIE_CIlength_mc)
  colnames(NIE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(NIE_result_mc) = paste0(rep(paste0("c",list_chi),each=l_zeta),rep(paste0("z",list_zeta),l_chi))
  
  NDE_result_mc = cbind(NDE_est_mc, NDE_sd_mc, NDE_quantile025_mc, NDE_quantile975_mc, NDE_CIlength_mc)
  colnames(NDE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(NDE_result_mc) = paste0(rep(paste0("c",list_chi),each=l_zeta),rep(paste0("z",list_zeta),l_chi))
  
  ATE_result_mc = cbind(ATE_est_mc, ATE_sd_mc, ATE_quantile025_mc, ATE_quantile975_mc, ATE_CIlength_mc)
  colnames(ATE_result_mc) = c("estimates","sd", "quantile025", "quantile975", "CIlength95")
  rownames(ATE_result_mc) = paste0(rep(paste0("c",list_chi),each=l_zeta),rep(paste0("z",list_zeta),l_chi))
  
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

CcondTSBBmediationPOST = function(object, # object from rBARTmediation
                                  C.test,	# matrix C to predict at
                                  V.test,	# matrix V to predict at
                                  cSindex, # C.cond.index (conditional)
                                  Uindex,
                                  esttype = "mean",
                                  saveall = FALSE,
                                  chi = 1, zeta = 0.5,
                                  CE = TRUE){
  # object = GLMfit
  # object = BARTfit
  # object = MODELfit
  # C.test = C
  # V.test = V
  # cSindex = S0C5
  # Uindex = Uindex
  
  # a constant to scale the desired rate of decay
  # chi = 1
  # a constant ratio of cluster-level confounders and individual-level confounders
  # zeta = 0.5
  
  # Significance Level alpha
  level = 0.05
  
  N = nrow(C.test)
  J = length(unique(Uindex))
  n_j = as.numeric(table(Uindex))
  p.test = ncol(C.test) + ncol(V.test)
  NJ = N * J
  
  UindexTSBB = numeric(NJ)
  X.testTSBB = matrix(nrow = NJ, ncol = p.test)
  for (l in 1:J) {
    ind_Uindex_l = which(Uindex == l)
    V_lunique = V.test[ind_Uindex_l[1],]
    Xtemp = matrix(nrow = N, ncol = p.test)
    for (j in 1:J) {
      ind_Uindex_j = which(Uindex == j)
      num_Uindex_j = length(ind_Uindex_j)
      C_j = matrix(C.test[ind_Uindex_j,], nrow = num_Uindex_j)
      V_l = matrix(rep(V_lunique,num_Uindex_j), nrow = num_Uindex_j, byrow = T)
      Xtemp[ind_Uindex_j,] = cbind(C_j, V_l)
    }
    # UindexTSBB[((N*(l-1)+1):(N*l))] = l
    X.testTSBB[((N*(l-1)+1):(N*l)),] = Xtemp
  }
  UindexTSBB = rep(1:J, each=N)
  cSindexTSBB = rep(cSindex, J)
  
  matX.testTSBB = cbind(1, X.testTSBB)
  
  # Parameters for cluster-level confounder
  a_rho = rep(1,J)
  
  # Define POST function for continuous Y and binary M using TSBB
  # Distances
  # \lambda_{lj}^{\oemga} = \alpha_{l}^{\oemga} \dist_{lj}: adding additional "pseudosubjects
  tau_ome = 1 # J # 100 # TSBB
  # alpha_ome = rep(tau_ome,J)
  # alpha_ome = (n_j/N) * tau_ome
  # alpha_ome = n_j * tau_ome
  # alpha_ome = (1/n_j) * tau_ome
  alpha_ome = (N/n_j) * tau_ome
  
  # j (row) times l (column)
  Vxi_dist = matrix_dist(scale(V.test),Uindex)
  Cxi_dist = matrix_dist(scale(C.test),Uindex)
  xi_dist = exp(-((zeta)*Vxi_dist + (1-zeta)*Cxi_dist)/chi)
  xi_distPOST = sapply(1:J, function(l) alpha_ome[l]*xi_dist[,l]) + diag(J)
  
  if (inherits(object,"rBARTmediation")) {
    object0 = object$object0
    object1 = object$object1
    
    # predict(rBARTmediation)
    # BARTfitPRED = predict(object, X.testTSBB, UindexTSBB)
    BARTfitPRED = predict(object0, object1, X.testTSBB, UindexTSBB)
    
    # Constants
    n_MCMC = nrow(BARTfitPRED$Yz0m0.test)
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = numeric(n_MCMC)
    E_Y_z1m0 = numeric(n_MCMC)
    E_Y_z1m1 = numeric(n_MCMC)
    for (post_reps in 1:n_MCMC) {
      Yhat_z0m0_mc = BARTfitPRED$Yz0m0.test[post_reps,]
      Yhat_z1m0_mc = BARTfitPRED$Yz1m0.test[post_reps,]
      Yhat_z1m1_mc = BARTfitPRED$Yz1m1.test[post_reps,]
      
      # ------------------------------------------------------------------------------
      # Update parameters for cluster-level confounder (rho: J X 1 vector)
      # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
      # Update parameters for link cluster-individual (omega: J X J matrix)
      rhoPars = c(rdirichlet_cpp(1,rep(1,J)))
      piPars = sapply(1:J, function(l) rdirichlet_cpp(1, rep(1,n_j[l])))
      omePars = sapply(1:J, function(l) rdirichlet_cpp(1, xi_distPOST[,l]))
      
      RhoOmegaPi = numeric(NJ)
      for (l in 1:J) {
        RhoOmegaPitemp = numeric(N)
        for (j in 1:J) {
          ind_Uindex_j = which(Uindex == j)
          RhoOmegaPitemp[ind_Uindex_j] = omePars[j,l]*c(piPars[[j]])
        }
        RhoOmegaPi[((N*(l-1)+1):(N*l))] = rhoPars[l]*RhoOmegaPitemp
      }
      RhoOmegaPi = RhoOmegaPi * cSindexTSBB
      
      E_Y_z0m0_mc = sum(RhoOmegaPi *  Yhat_z0m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m0_mc = sum(RhoOmegaPi *  Yhat_z1m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m1_mc = sum(RhoOmegaPi *  Yhat_z1m1_mc)/sum(RhoOmegaPi)
      
      # Causal Effects
      E_Y_z0m0[post_reps] = mean(E_Y_z0m0_mc)
      E_Y_z1m0[post_reps] = mean(E_Y_z1m0_mc)
      E_Y_z1m1[post_reps] = mean(E_Y_z1m1_mc)
      
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  } else if (inherits(object,"rGLMmediation")) {
    # object0 = object$object0
    # object1 = object$object1
    # n_MCMC = object0$constants$n_MCMC
    
    n_MCMC = object$constants$n_MCMC
    
    # 
    z0 = 0
    z1 = 1
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = numeric(n_MCMC)
    E_Y_z1m0 = numeric(n_MCMC)
    E_Y_z1m1 = numeric(n_MCMC)
    for (post_reps in 1:n_MCMC) {
      # ------------------------------------------------------------------------------
      YbetaPars = object$MCMCposteriors$YbetaLists[[post_reps]]
      Ysig2Pars = object$MCMCposteriors$Ysig2Lists[[post_reps]]
      # YpsiPars  = object$MCMCposteriors$YpsiLists[[post_reps]]
      MbetaPars = object$MCMCposteriors$MbetaLists[[post_reps]]
      Msig2Pars = object$MCMCposteriors$Msig2Lists[[post_reps]]
      # MpsiPars  = object$MCMCposteriors$MpsiLists[[post_reps]]
      YMpsiPars = object$MCMCposteriors$YMpsiLists[[post_reps]]
      meanYMpsi = object$MCMCposteriors$meanYMpsiLists[[post_reps]]
      covYMpsiPars = object$MCMCposteriors$covYMpsiLists[[post_reps]]

      # random-effects
      YMpsi_J = rmvn_cpp(J, meanYMpsi, covYMpsiPars) # *sqrt(diag(YMpsiPars)) # rbind(YpsiPars,MpsiPars)
      Ypsi_N = YMpsi_J[1,UindexTSBB] # rep(0, N)
      Mpsi_N = YMpsi_J[2,UindexTSBB] # rep(0, N)

      # Normal model for continuous M
      # draw M0_mc and M1_mc
      matX_z0_mc = cbind(1, z0, X.testTSBB)
      matX_z1_mc = cbind(1, z1, X.testTSBB)
      Mhat_z0_mc = matX_z0_mc %*% MbetaPars + Mpsi_N
      Mhat_z1_mc = matX_z1_mc %*% MbetaPars + Mpsi_N
      if (object$typeM == "continuous") {
        MsigPars = sqrt(Msig2Pars)
        M0_mc = rnorm(NJ, Mhat_z0_mc, MsigPars)
        M1_mc = rnorm(NJ, Mhat_z1_mc, MsigPars)
      } else if (object$typeM == "binary") {
      if (is.null(object$modelM)) {
        callM = plogis
      } else {
        callM = pnorm
      }
        M0_mc = rbinom(NJ,1, callM(Mhat_z0_mc))
        M1_mc = rbinom(NJ,1, callM(Mhat_z1_mc))
      }

      # Normal model for continuous Y
      # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      matM_z0m0_mc = cbind(1, M0_mc, z0, X.testTSBB)
      matM_z1m0_mc = cbind(1, M0_mc, z1, X.testTSBB)
      matM_z1m1_mc = cbind(1, M1_mc, z1, X.testTSBB)
      Yhat_z0m0_mc = c(matM_z0m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m0_mc = c(matM_z1m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m1_mc = c(matM_z1m1_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)

      if (object$typeY == "continuous") {
        YsigPars = sqrt(Ysig2Pars)
        Yhat_z0m0_mc = rnorm(NJ, Yhat_z0m0_mc, YsigPars)
        Yhat_z1m0_mc = rnorm(NJ, Yhat_z1m0_mc, YsigPars)
        Yhat_z1m1_mc = rnorm(NJ, Yhat_z1m1_mc, YsigPars)
      } else if (object$typeY == "binary") {
        Yhat_z0m0_mc = rbinom(NJ, 1, plogis(Yhat_z0m0_mc))
        Yhat_z1m0_mc = rbinom(NJ, 1, plogis(Yhat_z1m0_mc))
        Yhat_z1m1_mc = rbinom(NJ, 1, plogis(Yhat_z1m1_mc))
      }
      
      # # -----------------------------------------------------------------------------
      # YbetaPars0 = object0$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars0 = object0$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars0  = object0$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars0 = object0$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars0 = object0$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiParsC0  = object0$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars0 = object0$MCMCposteriors$YMpsiLists[[post_reps]]
      # meanYMpsi0 = object0$MCMCposteriors$meanYMpsiLists[[post_reps]]
      # covYMpsiPars0 = object0$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # YbetaPars1 = object1$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars1 = object1$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars1  = object1$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars1 = object1$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars1 = object1$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiPars1  = object1$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars1 = object1$MCMCposteriors$YMpsiLists[[post_reps]]
      # meanYMpsi1 = object1$MCMCposteriors$meanYMpsiLists[[post_reps]]
      # covYMpsiPars1 = object1$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # # random-effects
      # YMpsi_J0 = rmvn_cpp(J, meanYMpsi0, covYMpsiPars0) # * sqrt(diag(YMpsiPars0)) # rbind(YpsiPars0,MpsiPars0)
      # Ypsi_N0 = YMpsi_J0[1,UindexTSBB]
      # Mpsi_N0 = YMpsi_J0[2,UindexTSBB]
      # 
      # YMpsi_J1 = rmvn_cpp(J, meanYMpsi1, covYMpsiPars1) # * sqrt(diag(YMpsiPars1)) # rbind(YpsiPars1,MpsiPars1)
      # Ypsi_N1 = YMpsi_J1[1,UindexTSBB]
      # Mpsi_N1 = YMpsi_J1[2,UindexTSBB]
      # 
      # # Normal model for continuous M
      # # draw M0_mc and M1_mc
      # Mhat_z0_mc = matX.testTSBB %*% MbetaPars0 + Mpsi_N0
      # Mhat_z1_mc = matX.testTSBB %*% MbetaPars1 + Mpsi_N1
      # 
      # if (object0$typeM == "continuous") {
      #   MsigPars0  = sqrt(Msig2Pars0)
      #   MsigPars1  = sqrt(Msig2Pars1)
      #   M0_mc = rnorm(NJ, Mhat_z0_mc, MsigPars0)
      #   M1_mc = rnorm(NJ, Mhat_z1_mc, MsigPars1)
      # } else if (object0$typeM == "binary") {
      #   if (is.null(object0$modelM)) {
      #     callM = plogis
      #   } else {
      #     callM = pnorm
      #   }
      #   M0_mc = rbinom(NJ,1, callM(Mhat_z0_mc))
      #   M1_mc = rbinom(NJ,1, callM(Mhat_z1_mc))
      # }
      # 
      # # Normal model for continuous Y
      # # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      # matM_m0_mc = cbind(1, M0_mc, X.testTSBB)
      # matM_m1_mc = cbind(1, M1_mc, X.testTSBB)
      # Yhat_z0m0_mc = c(matM_m0_mc %*% YbetaPars0)
      # Yhat_z1m0_mc = c(matM_m0_mc %*% YbetaPars1)
      # Yhat_z1m1_mc = c(matM_m1_mc %*% YbetaPars1)
      # 
      # if (object0$typeY == "continuous") {
      #   YsigPars0 = sqrt(Ysig2Pars0)
      #   YsigPars1 = sqrt(Ysig2Pars1)
      #   Yhat_z0m0_mc = rnorm(NJ, Yhat_z0m0_mc, YsigPars0)
      #   Yhat_z1m0_mc = rnorm(NJ, Yhat_z1m0_mc, YsigPars1)
      #   Yhat_z1m1_mc = rnorm(NJ, Yhat_z1m1_mc, YsigPars1)
      # } else if (object0$typeY == "binary") {
      #   Yhat_z0m0_mc = rbinom(NJ, 1, plogis(Yhat_z0m0_mc))
      #   Yhat_z1m0_mc = rbinom(NJ, 1, plogis(Yhat_z1m0_mc))
      #   Yhat_z1m1_mc = rbinom(NJ, 1, plogis(Yhat_z1m1_mc))
      # }
      
      # ------------------------------------------------------------------------------
      # Update parameters for cluster-level confounder (rho: J X 1 vector)
      # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
      # Update parameters for link cluster-individual (omega: J X J matrix)
      rhoPars = c(rdirichlet_cpp(1,rep(1,J)))
      piPars = sapply(1:J, function(l) rdirichlet_cpp(1, rep(1,n_j[l])))
      omePars = sapply(1:J, function(l) rdirichlet_cpp(1, xi_distPOST[,l]))
      
      RhoOmegaPi = numeric(NJ)
      for (l in 1:J) {
        RhoOmegaPitemp = numeric(N)
        for (j in 1:J) {
          ind_Uindex_j = which(Uindex == j)
          RhoOmegaPitemp[ind_Uindex_j] = omePars[j,l]*c(piPars[[j]])
        }
        RhoOmegaPi[((N*(l-1)+1):(N*l))] = rhoPars[l]*RhoOmegaPitemp
      }
      RhoOmegaPi = RhoOmegaPi * cSindexTSBB
      
      E_Y_z0m0_mc = sum(RhoOmegaPi * Yhat_z0m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m0_mc = sum(RhoOmegaPi * Yhat_z1m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m1_mc = sum(RhoOmegaPi * Yhat_z1m1_mc)/sum(RhoOmegaPi)
      
      # Causal Effects
      E_Y_z0m0[post_reps] = mean(E_Y_z0m0_mc)
      E_Y_z1m0[post_reps] = mean(E_Y_z1m0_mc)
      E_Y_z1m1[post_reps] = mean(E_Y_z1m1_mc)
      
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  }
  if (CE) {
    NIE = E_Y_z1m1 - E_Y_z1m0
    NDE = E_Y_z1m0 - E_Y_z0m0
    ATE = E_Y_z1m1 - E_Y_z0m0
  } else {
    NIE = E_Y_z1m1
    NDE = E_Y_z1m0
    ATE = E_Y_z0m0
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
  
  # Calculate median of posterior for NIE, NDE, and ATE
  NIE_result_mc = POSTsummary(NIE, "median", quantile_alpha)
  NDE_result_mc = POSTsummary(NDE, "median", quantile_alpha)
  ATE_result_mc = POSTsummary(ATE, "median", quantile_alpha)
  
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

VcondTSBBmediationPOST = function(object, # object from rBARTmediation
                                  C.test,	# matrix C to predict at
                                  V.test,	# matrix V to predict at
                                  vSindex, # V.cond.index (conditional)
                                  Uindex,
                                  esttype = "mean",
                                  saveall = FALSE,
                                  chi = 1, zeta = 0.5,
                                  CE = TRUE){
  # object = GLMfit
  # object = BARTfit
  # object = MODELfit
  # C.test = C
  # V.test = V
  # vSindex = S_V3_0
  # Uindex = Uindex
  
  # a constant to scale the desired rate of decay
  # chi = 1
  # a constant ratio of cluster-level confounders and individual-level confounders
  # zeta = 0.5
  
  # Significance Level alpha
  level = 0.05
  
  N = nrow(C.test)
  J = length(unique(Uindex))
  n_j = as.numeric(table(Uindex))
  p.test = ncol(C.test) + ncol(V.test)
  NJ = N * J
  
  UindexTSBB = numeric(NJ)
  X.testTSBB = matrix(nrow = NJ, ncol = p.test)
  for (l in 1:J) {
    ind_Uindex_l = which(Uindex == l)
    V_lunique = V.test[ind_Uindex_l[1],]
    Xtemp = matrix(nrow = N, ncol = p.test)
    for (j in 1:J) {
      ind_Uindex_j = which(Uindex == j)
      num_Uindex_j = length(ind_Uindex_j)
      C_j = matrix(C.test[ind_Uindex_j,], nrow = num_Uindex_j)
      V_l = matrix(rep(V_lunique,num_Uindex_j), nrow = num_Uindex_j, byrow = T)
      Xtemp[ind_Uindex_j,] = cbind(C_j, V_l)
    }
    # UindexTSBB[((N*(l-1)+1):(N*l))] = l
    X.testTSBB[((N*(l-1)+1):(N*l)),] = Xtemp
  }
  UindexTSBB = rep(1:J, each=N)
  vSindexTSBB = rep(vSindex, each=J)
  
  matX.testTSBB = cbind(1, X.testTSBB)
  
  # Parameters for cluster-level confounder
  a_rho = rep(1,J)
  
  # Define POST function for continuous Y and binary M using TSBB
  # Distances
  # \lambda_{lj}^{\oemga} = \alpha_{l}^{\oemga} \dist_{lj}: adding additional "pseudosubjects
  tau_ome = 1 # J # 100 # TSBB
  # alpha_ome = rep(tau_ome,J)
  # alpha_ome = (n_j/N) * tau_ome
  # alpha_ome = n_j * tau_ome
  # alpha_ome = (1/n_j) * tau_ome
  alpha_ome = (N/n_j) * tau_ome

  # j (row) times l (column)
  Vxi_dist = matrix_dist(scale(V.test),Uindex)
  Cxi_dist = matrix_dist(scale(C.test),Uindex)
  xi_dist = exp(-((zeta)*Vxi_dist + (1-zeta)*Cxi_dist)/chi)
  xi_distPOST = sapply(1:J, function(l) alpha_ome[l]*xi_dist[,l]) + diag(J)
  
  if (inherits(object,"rBARTmediation")) {
    object0 = object$object0
    object1 = object$object1
    
    # predict(rBARTmediation)
    # BARTfitPRED = predict(object, X.testTSBB, UindexTSBB)
    BARTfitPRED = predict(object0, object1, X.testTSBB, UindexTSBB)
    
    # Constants
    n_MCMC = nrow(BARTfitPRED$Yz0m0.test)
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = numeric(n_MCMC)
    E_Y_z1m0 = numeric(n_MCMC)
    E_Y_z1m1 = numeric(n_MCMC)
    for (post_reps in 1:n_MCMC) {
      Yhat_z0m0_mc = BARTfitPRED$Yz0m0.test[post_reps,]
      Yhat_z1m0_mc = BARTfitPRED$Yz1m0.test[post_reps,]
      Yhat_z1m1_mc = BARTfitPRED$Yz1m1.test[post_reps,]
      
      # ------------------------------------------------------------------------------
      # Update parameters for cluster-level confounder (rho: J X 1 vector)
      # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
      # Update parameters for link cluster-individual (omega: J X J matrix)
      rhoPars = c(rdirichlet_cpp(1,rep(1,J)))
      piPars = sapply(1:J, function(l) rdirichlet_cpp(1, rep(1,n_j[l])))
      omePars = sapply(1:J, function(l) rdirichlet_cpp(1, xi_distPOST[,l]))
      
      RhoOmegaPi = numeric(NJ)
      for (l in 1:J) {
        RhoOmegaPitemp = numeric(N)
        for (j in 1:J) {
          ind_Uindex_j = which(Uindex == j)
          RhoOmegaPitemp[ind_Uindex_j] = omePars[j,l]*c(piPars[[j]])
        }
        RhoOmegaPi[((N*(l-1)+1):(N*l))] = rhoPars[l]*RhoOmegaPitemp
      }
      RhoOmegaPi = RhoOmegaPi * vSindexTSBB
      
      E_Y_z0m0_mc = sum(RhoOmegaPi *  Yhat_z0m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m0_mc = sum(RhoOmegaPi *  Yhat_z1m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m1_mc = sum(RhoOmegaPi *  Yhat_z1m1_mc)/sum(RhoOmegaPi)
      
      # Causal Effects
      E_Y_z0m0[post_reps] = mean(E_Y_z0m0_mc)
      E_Y_z1m0[post_reps] = mean(E_Y_z1m0_mc)
      E_Y_z1m1[post_reps] = mean(E_Y_z1m1_mc)
      
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  } else if (inherits(object,"rGLMmediation")) {
    # object0 = object$object0
    # object1 = object$object1
    # n_MCMC = object0$constants$n_MCMC
    
    n_MCMC = object$constants$n_MCMC
    
    # 
    z0 = 0
    z1 = 1
    
    # Define interation check
    iter_check = floor(n_MCMC/10)
    
    # Causal Effects
    E_Y_z0m0 = numeric(n_MCMC)
    E_Y_z1m0 = numeric(n_MCMC)
    E_Y_z1m1 = numeric(n_MCMC)
    for (post_reps in 1:n_MCMC) {
      # ------------------------------------------------------------------------------
      YbetaPars = object$MCMCposteriors$YbetaLists[[post_reps]]
      Ysig2Pars = object$MCMCposteriors$Ysig2Lists[[post_reps]]
      # YpsiPars  = object$MCMCposteriors$YpsiLists[[post_reps]]
      MbetaPars = object$MCMCposteriors$MbetaLists[[post_reps]]
      Msig2Pars = object$MCMCposteriors$Msig2Lists[[post_reps]]
      # MpsiPars  = object$MCMCposteriors$MpsiLists[[post_reps]]
      YMpsiPars = object$MCMCposteriors$YMpsiLists[[post_reps]]
      meanYMpsi = object$MCMCposteriors$meanYMpsiLists[[post_reps]]
      covYMpsiPars = object$MCMCposteriors$covYMpsiLists[[post_reps]]

      # random-effects
      YMpsi_J = rmvn_cpp(J, meanYMpsi, covYMpsiPars) # *sqrt(diag(YMpsiPars)) # rbind(YpsiPars,MpsiPars)
      Ypsi_N = YMpsi_J[1,UindexTSBB] # rep(0, N)
      Mpsi_N = YMpsi_J[2,UindexTSBB] # rep(0, N)

      # Normal model for continuous M
      # draw M0_mc and M1_mc
      matX_z0_mc = cbind(1, z0, X.testTSBB)
      matX_z1_mc = cbind(1, z1, X.testTSBB)
      Mhat_z0_mc = matX_z0_mc %*% MbetaPars + Mpsi_N
      Mhat_z1_mc = matX_z1_mc %*% MbetaPars + Mpsi_N
      if (object$typeM == "continuous") {
        MsigPars = sqrt(Msig2Pars)
        M0_mc = rnorm(NJ, Mhat_z0_mc, MsigPars)
        M1_mc = rnorm(NJ, Mhat_z1_mc, MsigPars)
      } else if (object$typeM == "binary") {
      if (is.null(object$modelM)) {
        callM = plogis
      } else {
        callM = pnorm
      }
        M0_mc = rbinom(NJ,1, callM(Mhat_z0_mc))
        M1_mc = rbinom(NJ,1, callM(Mhat_z1_mc))
      }

      # Normal model for continuous Y
      # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      matM_z0m0_mc = cbind(1, M0_mc, z0, X.testTSBB)
      matM_z1m0_mc = cbind(1, M0_mc, z1, X.testTSBB)
      matM_z1m1_mc = cbind(1, M1_mc, z1, X.testTSBB)
      Yhat_z0m0_mc = c(matM_z0m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m0_mc = c(matM_z1m0_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)
      Yhat_z1m1_mc = c(matM_z1m1_mc %*% YbetaPars + Ypsi_N) # rnorm(N) * Ypsi_N)

      if (object$typeY == "continuous") {
        YsigPars = sqrt(Ysig2Pars)
        Yhat_z0m0_mc = rnorm(NJ, Yhat_z0m0_mc, YsigPars)
        Yhat_z1m0_mc = rnorm(NJ, Yhat_z1m0_mc, YsigPars)
        Yhat_z1m1_mc = rnorm(NJ, Yhat_z1m1_mc, YsigPars)
      } else if (object$typeY == "binary") {
        Yhat_z0m0_mc = rbinom(NJ, 1, plogis(Yhat_z0m0_mc))
        Yhat_z1m0_mc = rbinom(NJ, 1, plogis(Yhat_z1m0_mc))
        Yhat_z1m1_mc = rbinom(NJ, 1, plogis(Yhat_z1m1_mc))
      }
      
      # # -----------------------------------------------------------------------------
      # YbetaPars0 = object0$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars0 = object0$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars0  = object0$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars0 = object0$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars0 = object0$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiParsC0  = object0$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars0 = object0$MCMCposteriors$YMpsiLists[[post_reps]]
      # meanYMpsi0 = object0$MCMCposteriors$meanYMpsiLists[[post_reps]]
      # covYMpsiPars0 = object0$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # YbetaPars1 = object1$MCMCposteriors$YbetaLists[[post_reps]]
      # Ysig2Pars1 = object1$MCMCposteriors$Ysig2Lists[[post_reps]]
      # # YpsiPars1  = object1$MCMCposteriors$YpsiLists[[post_reps]]
      # MbetaPars1 = object1$MCMCposteriors$MbetaLists[[post_reps]]
      # Msig2Pars1 = object1$MCMCposteriors$Msig2Lists[[post_reps]]
      # # MpsiPars1  = object1$MCMCposteriors$MpsiLists[[post_reps]]
      # YMpsiPars1 = object1$MCMCposteriors$YMpsiLists[[post_reps]]
      # meanYMpsi1 = object1$MCMCposteriors$meanYMpsiLists[[post_reps]]
      # covYMpsiPars1 = object1$MCMCposteriors$covYMpsiLists[[post_reps]]
      # 
      # # random-effects
      # YMpsi_J0 = rmvn_cpp(J, meanYMpsi0, covYMpsiPars0) # * sqrt(diag(YMpsiPars0)) # rbind(YpsiPars0,MpsiPars0)
      # Ypsi_N0 = YMpsi_J0[1,UindexTSBB]
      # Mpsi_N0 = YMpsi_J0[2,UindexTSBB]
      # 
      # YMpsi_J1 = rmvn_cpp(J, meanYMpsi1, covYMpsiPars1) # * sqrt(diag(YMpsiPars1)) # rbind(YpsiPars1,MpsiPars1)
      # Ypsi_N1 = YMpsi_J1[1,UindexTSBB]
      # Mpsi_N1 = YMpsi_J1[2,UindexTSBB]
      # 
      # # Normal model for continuous M
      # # draw M0_mc and M1_mc
      # Mhat_z0_mc = matX.testTSBB %*% MbetaPars0 + Mpsi_N0
      # Mhat_z1_mc = matX.testTSBB %*% MbetaPars1 + Mpsi_N1
      # 
      # if (object0$typeM == "continuous") {
      #   MsigPars0  = sqrt(Msig2Pars0)
      #   MsigPars1  = sqrt(Msig2Pars1)
      #   M0_mc = rnorm(NJ, Mhat_z0_mc, MsigPars0)
      #   M1_mc = rnorm(NJ, Mhat_z1_mc, MsigPars1)
      # } else if (object0$typeM == "binary") {
      #   if (is.null(object0$modelM)) {
      #     callM = plogis
      #   } else {
      #     callM = pnorm
      #   }
      #   M0_mc = rbinom(NJ,1, callM(Mhat_z0_mc))
      #   M1_mc = rbinom(NJ,1, callM(Mhat_z1_mc))
      # }
      # 
      # # Normal model for continuous Y
      # # draw Y_z0m0_mc, Y_z1m0_mc, Y_z1m1_mc
      # matM_m0_mc = cbind(1, M0_mc, X.testTSBB)
      # matM_m1_mc = cbind(1, M1_mc, X.testTSBB)
      # Yhat_z0m0_mc = c(matM_m0_mc %*% YbetaPars0)
      # Yhat_z1m0_mc = c(matM_m0_mc %*% YbetaPars1)
      # Yhat_z1m1_mc = c(matM_m1_mc %*% YbetaPars1)
      # 
      # if (object0$typeY == "continuous") {
      #   YsigPars0 = sqrt(Ysig2Pars0)
      #   YsigPars1 = sqrt(Ysig2Pars1)
      #   Yhat_z0m0_mc = rnorm(NJ, Yhat_z0m0_mc, YsigPars0)
      #   Yhat_z1m0_mc = rnorm(NJ, Yhat_z1m0_mc, YsigPars1)
      #   Yhat_z1m1_mc = rnorm(NJ, Yhat_z1m1_mc, YsigPars1)
      # } else if (object0$typeY == "binary") {
      #   Yhat_z0m0_mc = rbinom(NJ, 1, plogis(Yhat_z0m0_mc))
      #   Yhat_z1m0_mc = rbinom(NJ, 1, plogis(Yhat_z1m0_mc))
      #   Yhat_z1m1_mc = rbinom(NJ, 1, plogis(Yhat_z1m1_mc))
      # }
      
      # ------------------------------------------------------------------------------
      # Update parameters for cluster-level confounder (rho: J X 1 vector)
      # Update parameters for individual-level confounder (pi: n_j X 1 matrix, j=1,...J)
      # Update parameters for link cluster-individual (omega: J X J matrix)
      rhoPars = c(rdirichlet_cpp(1,rep(1,J)))
      piPars = sapply(1:J, function(l) rdirichlet_cpp(1, rep(1,n_j[l])))
      omePars = sapply(1:J, function(l) rdirichlet_cpp(1, xi_distPOST[,l]))
      
      RhoOmegaPi = numeric(NJ)
      for (l in 1:J) {
        RhoOmegaPitemp = numeric(N)
        for (j in 1:J) {
          ind_Uindex_j = which(Uindex == j)
          RhoOmegaPitemp[ind_Uindex_j] = omePars[j,l]*c(piPars[[j]])
        }
        RhoOmegaPi[((N*(l-1)+1):(N*l))] = rhoPars[l]*RhoOmegaPitemp
      }
      RhoOmegaPi = RhoOmegaPi * vSindexTSBB
      
      E_Y_z0m0_mc = sum(RhoOmegaPi *  Yhat_z0m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m0_mc = sum(RhoOmegaPi *  Yhat_z1m0_mc)/sum(RhoOmegaPi)
      E_Y_z1m1_mc = sum(RhoOmegaPi *  Yhat_z1m1_mc)/sum(RhoOmegaPi)
      
      # Causal Effects
      E_Y_z0m0[post_reps] = mean(E_Y_z0m0_mc)
      E_Y_z1m0[post_reps] = mean(E_Y_z1m0_mc)
      E_Y_z1m1[post_reps] = mean(E_Y_z1m1_mc)
      
      if (post_reps %% iter_check == 0){
        cat("Post-Processing",post_reps,"(",(post_reps/n_MCMC)*100,"%)","Time:",date(),"\n")
      }
    }
  }
  if (CE) {
    NIE = E_Y_z1m1 - E_Y_z1m0
    NDE = E_Y_z1m0 - E_Y_z0m0
    ATE = E_Y_z1m1 - E_Y_z0m0
  } else {
    NIE = E_Y_z1m1
    NDE = E_Y_z1m0
    ATE = E_Y_z0m0
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
  
  # Calculate median of posterior for NIE, NDE, and ATE
  NIE_result_mc = POSTsummary(NIE, "median", quantile_alpha)
  NDE_result_mc = POSTsummary(NDE, "median", quantile_alpha)
  ATE_result_mc = POSTsummary(ATE, "median", quantile_alpha)
  
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
# -----------------------------------------------------------------------------
# Define MCMC function for continuous Y and continuous M
GLMmediationMCMCconYconM = function(Y, M, Z, C, V, Uindex,
                                    gibbs_iter = gibbs_iter, 
                                    gibbs_burnin = gibbs_burnin, 
                                    gibbs_thin = gibbs_thin){
  
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
  
  # Set initial values ----------------------------------------------------------
  # with Z ----------------------------------------------------------------------
  matX = cbind(1, Z, C, V)
  matM = cbind(1, M, Z, C, V)
  # without Z -------------------------------------------------------------------
  # matX = cbind(1, C, V)
  # matM = cbind(1, M, C, V)
  
  pm = ncol(matX)
  py = ncol(matM)
  
  # -----------------------------------------------------------------------------
  # -----------------------------------------------------------------------------
  # --------------------------- Outcome & Mediator Fit --------------------------
  # -----------------------------------------------------------------------------
  # --------------------------------- Stan Code ---------------------------------
  stan_code = '
    data {
      // Data
      int<lower=1> J;  // Number of clusters
      int<lower=1> N;  // Number of observations
      int<lower=1> pm;  // Number of confounders for mediator glm
      int<lower=1> py;  // Number of confounders for outcome glm
      
      int<lower=1, upper=J> Uindex[N]; // Cluster assignment for each observation
      row_vector[pm] matX[N];  // design matrix for mediator glm
      row_vector[py] matM[N];  // design matrix for outcome glm
      real M[N];  // Mediator
      real Y[N];  // Outcome
      
      cov_matrix[2] B_covYMpsi;
    }
    parameters {
      real<lower=0> Msig2;
      vector[pm] Mbeta;
      real<lower=0> Ysig2;
      vector[py] Ybeta;
      
      array[J] vector[2] YMpsi;
      vector[2] meanYMpsi;
      // cov_matrix[2] covYMpsi;
      corr_matrix[2] corrYMpsi;
      vector<lower=0>[2] diagYMpsi;
    }
    transformed parameters {
      matrix[2, 2] covYMpsi = quad_form_diag(corrYMpsi, diagYMpsi);
    }
    model {
      // Priors for psi - cluster-level intercepts: mediator & outcome glm
      // hyper-priors
      meanYMpsi ~ multi_normal(rep_vector(0, 2), covYMpsi);
      // covYMpsi ~ inv_wishart(4, B_covYMpsi);
      // diagYMpsi ~ cauchy(0, 5);
      // corrYMpsi ~ lkj_corr(1);
      
      // priors
      YMpsi ~ multi_normal(meanYMpsi, covYMpsi);
      // Mbeta ~ normal(0, 5);
      // Msig2 ~ cauchy(0, 25); // , 2);
      // Ybeta ~ normal(0, 5);
      // Ysig2 ~ cauchy(0, 25); // , 2);
      
      // GLM likelihood
      vector[N] Mhat;
      vector[N] Yhat;
      for (i in 1:N) {
        Mhat[i] = matX[i] * Mbeta + YMpsi[Uindex[i],2];
        Yhat[i] = matM[i] * Ybeta + YMpsi[Uindex[i],1];
      }
      M ~ normal(Mhat, sqrt(Msig2));
      Y ~ normal(Yhat, sqrt(Ysig2));
    }
  '
  
  # -----------------------------------------------------------------------------
  lmeMtemp = lme(M~., random = ~ 1 |factor(Uindex), data.frame(M, C, V, Uindex))
  Msigest = summary(lmeMtemp)$sigma
  uM = c(lmeMtemp$coefficients$random[[1]])
  B_uM = sd(uM)
  # B_uM = sd(M)
  
  # -----------------------------------------------------------------------------
  lmeYtemp = lme(Y~., random = ~ 1 |factor(Uindex), data.frame(Y, M, C, V, Uindex))
  Ysigest = summary(lmeYtemp)$sigma
  uY = c(lmeYtemp$coefficients$random[[1]])
  B_uY = sd(uY)
  # B_uY = sd(Y)
  
  # -----------------------------------------------------------------------------
  # --------------------------------- Stan Data ---------------------------------
  stan_data = list(Y = Y, M = M, Uindex = Uindex, matX = matX, matM = matM)
  
  stan_data$J = length(unique(Uindex))
  stan_data$n_j = sapply(1:J, function(l) sum(Uindex==l))
  stan_data$N = sum(stan_data$n_j)
  stan_data$pm = ncol(stan_data$matX)
  stan_data$py = ncol(stan_data$matM)
  stan_data$B_covYMpsi = cbind(c(B_uM^{2},0),c(0,B_uY^{2}))
  
  # -----------------------------------------------------------------------------
  # ---------------------------------- Stan Fit ---------------------------------
  # Make vectors to store draws from Gibbs Sampler
  n_MCMC = floor(gibbs_iter/gibbs_thin)
  
  # Define the number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  gibbs_total = gibbs_iter + gibbs_burnin
  
  stanfit = stan(model_code = stan_code, data = stan_data,
                 iter = gibbs_total, warmup = gibbs_burnin, thin = gibbs_thin,
                 pars = c("Ybeta", "Ysig2", "Mbeta", "Msig2", "YMpsi", "meanYMpsi", "covYMpsi"),
                 chains = 1)
  # , control=list(adapt_delta=0.95)
  
  extractfit = rstan::extract(stanfit)
  
  # End initial values ----------------------------------------------------------
  
  # Make vectors to store draws from Gibbs Sampler
  YbetaLists = list(NA)
  Ysig2Lists = list(NA)
  YpsiLists  = list(NA)
  MbetaLists = list(NA)
  Msig2Lists = list(NA)
  MpsiLists  = list(NA)
  YMpsiLists = list(NA)
  meanYMpsiLists = list(NA)
  covYMpsiLists = list(NA)
  # End initial values ------------------------------------------------------------
  for (gibbs_reps in 1:n_MCMC) {
    # Save all parameters
    YbetaLists[[gibbs_reps]] = extractfit$Ybeta[gibbs_reps,]
    Ysig2Lists[[gibbs_reps]] = extractfit$Ysig2[gibbs_reps]
    YpsiLists[[gibbs_reps]]  = extractfit$Ypsi[gibbs_reps,]
    MbetaLists[[gibbs_reps]] = extractfit$Mbeta[gibbs_reps,]
    Msig2Lists[[gibbs_reps]] = extractfit$Msig2[gibbs_reps]
    MpsiLists[[gibbs_reps]]  = extractfit$Mpsi[gibbs_reps,]
    YMpsiLists[[gibbs_reps]] = extractfit$YMpsi[gibbs_reps,,]
    meanYMpsiLists[[gibbs_reps]] = extractfit$meanYMpsi[gibbs_reps,]
    covYMpsiLists[[gibbs_reps]] = extractfit$covYMpsi[gibbs_reps,,]
  }
  
  # constants
  constants = list(J = J, N = N, n_j = n_j,
                   py = py, pm = pm, 
                   n_MCMC = n_MCMC)
  
  # Observed Data
  obsdata = list(Y = Y, M = M, C = C, Z = Z, V = V, Uindex = Uindex)
  
  # MCMC Posteriors
  MCMCposteriors = list(Ysig2Lists = Ysig2Lists,
                        YbetaLists = YbetaLists,
                        YpsiLists = YpsiLists,
                        MbetaLists = MbetaLists,
                        Msig2Lists = Msig2Lists,
                        MpsiLists = MpsiLists,
                        YMpsiLists = YMpsiLists,
                        meanYMpsiLists = meanYMpsiLists,
                        covYMpsiLists = covYMpsiLists)
  
  # MCMC results
  object = list(obsdata = obsdata,
                constants = constants,
                MCMCposteriors = MCMCposteriors)
  
  object$typeM = "continuous"
  object$typeY = "continuous"
  class(object) = "rGLMmediation"
  
  return(object)
}

# Define MCMC function for continuous Y and binary M
GLMmediationMCMCconYbinM = function(Y, M, Z, C, V, Uindex,
                                    gibbs_iter = gibbs_iter, 
                                    gibbs_burnin = gibbs_burnin, 
                                    gibbs_thin = gibbs_thin,
                                    modelM = modelM){
  
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
  
  # Set initial values ----------------------------------------------------------
  # with Z ----------------------------------------------------------------------
  matX = cbind(1, Z, C, V)
  matM = cbind(1, M, Z, C, V)
  # without Z -------------------------------------------------------------------
  # matX = cbind(1, C, V)
  # matM = cbind(1, M, C, V)
  
  pm = ncol(matX)
  py = ncol(matM)
  
  # -----------------------------------------------------------------------------
  # -----------------------------------------------------------------------------
  # --------------------------- Outcome & Mediator Fit --------------------------
  # -----------------------------------------------------------------------------
  # --------------------------------- Stan Code ---------------------------------
  if (modelM == "probit") {
    stan_code = '
    data {
      // Data
      int<lower=1> J;  // Number of clusters
      int<lower=1> N;  // Number of observations
      int<lower=1> pm;  // Number of confounders for mediator glm
      int<lower=1> py;  // Number of confounders for outcome glm
      
      int<lower=1, upper=J> Uindex[N]; // Cluster assignment for each observation
      row_vector[pm] matX[N];  // design matrix for mediator glm
      row_vector[py] matM[N];  // design matrix for outcome glm
      int<lower=0, upper=1> M[N];  // Mediator
      real Y[N];  // Outcome
      
      cov_matrix[2] B_covYMpsi;
    }
    parameters {
      vector[pm] Mbeta;
      real<lower=0> Ysig2;
      vector[py] Ybeta;
      
      array[J] vector[2] YMpsi;
      vector[2] meanYMpsi;
      // cov_matrix[2] covYMpsi;
      corr_matrix[2] corrYMpsi;
      vector<lower=0>[2] diagYMpsi;
    }
    transformed parameters {
      matrix[2, 2] covYMpsi = quad_form_diag(corrYMpsi, diagYMpsi);
    }
    model {
      // Priors for psi - cluster-level intercepts: mediator & outcome glm
      // hyper-priors
      meanYMpsi ~ multi_normal(rep_vector(0, 2), covYMpsi);
      // covYMpsi ~ inv_wishart(4, B_covYMpsi);
      // diagYMpsi ~ cauchy(0, 5);
      // corrYMpsi ~ lkj_corr(1);
      
      // priors
      YMpsi ~ multi_normal(meanYMpsi, covYMpsi);
      // Mbeta ~ normal(0, 5);
      // Ybeta ~ normal(0, 5);
      // Ysig2 ~ cauchy(0, 25); // , 2);
      
      // GLM likelihood
      vector[N] Mhat;
      vector[N] Yhat;
      for (i in 1:N) {
        Mhat[i] = matX[i] * Mbeta + YMpsi[Uindex[i],2];
        Yhat[i] = matM[i] * Ybeta + YMpsi[Uindex[i],1];
      }
      M ~ bernoulli(Phi(Mhat));
      Y ~ normal(Yhat, sqrt(Ysig2));
    }
    '
  } else if (modelM == "logit") {
    stan_code = '
    data {
      // Data
      int<lower=1> J;  // Number of clusters
      int<lower=1> N;  // Number of observations
      int<lower=1> pm;  // Number of confounders for mediator glm
      int<lower=1> py;  // Number of confounders for outcome glm
      
      int<lower=1, upper=J> Uindex[N]; // Cluster assignment for each observation
      row_vector[pm] matX[N];  // design matrix for mediator glm
      row_vector[py] matM[N];  // design matrix for outcome glm
      int<lower=0, upper=1> M[N];  // Mediator
      real Y[N];  // Outcome
      
      cov_matrix[2] B_covYMpsi;
    }
    parameters {
      vector[pm] Mbeta;
      real<lower=0> Ysig2;
      vector[py] Ybeta;
      
      array[J] vector[2] YMpsi;
      vector[2] meanYMpsi;
      // cov_matrix[2] covYMpsi;
      corr_matrix[2] corrYMpsi;
      vector<lower=0>[2] diagYMpsi;
    }
    transformed parameters {
      matrix[2, 2] covYMpsi = quad_form_diag(corrYMpsi, diagYMpsi);
    }
    model {
      // Priors for psi - cluster-level intercepts: mediator & outcome glm
      // hyper-priors
      meanYMpsi ~ multi_normal(rep_vector(0, 2), covYMpsi);
      // covYMpsi ~ inv_wishart(4, B_covYMpsi);
      // diagYMpsi ~ cauchy(0, 5);
      // corrYMpsi ~ lkj_corr(1);
      
      // priors
      YMpsi ~ multi_normal(meanYMpsi, covYMpsi);
      // Mbeta ~ normal(0, 5);
      // Ybeta ~ normal(0, 5);
      // Ysig2 ~ cauchy(0, 25); // , 2);
      
      // GLM likelihood
      vector[N] Mhat;
      vector[N] Yhat;
      for (i in 1:N) {
        Mhat[i] = matX[i] * Mbeta + YMpsi[Uindex[i],2];
        Yhat[i] = matM[i] * Ybeta + YMpsi[Uindex[i],1];
      }
      M ~ bernoulli_logit(Mhat);
      Y ~ normal(Yhat, sqrt(Ysig2));
    }
    '
  }
  
  # -----------------------------------------------------------------------------
  # lmeMtemp = lme(M~., random = ~ 1 |factor(Uindex), data.frame(M, C, V, Uindex))
  # Msigest = summary(lmeMtemp)$sigma
  # uM = c(lmeMtemp$coefficients$random[[1]])
  # B_uM = sd(uM)
  B_uM = sd(M)
  
  # -----------------------------------------------------------------------------
  lmeYtemp = lme(Y~., random = ~ 1 |factor(Uindex), data.frame(Y, M, C, V, Uindex))
  Ysigest = summary(lmeYtemp)$sigma
  uY = c(lmeYtemp$coefficients$random[[1]])
  B_uY = sd(uY)
  # B_uY = sd(Y)
  
  # -----------------------------------------------------------------------------
  # --------------------------------- Stan Data ---------------------------------
  stan_data = list(Y = Y, M = M, Uindex = Uindex, matX = matX, matM = matM)
  
  stan_data$J = length(unique(Uindex))
  stan_data$n_j = sapply(1:J, function(l) sum(Uindex==l))
  stan_data$N = sum(stan_data$n_j)
  stan_data$pm = ncol(stan_data$matX)
  stan_data$py = ncol(stan_data$matM)
  stan_data$B_covYMpsi = cbind(c(B_uM^{2},0),c(0,B_uY^{2}))
  
  # -----------------------------------------------------------------------------
  # ---------------------------------- Stan Fit ---------------------------------
  # Make vectors to store draws from Gibbs Sampler
  n_MCMC = floor(gibbs_iter/gibbs_thin)
  
  # Define the number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  gibbs_total = gibbs_iter + gibbs_burnin
  
  stanfit = stan(model_code = stan_code, data = stan_data, 
                 iter = gibbs_total, warmup = gibbs_burnin, thin = gibbs_thin,
                 pars = c("Ybeta", "Ysig2", "Mbeta", "YMpsi", "meanYMpsi", "covYMpsi"),
                 chains = 1)
  # init = list(Ybeta = rep(0,py), Ysig2 = 1, Mbeta = rep(0,pm), 
  #             Ypsi = rep(0,J), Mpsi = rep(0,J), 
  #             diagYMpsi = rep(0,2), corrYMpsi = 0)
  # stanfit = stan(model_code = stan_code, data = stan_data, init = init, 
  #                iter = gibbs_total, warmup = gibbs_burnin, thin = gibbs_thin,
  #                pars = c("Ybeta", "Ysig2", "Ypsi", "Mbeta", "Mpsi", "YMpsi", "meanYMpsi", "covYMpsi"),
  #                chains = 1)
  # , control=list(adapt_delta=0.95)
  
  extractfit = rstan::extract(stanfit)
  
  # End initial values ----------------------------------------------------------
  
  # Make vectors to store draws from Gibbs Sampler
  YbetaLists = list(NA)
  Ysig2Lists = list(NA)
  YpsiLists  = list(NA)
  MbetaLists = list(NA)
  MpsiLists  = list(NA)
  YMpsiLists = list(NA)
  meanYMpsiLists = list(NA)
  covYMpsiLists = list(NA)
  # End initial values ------------------------------------------------------------
  for (gibbs_reps in 1:n_MCMC) {
    # Save all parameters
    YbetaLists[[gibbs_reps]] = extractfit$Ybeta[gibbs_reps,]
    Ysig2Lists[[gibbs_reps]] = extractfit$Ysig2[gibbs_reps]
    YpsiLists[[gibbs_reps]]  = extractfit$Ypsi[gibbs_reps,]
    MbetaLists[[gibbs_reps]] = extractfit$Mbeta[gibbs_reps,]
    MpsiLists[[gibbs_reps]]  = extractfit$Mpsi[gibbs_reps,]
    YMpsiLists[[gibbs_reps]] = extractfit$YMpsi[gibbs_reps,,]
    meanYMpsiLists[[gibbs_reps]] = extractfit$meanYMpsi[gibbs_reps,]
    covYMpsiLists[[gibbs_reps]] = extractfit$covYMpsi[gibbs_reps,,]
  }
  
  # constants
  constants = list(J = J, N = N, n_j = n_j,
                   py = py, pm = pm, 
                   n_MCMC = n_MCMC)
  
  # Observed Data
  obsdata = list(Y = Y, M = M, C = C, Z = Z, V = V, Uindex = Uindex)
  
  # MCMC Posteriors
  MCMCposteriors = list(Ysig2Lists = Ysig2Lists,
                        YbetaLists = YbetaLists,
                        YpsiLists = YpsiLists,
                        MbetaLists = MbetaLists,
                        MpsiLists = MpsiLists,
                        YMpsiLists = YMpsiLists,
                        meanYMpsiLists = meanYMpsiLists,
                        covYMpsiLists = covYMpsiLists)
  
  # MCMC results
  object = list(obsdata = obsdata,
                constants = constants,
                MCMCposteriors = MCMCposteriors)
  
  object$typeM = "binary"
  object$typeY = "continuous"
  object$modelM = modelM
  class(object) = "rGLMmediation"
  
  return(object)
}

# Define MCMC function for binary Y and binary M
GLMmediationMCMCbinYbinM = function(Y, M, Z, C, V, Uindex,
                                    gibbs_iter = gibbs_iter, 
                                    gibbs_burnin = gibbs_burnin, 
                                    gibbs_thin = gibbs_thin){
  
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
  
  # Set initial values ----------------------------------------------------------
  # with Z ----------------------------------------------------------------------
  matX = cbind(1, Z, C, V)
  matM = cbind(1, M, Z, C, V)
  # without Z -------------------------------------------------------------------
  # matX = cbind(1, C, V)
  # matM = cbind(1, M, C, V)
  
  pm = ncol(matX)
  py = ncol(matM)
  
  # -----------------------------------------------------------------------------
  # -----------------------------------------------------------------------------
  # --------------------------- Outcome & Mediator Fit --------------------------
  # -----------------------------------------------------------------------------
  # --------------------------------- Stan Code ---------------------------------
  stan_code = '
    data {
      // Data
      int<lower=1> J;  // Number of clusters
      int<lower=1> N;  // Number of observations
      int<lower=1> pm;  // Number of confounders for mediator glm
      int<lower=1> py;  // Number of confounders for outcome glm
      
      int<lower=1, upper=J> Uindex[N]; // Cluster assignment for each observation
      row_vector[pm] matX[N];  // design matrix for mediator glm
      row_vector[py] matM[N];  // design matrix for outcome glm
      int<lower=0, upper=1> M[N];  // Mediator
      int<lower=0, upper=1> Y[N];  // Outcome
      
      cov_matrix[2] B_covYMpsi;
    }
    parameters {
      vector[pm] Mbeta;
      vector[py] Ybeta;
      
      array[J] vector[2] YMpsi;
      vector[2] meanYMpsi;
      // cov_matrix[2] covYMpsi;
      corr_matrix[2] corrYMpsi;
      vector<lower=0>[2] diagYMpsi;
    }
    transformed parameters {
      matrix[2, 2] covYMpsi = quad_form_diag(corrYMpsi, diagYMpsi);
    }
    model {
      // Priors for psi - cluster-level intercepts: mediator & outcome glm
      // hyper-priors
      meanYMpsi ~ multi_normal(rep_vector(0, 2), covYMpsi);
      // covYMpsi ~ inv_wishart(4, B_covYMpsi);
      // diagYMpsi ~ cauchy(0, 5);
      // corrYMpsi ~ lkj_corr(1);
      
      // priors
      YMpsi ~ multi_normal(meanYMpsi, covYMpsi);
      // Mbeta ~ normal(0, 5);
      // Ybeta ~ normal(0, 5);
      
      // GLM likelihood
      vector[N] Mhat;
      vector[N] Yhat;
      for (i in 1:N) {
        Mhat[i] = matX[i] * Mbeta + YMpsi[Uindex[i],2];
        Yhat[i] = matM[i] * Ybeta + YMpsi[Uindex[i],1];
      }
      M ~ bernoulli_logit(Mhat);
      Y ~ bernoulli_logit(Yhat);
    }
  '
  
  # -----------------------------------------------------------------------------
  # lmeMtemp = lme(M~., random = ~ 1 |factor(Uindex), data.frame(M, C, V, Uindex))
  # Msigest = summary(lmeMtemp)$sigma
  # uM = c(lmeMtemp$coefficients$random[[1]])
  # B_uM = sd(uM)
  B_uM = sd(M)
  
  # -----------------------------------------------------------------------------
  # lmeYtemp = lme(Y~., random = ~ 1 |factor(Uindex), data.frame(Y, M, C, V, Uindex))
  # Ysigest = summary(lmeYtemp)$sigma
  # uY = c(lmeYtemp$coefficients$random[[1]])
  # B_uY = sd(uY)
  B_uY = sd(Y)
  
  # -----------------------------------------------------------------------------
  # --------------------------------- Stan Data ---------------------------------
  stan_data = list(Y = Y, M = M, Uindex = Uindex, matX = matX, matM = matM)
  
  stan_data$J = length(unique(Uindex))
  stan_data$n_j = sapply(1:J, function(l) sum(Uindex==l))
  stan_data$N = sum(stan_data$n_j)
  stan_data$pm = ncol(stan_data$matX)
  stan_data$py = ncol(stan_data$matM)
  stan_data$B_covYMpsi = cbind(c(B_uM^{2},0),c(0,B_uY^{2}))
  
  # -----------------------------------------------------------------------------
  # ---------------------------------- Stan Fit ---------------------------------
  # Make vectors to store draws from Gibbs Sampler
  n_MCMC = floor(gibbs_iter/gibbs_thin)
  
  # Define the number of Markov Chain Monte Carlo (MCMC) draws in Gibbs Sampler
  gibbs_total = gibbs_iter + gibbs_burnin
  
  stanfit = stan(model_code = stan_code, data = stan_data, 
                 iter = gibbs_total, warmup = gibbs_burnin, thin = gibbs_thin,
                 pars = c("Ybeta", "Mbeta", "YMpsi", "meanYMpsi", "covYMpsi"),
                 chains = 1)
  
  extractfit = rstan::extract(stanfit)
  
  # End initial values ----------------------------------------------------------
  
  # Make vectors to store draws from Gibbs Sampler
  YbetaLists = list(NA)
  YpsiLists  = list(NA)
  MbetaLists = list(NA)
  MpsiLists  = list(NA)
  YMpsiLists = list(NA)
  meanYMpsiLists = list(NA)
  covYMpsiLists = list(NA)
  # End initial values ------------------------------------------------------------
  for (gibbs_reps in 1:n_MCMC) {
    # Save all parameters
    YbetaLists[[gibbs_reps]] = extractfit$Ybeta[gibbs_reps,]
    YpsiLists[[gibbs_reps]]  = extractfit$Ypsi[gibbs_reps,]
    MbetaLists[[gibbs_reps]] = extractfit$Mbeta[gibbs_reps,]
    MpsiLists[[gibbs_reps]]  = extractfit$Mpsi[gibbs_reps,]
    YMpsiLists[[gibbs_reps]] = extractfit$YMpsi[gibbs_reps,,]
    meanYMpsiLists[[gibbs_reps]] = extractfit$meanYMpsi[gibbs_reps,]
    covYMpsiLists[[gibbs_reps]] = extractfit$covYMpsi[gibbs_reps,,]
  }
  
  # constants
  constants = list(J = J, N = N, n_j = n_j,
                   py = py, pm = pm, 
                   n_MCMC = n_MCMC)
  
  # Observed Data
  obsdata = list(Y = Y, M = M, C = C, Z = Z, V = V, Uindex = Uindex)
  
  # MCMC Posteriors
  MCMCposteriors = list(YbetaLists = YbetaLists,
                        YpsiLists = YpsiLists,
                        MbetaLists = MbetaLists,
                        MpsiLists = MpsiLists,
                        YMpsiLists = YMpsiLists,
                        meanYMpsiLists = meanYMpsiLists,
                        covYMpsiLists = covYMpsiLists)
  
  # MCMC results
  object = list(obsdata = obsdata,
                constants = constants,
                MCMCposteriors = MCMCposteriors)
  
  object$typeM = "binary"
  object$typeY = "binary"
  class(object) = "rGLMmediation"
  
  return(object)
}

# Define MCMC function for any Y and any M
GLMmediation = function(Y, M, Z, C, V, Uindex,
                        typeY = "continuous", typeM = "continuous", 
                        gibbs_iter = 2e3, gibbs_burnin = 2e3, gibbs_thin = 10,
                        modelM = "probit"){
  if (typeY == "continuous") {
    if (typeM == "continuous") {
      object = GLMmediationMCMCconYconM(Y, M, Z, C, V, Uindex,
                                        gibbs_iter,
                                        gibbs_burnin,
                                        gibbs_thin)
    } else if (typeM == "binary") {
      object = GLMmediationMCMCconYbinM(Y, M, Z, C, V, Uindex,
                                        gibbs_iter = gibbs_iter,
                                        gibbs_burnin = gibbs_burnin,
                                        gibbs_thin = gibbs_thin,
                                        modelM = modelM)
    } else if (typeM == "ordinal") {
      object = GLMmediationMCMCconYordM(Y, M, Z, C, V, Uindex,
                                        gibbs_iter = gibbs_iter,
                                        gibbs_burnin = gibbs_burnin,
                                        gibbs_thin = gibbs_thin)
    } else if (typeM == "count") {
      object = GLMmediationMCMCconYcouM(Y, M, Z, C, V, Uindex,
                                        gibbs_iter = gibbs_iter,
                                        gibbs_burnin = gibbs_burnin,
                                        gibbs_thin = gibbs_thin)
    }
  } else if (typeY == "binary") {
    if (typeM == "continuous") {
      #
    } else if (typeM == "binary") {
      object = GLMmediationMCMCbinYbinM(Y, M, Z, C, V, Uindex,
                                        gibbs_iter = gibbs_iter,
                                        gibbs_burnin = gibbs_burnin,
                                        gibbs_thin = gibbs_thin)
    } else if (typeM == "ordinal") {
      #
    } else if (typeM == "count") {
      #
    }
  }
  # else if (typeY == "ordinal") ㄴ{
  #
  # } else if (typeY == "count") {
  #
  # }
  
  # Z0 = 0
  # which_Z0 = which(Z == Z0)
  # Y0 = Y[which_Z0]
  # M0 = M[which_Z0]
  # C0 = C[which_Z0,]
  # V0 = V[which_Z0,]
  # Uindex0 = Uindex[which_Z0]
  # matX0 = cbind(C0, V0)
  # matM0 = cbind(M0, C0, V0)
  # N0 = length(Y0)
  # uniqueUindex0 = sort(unique(Uindex0))
  # J0 = length(uniqueUindex0)
  # Uindex0 = apply(sapply(1:J0, function(l) ifelse(Uindex0==uniqueUindex0[l],l,0)),1,sum)
  # n_j0 = as.numeric(table(Uindex0))
  # order_Uindex0 = order(Uindex0)
  # Y0 = Y0[order_Uindex0]
  # M0 = M0[order_Uindex0]
  # C0 = C0[order_Uindex0,]
  # V0 = V0[order_Uindex0,]
  # Uindex0 = Uindex0[order_Uindex0]
  # 
  # Z1 = 1
  # which_Z1 = which(Z == Z1)
  # Y1 = Y[which_Z1]
  # M1 = M[which_Z1]
  # C1 = C[which_Z1,]
  # V1 = V[which_Z1,]
  # Uindex1 = Uindex[which_Z1]
  # matX1 = cbind(C1, V1)
  # matM1 = cbind(M1, C1, V1)
  # N1 = length(Y1)
  # uniqueUindex1 = sort(unique(Uindex1))
  # J1 = length(uniqueUindex1)
  # Uindex1 = apply(sapply(1:J1, function(l) ifelse(Uindex1==uniqueUindex1[l],l,0)),1,sum)
  # n_j1 = as.numeric(table(Uindex1))
  # order_Uindex1 = order(Uindex1)
  # Y1 = Y1[order_Uindex1]
  # M1 = M1[order_Uindex1]
  # C1 = C1[order_Uindex1,]
  # V1 = V1[order_Uindex1,]
  # Uindex1 = Uindex1[order_Uindex1]
  # 
  # if (typeY == "continuous") {
  #   if (typeM == "continuous") {
  #     object0 = GLMmediationMCMCconYconM(Y0, M0, rep(Z0, N0), C0, V0, Uindex0,
  #                                        gibbs_iter,
  #                                        gibbs_burnin,
  #                                        gibbs_thin)
  #     object1 = GLMmediationMCMCconYconM(Y1, M1, rep(Z1, N1), C1, V1, Uindex1,
  #                                        gibbs_iter,
  #                                        gibbs_burnin,
  #                                        gibbs_thin)
  #   } else if (typeM == "binary") {
  #     object0 = GLMmediationMCMCconYbinM(Y0, M0, rep(Z0, N0), C0, V0, Uindex0,
  #                                        gibbs_iter = gibbs_iter,
  #                                        gibbs_burnin = gibbs_burnin,
  #                                        gibbs_thin = gibbs_thin,
  #                                        modelM = modelM)
  #     object1 = GLMmediationMCMCconYbinM(Y1, M1, rep(Z1, N1), C1, V1, Uindex1,
  #                                        gibbs_iter = gibbs_iter,
  #                                        gibbs_burnin = gibbs_burnin,
  #                                        gibbs_thin = gibbs_thin,
  #                                        modelM = modelM)
  #   } else if (typeM == "ordinal") {
  #     object0 = GLMmediationMCMCconYordM(Y0, M0, rep(Z0, N0), C0, V0, Uindex0,
  #                                        gibbs_iter = gibbs_iter,
  #                                        gibbs_burnin = gibbs_burnin,
  #                                        gibbs_thin = gibbs_thin)
  #     object1 = GLMmediationMCMCconYordM(Y1, M1, rep(Z1, N1), C1, V1, Uindex1,
  #                                        gibbs_iter = gibbs_iter,
  #                                        gibbs_burnin = gibbs_burnin,
  #                                        gibbs_thin = gibbs_thin)
  #   } else if (typeM == "count") {
  #     object0 = GLMmediationMCMCconYcouM(Y0, M0, rep(Z0, N0), C0, V0, Uindex0,
  #                                        gibbs_iter = gibbs_iter,
  #                                        gibbs_burnin = gibbs_burnin,
  #                                        gibbs_thin = gibbs_thin)
  #     object1 = GLMmediationMCMCconYcouM(Y1, M1, rep(Z1, N1), C1, V1, Uindex1,
  #                                        gibbs_iter = gibbs_iter,
  #                                        gibbs_burnin = gibbs_burnin,
  #                                        gibbs_thin = gibbs_thin)
  #   }
  # } else if (typeY == "binary") {
  #   if (typeM == "continuous") {
  #     #
  #   } else if (typeM == "binary") {
  #     object0 = GLMmediationMCMCbinYbinM(Y0, M0, rep(Z0, N0), C0, V0, Uindex0,
  #                                        gibbs_iter = gibbs_iter,
  #                                        gibbs_burnin = gibbs_burnin,
  #                                        gibbs_thin = gibbs_thin)
  #     object1 = GLMmediationMCMCbinYbinM(Y1, M1, rep(Z1, N1), C1, V1, Uindex1,
  #                                        gibbs_iter = gibbs_iter,
  #                                        gibbs_burnin = gibbs_burnin,
  #                                        gibbs_thin = gibbs_thin)
  #   } else if (typeM == "ordinal") {
  #     #
  #   } else if (typeM == "count") {
  #     #
  #   }
  # }
  # # else if (typeY == "ordinal") ㄴ{
  # #
  # # } else if (typeY == "count") {
  # #
  # # }
  # 
  # object = list(object0 = object0, object1 = object1)
  
  object$modelM = modelM
  class(object) = "rGLMmediation"
  
  return(object)
}

# -----------------------------------------------------------------------------
# End function definitions ----------------------------------------------------
# -----------------------------------------------------------------------------
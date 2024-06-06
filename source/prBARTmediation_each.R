prBARTmediation_each = function(object0, # object for Z=0 from rBARTmediation
                                object1, # object for Z=1 from rBARTmediation
                                X.test,  # matrix X to predict at
                                Uindex0,
                                Uindex1){
  
  # --------------------------------------------------
  mc.cores = 1
  
  N = nrow(X.test)
  J0 = ncol(object0$uMdraw)
  J1 = ncol(object1$uMdraw)
  maxJ = max(J0, J1)
  n_MCMC = nrow(object0$uMdraw)
  
  X.test <- t(bartModelMatrix(X.test))
  
  pm <- length(object0$matXtreedraws$cutpoints)
  if(pm!=nrow(X.test)) {
    stop(paste0('The number of columns in X.test must be equal to ', pm))
  }
  
  M0res = .Call("cprBART", object0$matXtreedraws, X.test, mc.cores)$yhat.test + object0$Moffset
  M1res = .Call("cprBART", object1$matXtreedraws, X.test, mc.cores)$yhat.test + object1$Moffset
  for (j in 1:maxJ) {
    whichUindex0 = which(Uindex0==j)
    whichUindex1 = which(Uindex1==j)
    if (length(whichUindex0)){
      M0res[,whichUindex0] = M0res[,whichUindex0] + object0$uMdraw[,j]
    } else if (length(whichUindex1)){
      M1res[,whichUindex1] = M1res[,whichUindex1] + object1$uMdraw[,j]
    }
  }
  M0.test = t(sapply(1:n_MCMC, function(l) rnorm(N, M0res[l,], object0$iMsigest[l])))
  M1.test = t(sapply(1:n_MCMC, function(l) rnorm(N, M1res[l,], object0$iMsigest[l])))
  M0.test = apply(M0.test, 2, mean)
  M1.test = apply(M1.test, 2, mean)
  
  # --------------------------------------------------
  matM0.test = rbind(M0.test,X.test)
  matM1.test = rbind(M1.test,X.test)
  
  Yz0m0res = .Call("cprBART", object0$matMtreedraws, matM0.test, mc.cores)$yhat.test + object0$Yoffset
  Yz1m0res = .Call("cprBART", object1$matMtreedraws, matM0.test, mc.cores)$yhat.test + object1$Yoffset
  Yz1m1res = .Call("cprBART", object1$matMtreedraws, matM1.test, mc.cores)$yhat.test + object1$Yoffset
  for (j in 1:maxJ) {
    whichUindex0 = which(Uindex0==j)
    whichUindex1 = which(Uindex1==j)
    if (length(whichUindex0)){
      M0res[,whichUindex0] = M0res[,whichUindex0] + object0$uMdraw[,j]
    } else if (length(whichUindex1)){
      M1res[,whichUindex1] = M1res[,whichUindex1] + object1$uMdraw[,j]
    }
  }
  Yz0m0.test = Yz0m0res
  Yz1m0.test = Yz1m0res
  Yz1m1.test = Yz1m1res
  # Yz0m0.test = t(sapply(1:n_MCMC, function(l) rnorm(N, Yz0m0.test[l,], object$iYsigest[l])))
  # Yz1m0.test = t(sapply(1:n_MCMC, function(l) rnorm(N, Yz1m0.test[l,], object$iYsigest[l])))
  # Yz1m1.test = t(sapply(1:n_MCMC, function(l) rnorm(N, Yz1m1.test[l,], object$iYsigest[l])))
  
  # quantile(Yz1m1.test-Yz1m0.test, c(0.025, 0.5, 0.975));mean(Yz1m1.test-Yz1m0.test)
  # quantile(Yz1m0.test-Yz0m0.test, c(0.025, 0.5, 0.975));mean(Yz1m0.test-Yz0m0.test)
  # quantile(Yz1m1.test-Yz0m0.test, c(0.025, 0.5, 0.975));mean(Yz1m1.test-Yz0m0.test)
  # c(E_true[1]-E_true[2],E_true[2]-E_true[3],E_true[1]-E_true[3])
  
  return(list(Yz0m0.test=Yz0m0.test,
              Yz1m0.test=Yz1m0.test,
              Yz1m1.test=Yz1m1.test))
}
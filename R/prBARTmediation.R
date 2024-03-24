
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2018 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

prBARTmediation = function(object,  # object from rBARTmediation
                           X.test,  # matrix X to predict at
                           Uindex){
  # --------------------------------------------------
  mc.cores = 1
  
  z0 = 0
  z1 = 1
  
  N = nrow(X.test)
  J = ncol(object$uMdraw)
  n_MCMC = nrow(object$uMdraw)
  
  matXz0.test <- t(bartModelMatrix(cbind(z0,X.test)))
  matXz1.test <- t(bartModelMatrix(cbind(z1,X.test)))
  
  pm <- length(object$matXtreedraws$cutpoints)
  if(pm!=nrow(matXz0.test)) {
    stop(paste0('The number of columns in matX.test must be equal to ', pm))
  }
  
  M0res = .Call("cprBART", object$matXtreedraws, matXz0.test, mc.cores)$yhat.test + object$Moffset
  M1res = .Call("cprBART", object$matXtreedraws, matXz1.test, mc.cores)$yhat.test + object$Moffset
  for (j in 1:J) {
    whichUindex = which(Uindex==j)
    if(length(whichUindex)>0){
      # Mreff_tmp = rnorm(1)*(object$sd.uM)
      Mreff_tmp = object$uMdraw[,j]
      M0res[,whichUindex] = M0res[,whichUindex] + Mreff_tmp
      M1res[,whichUindex] = M1res[,whichUindex] + Mreff_tmp
    }
  }
  if(object$typeM == "continuous"){
    M0.test = apply(M0res, 2, mean)
    M1.test = apply(M1res, 2, mean)
    M0.test = rnorm(n_MCMC, M0.test, object$iMsigest)
    M1.test = rnorm(n_MCMC, M1.test, object$iMsigest)
  } else if(object$typeM == "binary"){
    M0.test = apply(pnorm(M0res), 2, mean)
    M1.test = apply(pnorm(M1res), 2, mean)
    M0.test = rbinom(n_MCMC, 1, M0.test)
    M1.test = rbinom(n_MCMC, 1, M1.test)
  } else if(object$typeM == "multinomial"){
    #
  }
  
  # --------------------------------------------------
  matM0z0.test = rbind(M0.test, matXz0.test)
  matM0z1.test = rbind(M0.test, matXz1.test)
  matM1z1.test = rbind(M1.test, matXz1.test)
  
  py <- length(object$matMtreedraws$cutpoints)
  if(py!=nrow(matM0z0.test)) {
    stop(paste0('The number of columns in matM.test must be equal to ', py))
  }
  
  Yz0m0res = .Call("cprBART", object$matMtreedraws, matM0z0.test, mc.cores)$yhat.test + object$Yoffset
  Yz1m0res = .Call("cprBART", object$matMtreedraws, matM0z1.test, mc.cores)$yhat.test + object$Yoffset
  Yz1m1res = .Call("cprBART", object$matMtreedraws, matM1z1.test, mc.cores)$yhat.test + object$Yoffset
  for (j in 1:J) {
    whichUindex = which(Uindex==j)
    if(length(whichUindex)>0){
      # Yreff_tmp = rnorm(1)*(object$sd.uY)
      Yreff_tmp = object$uYdraw[,j]
      Yz0m0res[,whichUindex] = Yz0m0res[,whichUindex] + Yreff_tmp
      Yz1m0res[,whichUindex] = Yz1m0res[,whichUindex] + Yreff_tmp
      Yz1m1res[,whichUindex] = Yz1m1res[,whichUindex] + Yreff_tmp
    }
  }
  if(object$typeY == "continuous"){
    Yz0m0.test = sapply(1:N, function(i) rnorm(n_MCMC, Yz0m0res[,i], object$iYsigest))
    Yz1m0.test = sapply(1:N, function(i) rnorm(n_MCMC, Yz1m0res[,i], object$iYsigest))
    Yz1m1.test = sapply(1:N, function(i) rnorm(n_MCMC, Yz1m1res[,i], object$iYsigest))
  } else if(object$typeY == "binary"){
    Yz0m0.test = pnorm(Yz0m0res)
    Yz1m0.test = pnorm(Yz1m0res)
    Yz1m1.test = pnorm(Yz1m1res)
    Yz0m0.test = sapply(1:N, function(i) rbinom(n_MCMC, 1, Yz0m0.test[,i]))
    Yz1m0.test = sapply(1:N, function(i) rbinom(n_MCMC, 1, Yz1m0.test[,i]))
    Yz1m1.test = sapply(1:N, function(i) rbinom(n_MCMC, 1, Yz1m1.test[,i]))
  } else if(object$typeY == "multinomial"){
    #
  }
  
  return(list(Yz0m0.test=Yz0m0.test,
              Yz1m0.test=Yz1m0.test,
              Yz1m1.test=Yz1m1.test))
}

## BART: Bayesian Additive Regression Trees
## Copyright (C) 2019 Rodney Sparapani

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

rBART = function(matX, M, 
                 Uindex=NULL, B_u=NULL,
                 sparse=FALSE, theta=0, omega=1,
                 a=0.5, b=1, augment=FALSE, rho=NULL,
                 xinfo=matrix(0,0,0), usequants=FALSE,
                 rm.const=TRUE,
                 sigest=NA, sigdf=3, sigquant=0.90,
                 k=2, power=2, base=0.95,
                 lambda=NA, tau.num=NA,
                 offset = mean(M),
                 ntree=200L, numcut=100L,
                 ndpost=1000L, nskip=100L, keepevery=1L,
                 printevery=100L, transposed=FALSE,
                 hostname=FALSE,
                 mc.cores = 1L, nice = 19L, seed = 99L){
  
  #--------------------------------------------------
  # data
  n = length(M)
  
  if(length(Uindex)==0)
    stop("the random effects indices must be provided")
  
  c.index=integer(n) ## changing from R/arbitrary indexing to C/0
  u.index=unique(Uindex)
  for(i in 1:n) {
    c.index[i]=which(Uindex[i]==u.index)-1
  }
  u.index=unique(c.index)
  J=length(u.index)
  n.j.vec=integer(J) ## n_j for each j
  for(j in 1:J) {
    n.j.vec[j]=length(which(u.index[j]==c.index))
  }
  
  if(!transposed) {
    temp = bartModelMatrix(matX, numcut, usequants=usequants,
                           xinfo=xinfo, rm.const=rm.const)
    matX = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    
    rm.const <- temp$rm.const
    rm(temp)
  } else {
    rm.const <- NULL
  }
  
  if(n!=ncol(matX)){
    stop('The length of M and the number of rows in matX must be identical')
  }
  
  p = nrow(matX)
  if(length(rho)==0) rho=p
  if(length(rm.const)==0) rm.const <- 1:p
  
  u=double(J)
  u[1]=NaN
  #--------------------------------------------------
  #prior
  nu = sigdf
  if(is.na(lambda)) {
    if(is.na(sigest)) {
      if(p < n) {
        temp = lme(M~., random=~1|factor(Uindex),
                   data.frame(t(matX),Uindex,M))
        sigest = summary(temp)$sigma
        u = c(temp$coefficients$random[[1]])
        if(length(B_u)==0) {
          B_u=2*sd(u)
        }
      } else {
        sigest = sd(M)
      }
    }
    qchi = qchisq(1.0-sigquant,nu)
    lambda = (sigest*sigest*qchi)/nu # lambda parameter for sigma prior
  } else {
    sigest=sqrt(lambda)
  }
  if(is.na(tau.num)) {
    tau=(max(M)-min(M))/(2*k*sqrt(ntree))
  } else {
    tau = tau.num/sqrt(ntree)
  }
  #--------------------------------------------------
  if(length(B_u)==0) {
    B_u <- sigest
  }
  
  if(.Platform$OS.type!='unix') {
    hostname <- FALSE
  } else if(hostname) {
    hostname <- system('hostname', intern=TRUE)
  }
  #--------------------------------------------------
  ptm <- proc.time()
  res = .Call("crBART",
              n,         # number of observations in training data
              p,         # dimension of x
              matX,   # p x n training data matX
              M,      # 1 x n training data M
              c.index,
              n.j.vec,
              u,         # random effects, if estimated
              J,
              B_u,
              ntree,
              numcut,
              ndpost*keepevery,
              nskip,
              keepevery,
              power,
              base,
              offset,
              tau,
              nu,
              lambda,
              sigest,
              sparse,
              theta,
              omega,
              a,
              b,
              rho,
              augment,
              printevery,
              xinfo)
  res$proc.time <- proc.time()-ptm
  #--------------------------------------------------
  
  res$mhat.mean <- apply(res$m.draw, 2, mean)
  
  names(res$treedraws$cutpoints) = dimnames(matX)[[1]]
  dimnames(res$varcount)[[2]] = as.list(dimnames(matX)[[1]])
  dimnames(res$varprob)[[2]] = as.list(dimnames(matX)[[1]])
  res$varcount.mean <- apply(res$varcount, 2, mean)
  res$varprob.mean <- apply(res$varprob, 2, mean)
  res$hostname <- hostname
  
  res$offset = offset
  res$rm.const <- rm.const
  res$sigest <- sigest
  res$B_u <- B_u
  res$u <- u
  attr(res, 'class') <- "rBART"
  
  return(res)
}

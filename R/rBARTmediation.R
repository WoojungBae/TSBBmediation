
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

rBARTmediation = function(Y, M, Z, C, V, Uindex=NULL, 
                          typeY = "continuous",
                          typeM = "continuous",
                          B_uM=NULL, B_uY=NULL,
                          sparse=FALSE, theta=0, omega=1,
                          a=0.5, b=1, augment=FALSE, 
                          matXrho=NULL, matMrho=NULL,
                          xinfo=matrix(0,0,0), usequants=FALSE,
                          matXrm.const=TRUE, matMrm.const=TRUE,
                          Msigest=NA, Ysigest=NA, 
                          sigdf=3, sigquant=0.90,
                          k=2, power=2, base=0.95,
                          Mlambda=NA, Ylambda=NA,  
                          Mtau.num=NA, Ytau.num=NA,
                          Moffset = NULL, Yoffset = NULL,
                          ntree=200L, 
                          matXnumcut=100L, matMnumcut=100L,
                          ndpost=1e3, nskip=1e4, keepevery=1e1,
                          printevery = (ndpost*keepevery)/10,
                          transposed=FALSE, hostname=FALSE,
                          mc.cores = 1L, nice = 19L, seed = 99L){
  # --------------------------------------------------
  # data
  ntypeY = as.integer(factor(typeY, levels = c("continuous", "binary", "multinomial")))
  ntypeM = as.integer(factor(typeM, levels = c("continuous", "binary", "multinomial")))
  if(is.na(ntypeY) || is.na(ntypeM)){
    stop("type argument must be set to either 'continuous', 'binary' or 'multinomial'")
  }
  
  if(is.null(Yoffset)){
    if(typeY == "continuous"){
      Yoffset <- mean(Y)
    } else if(typeY == "binary"){
      Yoffset <- qnorm(mean(Y))
      # Yoffset <- qlogis(Yoffset)
    } else if(typeY == "multinomial"){
      Yoffset <- 0
      stop("in development: 'multinomial'")
    }
  }
  if(is.null(Moffset)){
    if(typeM == "continuous"){
      Moffset <- mean(M)
    } else if(typeM == "binary"){
      Moffset <- qnorm(mean(M))
      # Moffset <- qlogis(Moffset)
    } else if(typeM == "multinomial"){
      Moffset <- 0
    }
  }
  
  n = length(Y)
  X = cbind(C, V)
  matX = cbind(Z, X)
  matM = cbind(M, matX)
  
  if(length(Uindex)==0){
    stop("the random effects indices must be provided")
  }
  
  u.index=unique(Uindex)
  u0.index=integer(n) ## changing from R/arbitrary indexing to C/0
  for(i in 1:n) {
    u0.index[i]=which(Uindex[i]==u.index)-1
  }
  u.index=unique(u0.index)
  J=length(u.index)
  n.j.vec=integer(J) ## n_j for each j
  for(j in 1:J) {
    n.j.vec[j]=length(which(u.index[j]==u0.index))
  }
  
  # --------------------------------------------------
  tmp.order = order(u0.index)
  matX = matX[tmp.order,]
  matM = matM[tmp.order,]
  M = M[tmp.order]
  Y = Y[tmp.order]
  u0.index = u0.index[tmp.order]
  
  # --------------------------------------------------
  if(!transposed) {
    matXtemp = bartModelMatrix(matX, matXnumcut, usequants=usequants,
                               xinfo=xinfo, rm.const=matXrm.const)
    matX = t(matXtemp$X)
    matXnumcut = matXtemp$numcut
    matXinfo = matXtemp$xinfo
    matXrm.const <- matXtemp$rm.const
    rm(matXtemp)
    
    matMtemp = bartModelMatrix(matM, matMnumcut, usequants=usequants,
                               xinfo=xinfo, rm.const=matMrm.const)
    matM = t(matMtemp$X)
    matMnumcut = matMtemp$numcut
    matMinfo = matMtemp$xinfo
    matMrm.const <- matMtemp$rm.const
    rm(matMtemp)
  } else {
    matXrm.const <- NULL
    matMrm.const <- NULL
  }
  
  if(n!=ncol(matX)){
    stop('The length of M and the number of rows in matX must be identical')
  }
  
  pm = nrow(matX)
  if(length(matXrho)==0) matXrho=pm
  if(length(matXrm.const)==0) matXrm.const <- 1:pm
  
  py = nrow(matM)
  if(length(matMrho)==0) matMrho=py
  if(length(matMrm.const)==0) matMrm.const <- 1:py
  
  uM=double(J);uM[1]=NaN
  uY=double(J);uY[1]=NaN
  # --------------------------------------------------
  #prior
  nu = sigdf
  
  if(typeM == "continuous") {
    if(is.na(Mlambda)) {
      if(is.na(Msigest)) {
        if(pm < n) {
          lmeMtemp = lme(M~., random=~1|factor(Uindex), data.frame(t(matX),Uindex,M))
          Msigest = summary(lmeMtemp)$sigma
          uM = c(lmeMtemp$coefficients$random[[1]])
          if(length(B_uM)==0) {
            B_uM=2*sd(uM)
          }
        } else {
          Msigest = sd(M)
        }
      }
      qchi = qchisq(1.0-sigquant,nu)
      Mlambda = (Msigest*Msigest*qchi)/nu # Mlambda parameter for sigma prior
    } else {
      Msigest=sqrt(Mlambda)
    }
  } else {
    Mlambda = 1
    Msigest = 1
  }
  if(is.na(Mtau.num)) {
    Mtau=(max(M)-min(M))/(2*k*sqrt(ntree))
  } else {
    Mtau = Mtau.num/sqrt(ntree)
  }
  
  if(typeY == "continuous") {
    if(is.na(Ylambda)) {
      if(is.na(Ysigest)) {
        if(py < n) {
          lmeYtemp = lme(Y~., random=~1|factor(Uindex), data.frame(t(matM),Uindex,Y))
          Ysigest = summary(lmeYtemp)$sigma
          uY = c(lmeYtemp$coefficients$random[[1]])
          if(length(B_uY)==0) {
            B_uY=2*sd(uY)
          }
        } else {
          Ysigest = sd(Y)
        }
      }
      qchi = qchisq(1.0-sigquant,nu)
      Ylambda = (Ysigest*Ysigest*qchi)/nu # Mlambda parameter for sigma prior
    } else {
      Ysigest=sqrt(Ylambda)
    }
  } else {
    Ylambda = 1
    Ysigest = 1
  }
  if(is.na(Ytau.num)) {
    Ytau=(max(Y)-min(Y))/(2*k*sqrt(ntree))
  } else {
    Ytau = Ytau.num/sqrt(ntree)
  }
  
  # --------------------------------------------------
  if(length(B_uM)==0) {
    B_uM <- Msigest
  }
  if(length(B_uY)==0) {
    B_uY <- Ysigest
  }
  
  # --------------------------------------------------
  ptm <- proc.time()
  res = .Call("crBARTmediation",
              ntypeM,   # 1:continuous, 2:binary, 3:multinomial
              ntypeY,   # 1:continuous, 2:binary, 3:multinomial
              n,        # number of observations in training data
              pm,       # dimension of matX
              matX,     # pm x n training data matX
              M,        # 1 x n training data M
              py,       # dimension of matM
              matM,     # py x n training data matM
              Y,        # 1 x n training data Y
              u0.index,
              n.j.vec,
              uM,       # random effects for M, if estimated
              uY,       # random effects for Y, if estimated
              J,
              B_uM,
              B_uY,
              ntree,
              matXnumcut,
              matMnumcut,
              ndpost*keepevery,
              nskip,
              keepevery,
              power,
              base,
              Moffset,
              Yoffset,
              Mtau,
              Ytau,
              nu,
              Mlambda,
              Ylambda,
              Msigest,
              Ysigest,
              sparse,
              theta,
              omega,
              a,
              b,
              matXrho,
              matMrho,
              augment,
              printevery,
              matXinfo,
              matMinfo)
  res$proc.time <- proc.time()-ptm
  # --------------------------------------------------
  res$Mdraw.mean <- apply(res$Mdraw, 2, mean)
  res$Mdraw.mean <- apply(res$Ydraw, 2, mean)
  
  names(res$matXtreedraws$cutpoints) = dimnames(matX)[[1]]
  dimnames(res$matXvarcount)[[2]] = as.list(dimnames(matX)[[1]])
  dimnames(res$matXvarprob)[[2]] = as.list(dimnames(matX)[[1]])
  res$matXvarcount.mean <- apply(res$matXvarcount, 2, mean)
  res$matXvarprob.mean <- apply(res$matXvarprob, 2, mean)
  
  names(res$matMtreedraws$cutpoints) = dimnames(matM)[[1]]
  dimnames(res$matMvarcount)[[2]] = as.list(dimnames(matM)[[1]])
  dimnames(res$matMvarprob)[[2]] = as.list(dimnames(matM)[[1]])
  res$matMvarcount.mean <- apply(res$matMvarcount, 2, mean)
  res$matMvarprob.mean <- apply(res$matMvarprob, 2, mean)
  
  res$typeM = typeM
  res$typeY = typeY
  
  res$Moffset = Moffset
  res$Yoffset = Yoffset
  res$matXrm.const <- matXrm.const
  res$matMrm.const <- matMrm.const
  res$Msigest <- Msigest
  res$Ysigest <- Ysigest
  res$B_uM <- B_uM
  res$B_uY <- B_uY
  
  res$uM <- uM
  res$uY <- uY
  
  if(.Platform$OS.type!='unix') {
    hostname <- FALSE
  } else if(hostname) {
    hostname <- system('hostname', intern=TRUE)
  }
  res$hostname <- hostname
  
  attr(res, 'class') <- "rBARTmediation"
  
  return(res)
}


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

mc.rBART <- function(Y, matX, Uindex=NULL, 
                     typeY = "continuous",
                     B_u=NULL,
                     sparse=FALSE, theta=0, omega=1,
                     a=0.5, b=1, augment=FALSE, rho=NULL,
                     xinfo=matrix(0,0,0), usequants=FALSE,
                     rm.const=TRUE,
                     sigest=NA, sigdf=3, sigquant=0.90,
                     k=2, power=2, base=0.95,
                     lambda=NA, tau.num=NA,
                     offsetY=mean(Y), 
                     ntree=200L, numcut=100L,
                     ndpost=1e3, nskip=1e4, keepevery=1e1,
                     printevery = (ndpost*keepevery)/10,
                     transposed=FALSE, hostname=FALSE,
                     mc.cores = 1L, nice = 19L, seed = 99L){
  
  #--------------------------------------------------
  # data  
  if(length(Uindex)==0)
    stop("the random effects indices must be provided")
  
  if(.Platform$OS.type!='unix')
    stop('parallel::mcparallel/mccollect do not exist on windows')
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  parallel::mc.reset.stream()
  
  if(!transposed) {
    temp = bartModelMatrix(matX, numcut, usequants=usequants,
                           xinfo=xinfo, rm.const=rm.const)
    matX = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    
    rm.const <- temp$rm.const
    rm(temp)
  }
  
  mc.cores.detected <- detectCores()
  
  if(mc.cores>mc.cores.detected) {
    mc.cores <- mc.cores.detected
  }
  
  mc.ndpost <- ceiling(ndpost/mc.cores)
  
  for(i in 1:mc.cores) {
    parallel::mcparallel({psnice(value=nice);
      rBART(Y=Y, matX=matX, Uindex=Uindex, 
            typeY = typeY,
            B_u=B_u,
            sparse=sparse, theta=theta, omega=omega,
            a=a, b=b, augment=augment, rho=rho,
            xinfo=xinfo, usequants=usequants,
            rm.const=rm.const,
            sigest=sigest, sigdf=sigdf, sigquant=sigquant,
            k=k, power=power, base=base,
            lambda=lambda, tau.num=tau.num,
            offsetY=offsetY, 
            ntree=ntree, numcut=numcut,
            ndpost=mc.ndpost, nskip=nskip,
            keepevery=keepevery, printevery=printevery,
            transposed=TRUE, hostname=hostname)},
      silent=(i!=1))
    ## to avoid duplication of output
    ## capture stdout from first posterior only
  }
  
  post.list <- parallel::mccollect()
  
  post <- post.list[[1]]
  
  if(mc.cores==1 | attr(post, 'class')!="rBART") {
    return(post)
  } else {
    if(class(rm.const)[1]!='logical') {
      post$rm.const <- rm.const
    }
    
    post$ndpost <- mc.cores*mc.ndpost
    
    p <- nrow(matX[post$rm.const, ])
    
    old.text <- paste0(as.character(mc.ndpost), ' ', as.character(ntree), ' ', as.character(p))
    old.stop <- nchar(old.text)
    
    post$treedraws$trees <- sub(old.text,
                                paste0(as.character(post$ndpost), ' ',
                                       as.character(ntree), ' ',
                                       as.character(p)),
                                post$treedraws$trees)
    
    for(i in 2:mc.cores) {
      post$hostname[i] <- post.list[[i]]$hostname
      post$yhat <- rbind(post$yhat, post.list[[i]]$yhat)
      post$sigma <- cbind(post$sigma, post.list[[i]]$sigma)
      post$varcount <- rbind(post$varcount, post.list[[i]]$varcount)
      post$varprob <- rbind(post$varprob, post.list[[i]]$varprob)
      post$sd.u <- cbind(post$sd.u, post.list[[i]]$sd.u)
      post$treedraws$trees <- paste0(post$treedraws$trees,
                                     substr(post.list[[i]]$treedraws$trees, old.stop+2,
                                            nchar(post.list[[i]]$treedraws$trees)))
      
      post$proc.time['elapsed'] <- max(post$proc.time['elapsed'],
                                       post.list[[i]]$proc.time['elapsed'])
      for(j in 1:5)
        if(j!=3)
          post$proc.time[j] <- post$proc.time[j]+post.list[[i]]$proc.time[j]
    }
    
    post$an <- apply(post$yhat, 1, function(x) !any(is.nan(x)))
    post$yhat.mean <- apply(post$yhat[post$an, ], 2, mean)
    post$varcount.mean <- apply(post$varcount, 2, mean)
    post$varprob.mean <- apply(post$varprob, 2, mean)
    
    if(all(post$an)) post$an <- NULL ## no NAs
    
    attr(post, 'class') <- "rBART"
    
    return(post)
  }
}
################################################################################
##
##   Copyright (C) 2015 - 2016 Alfred Galichon
##
##   This file is part of the R package TraME.
##
##   The R package TraME is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 2 of the License, or
##   (at your option) any later version.
##
##   The R package TraME is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################
#
################################################################################
########################       default methods           #######################
################################################################################
initparam.default <- function(model)
{
  ret = list(param=rep(0,model$nbParams),
             lb = NULL, ub = NULL)
}

estimate.default = mle
#
################################################################################
########################      ETU-logit model            #######################
################################################################################
#
buildModel_ETU_logit <- function(Xvals, Yvals, n=NULL, m=NULL)
{
  nbX = dim(t(t(Xvals)))[1]
  nbY = dim(t(t(Yvals)))[1]
  #
  dX = dim(t(t(Xvals)))[2]
  dY = dim(t(t(Yvals)))[2]
  #
  eX = matrix(rep(1,nbX),ncol=1)
  eY = matrix(rep(1,nbY),ncol=1)
  #
  diff = abs(kronecker(eY,Xvals)-kronecker(Yvals,eX))
  #
  if(is.null(n)){
    n=rep(1,nbX)
  }
  if(is.null(m)){
    m=rep(1,nbY)
  }
  #
  ret = list(types = c("itu-rum", "mfe"),
             diff=diff,
             nbParams=2*dim(t(t(diff)))[2]+1,
             nbX=nbX, nbY=nbY,
             dX=dX, dY=dY,
             n=n, m=m,
             isLinear=TRUE)
  class(ret) = "ETU_logit"
  #
  return(ret)
}

parametricMarket.etu <- function(model, theta)
  # the theta are the parameters for alpha, gamma and tau
{
  theta1 = theta[1:model$dX]
  theta2 = theta[(model$dX+1):(model$dX+model$dY)]
  theta3 = theta[length(theta)]
  #
  alpha = matrix(model$diff %*% theta1,nrow=model$nbX)
  gamma = matrix(model$diff %*% theta2,nrow=model$nbX)
  
  tau = matrix(theta3, model$nbX, model$nbY)
  #
  ret = build_market_ETU_logit(model$n,model$m,alpha,gamma,tau,sigma=1)
  #
  return(ret)
}

dparam.etu <- function(model, dparams=diag(model$nbParams))
  # params is simply the affinity matrix
{
  zero1 = matrix(0,model$nbX*model$nbY,model$dX)
  zero2 = matrix(0,model$nbX*model$nbY,model$dY)
  zero3 = matrix(0,model$nbX*model$nbY,1)
  #
  dparamsPsi = rbind(cbind(model$diff,zero2,zero3),
                     cbind(zero1,model$diff,zero3),
                     cbind(zero1,zero2,rep(1,model$nbX*model$nbY)))
  dparamsG = matrix(0,nrow=1,ncol=dim(dparams)[2])
  dparamsH = matrix(0,nrow=1,ncol=dim(dparams)[2])
  #
  ret = list(dparamsPsi = dparamsPsi,
             dparamsG   = dparamsG,
             dparamsH   = dparamsH)
  #
  return(ret)
}

initparam.etu <- function(model)
{
  ret = list(param=c(rep(0,model$nbParams-1),1),
             lb=NULL,ub=NULL)
  #
  return(ret)
}
#
################################################################################
########################      TU-empirical model            ####################
################################################################################
#
buildModel_TU_empirical = function(n,m,phi_xyk, hetG, hetH, nbDraws,seed=777) {
  dims = dim(phi_xyk)
  nbParams = dims[1]
  nbX = dims[2]
  nbY = dims[3]
  ret = list(  phi_xyk = phi_xyk,
               nbParams = nbParams,                         
               nbX=nbX,
               nbY=nbY,
               n = n,
               m = m,
               hetG=hetG,
               hetH=hetH)
  
  class(ret) =   "TU_empirical"
  return(ret)
  
}
#
parametricMarket.TU_empirical <- function(model, theta)
{
  phimat = matrix(model$phi_xyk,ncol = model$nbParams)
  Phimat = apply(phimat,1,sum)
  Phi = matrix(Phi,model$nbX,model$nbY)
return(  build_market_TU_general(model$n,model$m,Phi,model$hetG,model$hetH))
}
#
dparam.TU_empirical  <- function(model,dparams=diag(model$nbParams))
  # params is (lambda1,lambda2)
{
  dparamsPsi = matrix(model$phi_xyk,ncol = model$nbParams) %*% dparams
  dparamsG = matrix(0,nrow=0,ncol=dim(dparams)[2])
  dparamsH = matrix(0,nrow=0,ncol=dim(dparams)[2])
  return( list(dparamsPsi=dparamsPsi,
               dparamsG = dparamsG,
               dparamsH = dparamsH))
}


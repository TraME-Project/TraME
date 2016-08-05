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
########################       affinity model            #######################
################################################################################
# The TU_logit and the affinity models should be merged 
#
buildModel_affinity <- function(Xvals, Yvals, n=NULL, m=NULL, noSingles=FALSE)
{
  nbX = dim(Xvals)[1]
  nbY = dim(Yvals)[1]
  #
  dX = dim(Xvals)[2]
  dY = dim(Yvals)[2]
  #
  if(is.null(n)){
    n = rep(1,nbX)
  }
  if(is.null(m)){
    m = rep(1,nbY)
  }
  #
  neededNorm = defaultNorm(noSingles)
  #
  ret = list(types = c("itu-rum", "mfe"),
             kron=kronecker(Yvals,Xvals),
             nbParams=dX*dY,
             dX=dX, dY=dY,
             nbX = nbX, nbY = nbY,
             n=n, m=m,
             neededNorm = neededNorm,
             isLinear=TRUE)
  class(ret) = "affinity"
  #
  return(ret)
}
#
parametricMarket.affinity <- function(model, theta)
  # theta is simply the affinity matrix
{
  phi = matrix(model$kron %*% c(theta),nrow=model$nbX)
  #
  ret = build_market_TU_logit(model$n,model$m,phi,
                              neededNorm=model$neededNorm)
  #
  return(ret)
}


dparam.affinity <- function(model, dparams=diag(model$nbParams))
{
  dparamsPsi = model$kron %*% dparams
  dparamsG = matrix(0,nrow=0,ncol=dim(dparams)[2])
  dparamsH = matrix(0,nrow=0,ncol=dim(dparams)[2])
  #
  ret = list(dparamsPsi = dparamsPsi,
             dparamsG   = dparamsG,
             dparamsH   = dparamsH)
  #
  return(ret)
}
#
#
mme.affinity <- function(model, muhat, xtol_rel=1e-4, maxeval=1e5, print_level=0)
  # mme.affinity should be improved as one should make use of the logit structure and use the ipfp
{
  if (print_level>0){
    message(paste0("Moment Matching Estimation of TU_rum model via optimization."))
  }
  
  theta0 = initparam(model)$param
  dtheta = dparam(model)
  market = parametricMarket(model,theta0)
  
  kron = dtheta$dparamsPsi
  Chat = c(c(muhat) %*% kron)
  #
  if (print_level>0){
    message(paste0("Moment Matching Estimation of ",class(model)," model via BFGS optimization."))
  }
  #
  nbX = length(model$n)
  nbY = length(model$m)
  
  nbParams = dim(kron)[2]
  #
  eval_f <- function(thearg){
    theU = matrix(thearg[1:(nbX*nbY)],nbX,nbY)
    thetheta = thearg[(1+nbX*nbY):(nbParams+nbX*nbY)]
    
    phi = kron %*% thetheta
    phimat = matrix(phi,nbX,nbY)
    #
    resG = G(market$arumsG,theU,model$n)
    resH = G(market$arumsH,t(phimat-theU),model$m)
    #
    Ehatphi = sum(thetheta * Chat)
    val = resG$val + resH$val - Ehatphi
    
    tresHmu = t(resH$mu)
    
    gradU = c(resG$mu - tresHmu)
    gradtheta = c( c(tresHmu) %*% kron ) - Chat
    #
    ret = list(objective = val,
               gradient = c(gradU,gradtheta))
    #
    return(ret)
  }
  #
  resopt = nloptr(x0=c( (kron %*% theta0) / 2,theta0 ),
                  eval_f=eval_f,
                  opt=list("algorithm" = "NLOPT_LD_LBFGS",
                           "xtol_rel"=xtol_rel,
                           "maxeval"=maxeval,
                           "print_level"=print_level))
  #
  #if(print_level>0){print(resopt, show.controls=((1+nbX*nbY):(nbParams+nbX*nbY)))}
  #
  U = matrix(resopt$solution[1:(nbX*nbY)],nbX,nbY)  
  thetahat = resopt$solution[(1+nbX*nbY):(nbParams+nbX*nbY)]
  V = matrix(kron %*% thetahat,nbX,nbY) - U
  #
  ret =list(thetahat=thetahat,
            U=U, V=V,
            val=resopt$objective)
  #
  return(ret)
}
#
################################################################################
########################       TU_logit model            #######################
################################################################################
# The TU_logit and the affinity models should be merged
buildModel_TU_logit <- function(phi_xyk, n=NULL, m=NULL,noSingles=FALSE)
{
  dims = dim(phi_xyk)
  nbX = dims[1]
  nbY = dims[2]
  nbParams = dims[3]
  #
  if(is.null(n)){
    n = rep(1,nbX)
  }
  if(is.null(m)){
    m = rep(1,nbY)
  }
  #  
  neededNorm = defaultNorm(noSingles)
  #
  ret = list(types = c("itu-rum", "mfe"),
             phi_xyk = phi_xyk,
             nbParams = nbParams,
             nbX = nbX, nbY = nbY,
             n=n, m=m,
             neededNorm = neededNorm,
             isLinear=TRUE)
  class(ret) = "TU_logit"
  #
  return(ret)
}
#
parametricMarket.TU_logit<- function(model, theta)
  # theta is the parameter vector for phi
{
  phi_xy_vec = matrix(model$phi_xyk,ncol = model$nbParams) %*% theta
  phi_xy_mat = matrix(phi_xy_vec,model$nbX,model$nbY)
  return( build_market_TU_logit(model$n,model$m,phi_xy_mat,
                                neededNorm=model$neededNorm) )
}

dparam.TU_logit <- function(model, dparams=diag(model$nbParams))
{
  dparamsPsi = matrix(model$phi_xyk,ncol = model$nbParams) %*% dparams
  dparamsG = matrix(0,nrow=0,ncol=dim(dparams)[2])
  dparamsH = matrix(0,nrow=0,ncol=dim(dparams)[2])
  return( list(dparamsPsi=dparamsPsi,
               dparamsG = dparamsG,
               dparamsH = dparamsH))
}
#
mme.TU_logit <- mme.affinity 
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
buildModel_TU_empirical = function(phi_xyk, n=NULL, m=NULL, arumsG, arumsH) {
  if  (class(arumsG)!="empirical" )
  {stop("arumsG provided to buildModel_TU_empirical is not of class empirical.")}
  if  (class(arumsH)!="empirical" )
  {stop("arumsH provided to buildModel_TU_empirical is not of class empirical.")}
  dims = dim(phi_xyk)
  nbX = dims[1]
  nbY = dims[2]
  #
  if(is.null(n)){
    n = rep(1,nbX)
  }
  if(is.null(m)){
    m = rep(1,nbY)
  }
  #
  nbParams = dims[3]
  ret = list(  phi_xyk = phi_xyk,
               nbParams = nbParams,                         
               nbX=nbX,
               nbY=nbY,
               n = n,
               m = m,
               arumsG=arumsG,
               arumsH=arumsH)
  
  class(ret) =   "TU_empirical"
  return(ret)
  
}
#
parametricMarket.TU_empirical <- function(model, theta)
{
  phi_xy_vec = matrix(model$phi_xyk,ncol = model$nbParams) %*% theta
  phi_xy_mat = matrix(phi_xy_vec,model$nbX,model$nbY)
return(  build_market_TU_general(model$n,model$m,phi_xy_mat,model$arumsG,model$arumsH))
}
#
dparam.TU_empirical  <- function(model,dparams=diag(model$nbParams))
{
  dparamsPsi = matrix(model$phi_xyk,ncol = model$nbParams) %*% dparams
  dparamsG = matrix(0,nrow=0,ncol=dim(dparams)[2])
  dparamsH = matrix(0,nrow=0,ncol=dim(dparams)[2])
  return( list(dparamsPsi=dparamsPsi,
               dparamsG = dparamsG,
               dparamsH = dparamsH))
}
#
mme.TU_empirical <- function(model, muhat, xtol_rel=1e-4, maxeval=1e5, print_level=0)
{
  if (print_level>0){
    message(paste0("Moment Matching Estimation of TU_empirical model via LP optimization."))
  }
  kron = matrix(model$phi_xyk,ncol = model$nbParams)
  Chat = c(c(muhat) %*% kron)
  #
  nbX = length (model$n)
  nbY = length (model$m)
  #
  res1 = build_disaggregate_epsilon(model$n,nbX,nbY,model$arumsG)
  res2 = build_disaggregate_epsilon(model$m,nbY,nbX,model$arumsH)
  #
  epsilon_iy = res1$epsilon_iy
  epsilon0_i = c(res1$epsilon0_i)
  I_ix = res1$I_ix
  #
  eta_xj = t(res2$epsilon_iy)
  eta0_j = c(res2$epsilon0_i)  
  I_yj = t(res2$I_ix)
  #
  ni = c(I_ix %*% model$n)/res1$nbDraws
  mj = c( model$m %*% I_yj)/res2$nbDraws
  #
  nbI = length(ni)
  nbJ = length(mj)
  #
  # based on this, can compute aggregated equilibrium in LP 
  #
  A_11 = kronecker(matrix(1,nbY,1),sparseMatrix(1:nbI,1:nbI,x=1))
  A_12 = sparseMatrix(i=NULL,j=NULL,dims=c(nbI*nbY,nbJ),x=0)
  A_13 = kronecker(sparseMatrix(1:nbY,1:nbY,x=-1),I_ix)
  A_14 = sparseMatrix(i=NULL,j=NULL,dims=c(nbI*nbY,model$nbParams),x=0)
  #
  A_21 = sparseMatrix(i=NULL,j=NULL,dims=c(nbX*nbJ,nbI),x=0)
  A_22 = kronecker(sparseMatrix(1:nbJ,1:nbJ,x=1),matrix(1,nbX,1))
  A_23 = kronecker(t(I_yj),sparseMatrix(1:nbX,1:nbX,x=1))
  A_24 = -t(matrix(matrix(t(kron),model$nbParams*nbX,nbY) %*% I_yj, model$nbParams, nbX*nbJ))
  #
  A_1  = cbind(A_11,A_12,A_13, A_14)
  A_2  = cbind(A_21,A_22,A_23, A_24)
  #
  A    = rbind(A_1,A_2)
  #
  nbconstr = dim(A)[1]
  nbvar = dim(A)[2]
  #
  lb  = c(epsilon0_i,t(eta0_j), rep(-Inf,nbX*nbY+model$nbParams))
  rhs = c(epsilon_iy, eta_xj)
  obj = c(ni,mj,rep(0,nbX*nbY),c(-Chat))
  #
  result = genericLP(obj=obj,A=A,modelsense="min",rhs=rhs,sense=rep(">=",nbconstr),lb=lb)
  #
  U = matrix(result$solution[(nbI+nbJ+1):(nbI+nbJ+nbX*nbY)],nrow=nbX)
  thetahat = result$solution[(nbI+nbJ+nbX*nbY+1):(nbI+nbJ+nbX*nbY+model$nbParams)]
  V = matrix(kron %*% thetahat,nbX,nbY) - U
  #
  muiy = matrix(result$pi[1:(nbI*nbY)],nrow=nbI)
  mu = t(I_ix) %*% muiy
  #
  val = result$objval
  #
  ret = list(thetahat=thetahat,
             U=U, V=V,
             val=val)
  #
  return(ret)
}
#
################################################################################
########################          TU-none model             ####################
################################################################################
#
buildModel_TU_none = function(phi_xyk, n=NULL, m=NULL,seed=777) {
  dims = dim(phi_xyk)
  nbX = dims[1]
  nbY = dims[2]
  nbParams = dims[3]
  ret = list(  phi_xyk = phi_xyk,
               nbParams = nbParams,                         
               nbX=nbX,
               nbY=nbY,
               n = n,
               m = m)
  class(ret) =   "TU_none"
  return(ret)
}
#
parametricMarket.TU_none<- function(model, theta)
{
  phi_xy_vec = matrix(model$phi_xyk,ncol = model$nbParams) %*% theta
  phi_xy_mat = matrix(phi_xy_vec,model$nbX,model$nbY)
  return( build_market_TU_none(model$n,model$m,phi_xy_mat) )
}
#
dparam.TU_none  <- function(model,dparams=diag(model$nbParams))
{
  dparamsPsi = matrix(model$phi_xyk,ncol = model$nbParams) %*% dparams
  dparamsG = matrix(0,nrow=0,ncol=dim(dparams)[2])
  dparamsH = matrix(0,nrow=0,ncol=dim(dparams)[2])
  return( list(dparamsPsi=dparamsPsi,
               dparamsG = dparamsG,
               dparamsH = dparamsH))
}
#
mme.TU_none <- function(model, muhat, xtol_rel=1e-4, maxeval=1e5, print_level=0)
# MomentMatchingTUNone <- function(n, m, kron, Chat, print_level=0)
{
  if (print_level>0){
    message(paste0("Moment Matching Estimation of TU_none model via LP optimization."))
  }
  kron = matrix(model$phi_xyk,ncol = model$nbParams)
  Chat = c(c(muhat) %*% kron)
  #
  nbX = length (model$n)
  nbY = length (model$m)
  #
  A_1 = kronecker(matrix(1,nbY,1),sparseMatrix(1:nbX,1:nbX))
  A_2 = kronecker(sparseMatrix(1:nbY,1:nbY),matrix(1,nbX,1))
  A_3 = -kron
  #
  A   = cbind(A_1,A_2,A_3)
  #
  nbconstr = dim(A)[1]
  nbvar = dim(A)[2]
  #
  rhs = rep(0,nbX*nbY)
  obj = c(model$n,model$m,c(-Chat))
  lb =c(rep(0,nbX+nbY),rep(-Inf,model$nbParams))
  #
  result = genericLP(obj=obj,A=A,modelsense="min",rhs=rhs,sense=rep(">=",nbconstr),lb=lb)
  #
  u = result$solution[1:nbX]
  v = result$solution[(nbX+1):(nbX+nbY)]
  thetahat = result$solution[(1+nbX+nbY):(model$nbParams+nbX+nbY)]
  mu = matrix(result$pi,nbX,nbY)
  val = result$objval
  #
  ret = list(thetahat=thetahat,
             u=u, v=v,
             val=val)
  #
  return(ret)
}
#
################################################################################
########################            TU-rum model            ####################
################################################################################
#
buildModel_TU_rum = function(phi_xyk, n=NULL, m=NULL, arumsG, arumsH) {
  dims = dim(phi_xyk)
  nbX = dims[1]
  nbY = dims[2]
  #
  if(is.null(n)){
    n = rep(1,nbX)
  }
  if(is.null(m)){
    m = rep(1,nbY)
  }
  #
  nbParams = dims[3]
  ret = list(  phi_xyk = phi_xyk,
               nbParams = nbParams,                         
               nbX=nbX,
               nbY=nbY,
               n = n,
               m = m,
               arumsG=arumsG,
               arumsH=arumsH)
  
  class(ret) =   "TU_rum"
  return(ret)
  
}
#
parametricMarket.TU_rum <- function(model, theta)
{
  phi_xy_vec = matrix(model$phi_xyk,ncol = model$nbParams) %*% theta
  phi_xy_mat = matrix(phi_xy_vec,model$nbX,model$nbY)
  return(  build_market_TU_general(model$n,model$m,phi_xy_mat,model$arumsG,model$arumsH))
}
#
dparam.TU_rum  <- function(model,dparams=diag(model$nbParams))
{
  dparamsPsi = matrix(model$phi_xyk,ncol = model$nbParams) %*% dparams
  dparamsG = matrix(0,nrow=0,ncol=dim(dparams)[2])
  dparamsH = matrix(0,nrow=0,ncol=dim(dparams)[2])
  return( list(dparamsPsi=dparamsPsi,
               dparamsG = dparamsG,
               dparamsH = dparamsH))
}
#
mme.TU_rum <- function(model, muhat, xtol_rel=1e-4, maxeval=1e5, print_level=0)
{
  if (print_level>0){
    message(paste0("Moment Matching Estimation of TU_rum model via optimization."))
  }

  theta0 = initparam(model)$param
  dtheta = dparam(model)
  market = parametricMarket(model,theta0)
  
  kron = dtheta$dparamsPsi
  Chat = c(c(muhat) %*% kron)
  #
  if (print_level>0){
    message(paste0("Moment Matching Estimation of ",class(model)," model via BFGS optimization."))
  }
  #
  nbX = length(model$n)
  nbY = length(model$m)
  
  nbParams = dim(kron)[2]
  #
  eval_f <- function(thearg){
    theU = matrix(thearg[1:(nbX*nbY)],nbX,nbY)
    thetheta = thearg[(1+nbX*nbY):(nbParams+nbX*nbY)]
    
    phi = kron %*% thetheta
    phimat = matrix(phi,nbX,nbY)
    #
    resG = G(market$arumsG,theU,model$n)
    resH = G(market$arumsH,t(phimat-theU),model$m)
    #
    Ehatphi = sum(thetheta * Chat)
    val = resG$val + resH$val - Ehatphi
    
    tresHmu = t(resH$mu)
    
    gradU = c(resG$mu - tresHmu)
    gradtheta = c( c(tresHmu) %*% kron ) - Chat
    #
    ret = list(objective = val,
               gradient = c(gradU,gradtheta))
    #
    return(ret)
  }
  #
  resopt = nloptr(x0=c( (kron %*% theta0) / 2,theta0 ),
                  eval_f=eval_f,
                  opt=list("algorithm" = "NLOPT_LD_LBFGS",
                           "xtol_rel"=xtol_rel,
                           "maxeval"=maxeval,
                           "print_level"=print_level))
  #
  #if(print_level>0){print(resopt, show.controls=((1+nbX*nbY):(nbParams+nbX*nbY)))}
  #
  U = matrix(resopt$solution[1:(nbX*nbY)],nbX,nbY)  
  thetahat = resopt$solution[(1+nbX*nbY):(nbParams+nbX*nbY)]
  V = matrix(kron %*% thetahat,nbX,nbY) - U
  #
  ret =list(thetahat=thetahat,
            U=U, V=V,
            val=resopt$objective)
  #
  return(ret)
}


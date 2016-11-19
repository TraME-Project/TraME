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
#
buildModel_affinity <- function(Xvals, Yvals, n=NULL, m=NULL, sigma = 1 )
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
  if ( sum(n) != sum(m) ) {stop("Unequal mass of individuals in an affinity model.")}
  #
  neededNorm = defaultNorm(T) # Keith: why set this but then set 'neededNorm = F' in the list?
  #
  ret = list(types = c("itu-rum", "mfe"),
             nbParams=dX*dY,
             dX=dX, dY=dY,
             nbX = nbX, nbY = nbY,
             n=n, m=m,
             sigma = sigma,
             neededNorm = neededNorm,
             phi_xyk_aux = kronecker(Yvals,Xvals),
             Phi_xyk = function(model) 
               (model$phi_xyk_aux),
             Phi_xy = function(model ,lambda) 
               ( array(model$phi_xyk_aux,dim=c(model$nbX*model$nbY,model$nbParams)) %*% lambda ), # Keith: lambda is a matrix here but a scalar elsewhere?
             Phi_k = function(model, muhat) 
               (c(c(muhat) %*% array(model$phi_xyk_aux,dim=c(model$nbX*model$nbY,model$nbParams))))
               )
  class(ret) = "affinity"
  #
  return(ret)
}
#
parametricMarket.affinity <- function(model, theta) 
# Keith: shouldn't this have sigma from model$  and needed_norm? 
  (build_market_TU_logit(model$n,model$m,
                         matrix(model$Phi_xy(model,c(theta)), nrow=model$nbX),
                         neededNorm=F))

dparam.affinity <- function(model, dparams=diag(model$nbParams))
  (list(dparamsPsi = model$Phi_xy(model, dparams),
        dparamsG   = matrix(0,nrow=0,ncol=dim(dparams)[2]),
        dparamsH   = matrix(0,nrow=0,ncol=dim(dparams)[2]))
  )
#
# Keith: this function should be specific to the affinity class
mmeaffinityNoRegul <- function(model, muhat, xtol_rel=1e-4, maxeval=1e5, tolIpfp=1E-14, maxiterIpfp = 1e5, print_level=0)
  # mmeaffinityNoRegul should be improved as one should make use of the logit structure and use the ipfp
{
  #
  if (print_level>0){
    message(paste0("Moment Matching Estimation of Affinity model via IPFP+BFGS optimization."))
  }
  #
  theta0 = initparam(model)$param
  market = parametricMarket(model,theta0) # Keith: why assign this if it's never used?
  Chat = model$Phi_k(model, muhat)
  nbX = model$nbX
  nbY = model$nbY
  dX = model$dX
  dY = model$dY
  nbParams = model$nbParams
  sigma = model$sigma
  #
  totmass = sum(model$n)
  
  if ( sum(model$n) != totmass ) {stop("Unequal mass of individuals in an affinity model.")}
  if (sum(muhat) !=  totmass) {stop("Total number of couples does not coincide with margins.")}
  p = model$n / totmass
  q = model$m / totmass
  IX=rep(1,nbX)
  tIY=matrix(rep(1,nbY),nrow=1)
  f = p %*% tIY
  g = IX %*% t(q)
  pihat = muhat / totmass
  v=rep(0,nbY)
  #
  #
  valuef = function(A)
  {
    Phi = matrix(model$Phi_xy(model,c(A)) ,nbX,nbY)
    # Phi = Xvals %*% matrix(A,nrow=dX) %*% t(Yvals)
    contIpfp = TRUE
    iterIpfp = 0
    while(contIpfp)
    {
      iterIpfp = iterIpfp+1
      u = sigma*log(apply(g * exp( ( Phi - IX %*% t(v) ) / sigma ),1,sum))
      vnext = sigma*log(apply(f * exp( ( Phi - u %*% tIY ) / sigma ),2,sum))
      error = max(abs(apply(g * exp( ( Phi - IX %*% t(vnext) - u %*% tIY ) / sigma ),1,sum)-1))
      if( (error<tolIpfp) | (iterIpfp >= maxiterIpfp)) {contIpfp=FALSE}
      v=vnext
    }
    #print(c("Converged in ", iterIpfp, " iterations."))
    pi = f * g * exp( ( Phi - IX %*% t(v) - u %*% tIY ) / model$sigma ) # Keith: don't use names like 'pi' in a function; use 'Pi' instead
    if (iterIpfp >= maxiterIpfp ) {stop('maximum number of iterations reached')} 
    v <<- vnext # Keith: Global assignment is not allowed in the package!
    #thegrad =  c(    c(pi - pihat) %*% phis)
    thegrad = model$Phi_k(model,pi - pihat)
    theval = sum(thegrad * c(A)) - sigma* sum(pi*log(pi))
    
    return(list(objective = theval,gradient = thegrad))
  }
  
  A0 = rep(0,dX*dY)
  resopt = nloptr(x0=A0, 
                  eval_f = valuef, 
                  opt = list(algorithm = 'NLOPT_LD_LBFGS',
                             xtol_rel = xtol_rel,
                             maxeval=maxeval,
                             print_level = print_level))
  #  AffinityMatrix = matrix(res$solution,nrow=dX)
  thetahat = resopt$solution
  ret =list(thetahat=thetahat,
            val=resopt$objective)
  #
  return(ret)
}
#
#
mmeaffinityWithRegul <- function(model, muhat, lambda, xtol_rel=1e-4, maxeval=1e5, tolIpfp=1E-14, maxiterIpfp = 1e5, print_level=0)
  # Reference: Arnaud Dupuy, Alfred Galichon, Yifei Sun (2016). "Learning Optimal Transport Costs under Low-Rank Constraints."
  # Implementation by Yifei Sun.
{
  #
  if (print_level>0){
    message(paste0("Moment Matching Estimation of Affinity model with regularization via proximal gradient."))
  }
  #
  theta0 = initparam(model)$param
  market = parametricMarket(model,theta0) # Keith: why assign this if it's never used?
  Chat = model$Phi_k(model, muhat)
  nbX = model$nbX
  nbY = model$nbY
  dX = model$dX
  dY = model$dY
  nbParams = model$nbParams
  sigma = model$sigma
  #
  totmass = sum(model$n)
  #
  if ( sum(model$n) != totmass ) {stop("Unequal mass of individuals in an affinity model.")}
  if (sum(muhat) !=  totmass) {stop("Total number of couples does not conicide with margins.")}
  p = model$n / totmass
  q = model$m / totmass
  IX=rep(1,nbX)
  tIY=matrix(rep(1,nbY),nrow=1)
  f = p %*% tIY
  g = IX %*% t(q)
  pihat = muhat / totmass
  v=rep(0,nbY)
  A = rep(0,dX*dY)
  t_k = .3   # step size for the prox grad algorithm (or grad descent when lambda=0)
  iterCount = 0
  while (1)
  {
    # compute pi^A
    # Phi = Xvals %*% matrix(A,nrow=dX) %*% t(Yvals)
    Phi = matrix(model$Phi_xy(model,c(A)) ,nbX,nbY)
    contIpfp = TRUE
    iterIpfp = 0
    while(contIpfp)
    {
      iterIpfp = iterIpfp+1
      u = sigma*log(apply(g * exp( ( Phi - IX %*% t(v) ) / sigma ),1,sum))
      vnext = sigma*log(apply(f * exp( ( Phi - u %*% tIY ) / sigma ),2,sum))
      error = max(abs(apply(g * exp( ( Phi - IX %*% t(vnext) - u %*% tIY ) / sigma ),1,sum)-1))
      if( (error<tolIpfp) | (iterIpfp >= maxiterIpfp)) {contIpfp=FALSE}
      v=vnext
    }
    
    pi = f * g * exp( ( Phi - IX %*% t(v) - u %*% tIY ) / sigma )
    
    if (iterIpfp >= maxiterIpfp ) {stop('maximum number of iterations reached')} 
    
    # do prox grad descent
    # thegrad = c(phis %*% c(pi - pihat))
    thegrad = model$Phi_k(model,pi - pihat)
    
    
    # take one gradient step
    A = A - t_k*thegrad
    
    if (lambda > 0)
    {
      # compute the proximal operator
      SVD = svd(matrix(A,nrow=dX))
      U = SVD$u
      D = SVD$d
      V = SVD$v
      
      D = pmax(D - lambda*t_k, 0.) # Keith: typo: 0.?
      A = c(U %*% diag(D) %*% t(V))
    } # if lambda = 0 then we are just taking one step of gradient descent
    ### testing optimality
    if (iterCount %% 10 == 0)
    {
      alpha = 1. # Keith: typo: 1.? why is this not set outside the loop?
      tmp = svd(matrix(A - alpha * thegrad, nrow=dX))
      tmp_second = sum((A - c(tmp$u %*% diag(pmax(tmp$d - alpha*lambda, 0.)) %*% t(tmp$v)))^2)
      cat("testing optimality ", tmp_second, "\n")
    }
    
    if (lambda > 0)
    {
      theval = sum(thegrad * c(A)) - sigma * sum(pi*log(pi)) + lambda * sum(D)
    } else
    {
      theval = sum(thegrad * c(A)) - sigma * sum(pi*log(pi))
    }
    
    iterCount = iterCount + 1
    
    if (iterCount>1 && abs(theval - theval_old) < 1E-6) { break }
    
    theval_old = theval
    
  }
  ret =list(thetahat=c(A),
            val=theval)
  
  #
  return(ret)
}
#
#
mme.affinity <- function(model, muhat, lambda = NULL, xtol_rel=1e-4, maxeval=1e5, tolIpfp=1E-14, maxiterIpfp = 1e5, print_level=0)
{
  if (is.null(lambda))
  {return(mmeaffinityNoRegul(model,muhat, xtol_rel , maxeval , tolIpfp , maxiterIpfp , print_level))}
  else
  { return(mmeaffinityWithRegul(model,muhat, lambda, xtol_rel , maxeval , tolIpfp , maxiterIpfp , print_level)) }
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
             neededNorm = neededNorm)
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
mme.TU_logit <-  function(model, muhat, xtol_rel=1e-4, maxeval=1e5, print_level=0)
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
#
################################################################################
########################      ETU-logit model            #######################
################################################################################
#
buildModel_ETU_logit <- function(Xvals, Yvals, n=NULL, m=NULL)
{
  nbX = dim(t(t(Xvals)))[1] # Keith: double transpose?
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
             n=n, m=m)
  class(ret) = "ETU_logit"
  #
  return(ret)
}
# Keith: something is wrong here; class should be ETU_logit instead of etu? 
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
mme.TU_none <- function(model, muhat, method = 1, xtol_rel=1e-4, maxeval=1e5, print_level=0)
{
  if (method==1)
{return(MARP_proj(model,muhat,xtol_rel ,maxeval,print_level ))}
  else
  {return(MARP_min(model,muhat,xtol_rel ,maxeval,print_level ))}

  }

  
MARP_min <-function(model, muhat, xtol_rel=1e-4, maxeval=1e5, print_level=0)   
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
  A_11 = kronecker(matrix(1,nbY,1),sparseMatrix(1:nbX,1:nbX))
  A_12 = kronecker(sparseMatrix(1:nbY,1:nbY),matrix(1,nbX,1))
  A_13 = -kron
  A_1   = cbind(A_11,A_12,A_13)
  A_2 = matrix(c(rep(0,nbX+nbY),rep(1,model$nbParams)),nrow=1)
  #
  A = rbind(A_1,A_2)
  #
  nbconstr = dim(A)[1]
  nbvar = dim(A)[2]
  #
  rhs = c(rep(0,nbX*nbY),1)
  obj = c(model$n,model$m,c(-Chat))
  result = genericLP(obj=obj,A=A,modelsense="min",rhs=rhs,sense=c(rep(">",nbconstr-1),"="),lb=c(rep(0,nbX+nbY),rep(-Inf, model$nbParams) ))
  u = result$solution[1:nbX]
  v = result$solution[(nbX+1):(nbX+nbY)]
  thetahat = result$solution[(1+nbX+nbY):(model$nbParams+nbX+nbY)]
  mu = matrix(result$pi[1:(nbX*nbY)],nbX,nbY)
  val = result$objval
  #
  ret = list(thetahat=thetahat,
             u=u, v=v,
             val=val)
  #
  return(ret)
}


MARP_proj <- function(model, muhat, xtol_rel=1e-4, maxeval=1e5, print_level=0)
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
  A_11 = kronecker(matrix(1,nbY,1),sparseMatrix(1:nbX,1:nbX))
  A_12 = kronecker(sparseMatrix(1:nbY,1:nbY),matrix(1,nbX,1))
  A_13 = -kron
  A_1   = cbind(A_11,A_12,A_13)
  A_2 = matrix(c(rep(0,nbX+nbY),c(Chat)),nrow=1)
  #
  A = rbind(A_1,A_2)
  #
  nbconstr = dim(A)[1]
  nbvar = dim(A)[2]
  #
  rhs = c(rep(0,nbX*nbY),1)
  obj = c(model$n,model$m,rep(0,model$nbParams))
  result = genericLP(obj=obj,A=A,modelsense="min",rhs=rhs,sense=c(rep(">=",nbconstr-1),"="),lb=c(rep(0,nbX+nbY),rep(-Inf, model$nbParams) ))
  val = result$objval
  u = result$solution[1:nbX]
  v = result$solution[(nbX+1):(nbX+nbY)]
  thetahat =result$solution[(1+nbX+nbY):(model$nbParams+nbX+nbY)]
  gamma = result$pi[nbX*nbY+1]
  #
  ret = list(thetahat=thetahat,
             u=u, v=v,
             gamma =gamma,
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


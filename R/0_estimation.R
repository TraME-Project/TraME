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

dtheta_mu_default <- function(model, market, theta, dtheta=diag(length(theta)))
{
    outcome = solveEquilibrium(market,notifications=FALSE)
    #
    mu = outcome$mu
    mux0s = market$n-apply(mu,1,sum)
    mu0ys = market$m-apply(mu,2,sum)
    
    U = outcome$U
    V = outcome$V
    
    dthetaPsiGH = dparam(model, dtheta)
    dthetaPsi = dthetaPsiGH$dparamsPsi
    dthetaG = dthetaPsiGH$dparamsG
    dthetaH = dthetaPsiGH$dparamsH
    
    tr = market$transfers
    hetG = market$hetG
    hetH = market$hetH
    
    duPsimat = du_Psi(tr,U,V)
    dvPsimat = 1 - duPsimat
    duPsivec = c(duPsimat)
    dvPsivec = c(dvPsimat)
    
    #   duPsi = Diagonal(x = duPsimat) 
    #   dvPsi = Diagonal( x = dvPsimat)
    #
    HessGstar = D2Gstar(market$hetG,mu,market$n,xFirst=T)  
    HessHstar = D2Gstar(market$hetH,t(mu),market$m,xFirst=F)
    #
    denom = duPsivec * HessGstar + dvPsivec * HessHstar
    
    num1 = dtheta_Psi(tr,U,V,dthetaPsi)
    num2 = duPsivec * dtheta_NablaGstar(hetG,mu,market$n,dthetaG,xFirst=TRUE)
    num3 = dvPsivec * dtheta_NablaGstar(hetH,t(mu),market$m,dthetaH,xFirst=FALSE)
    #
    dmu = -solve(denom, num1)
    num = cbind(num1, num2, num3)
    #
    return(list(mu= c(mu),mux0s=mux0s, mu0ys=mu0ys,dmu=dmu))
}

dtheta_mu_logit <- function(model, market, theta, dtheta=diag(length(theta)))
{
    rangeParams = dim(dtheta)[2]
    sigma = market$hetG$sigma
    
    dthetaPsiGH = dparam(model, dtheta)
    dthetaPsi = dthetaPsiGH$dparamsPsi
    
    tr = market$transfers
    #
    outcome = solveEquilibrium(market,notifications=FALSE,debugmode=FALSE)
    #
    mu = outcome$mu
    mux0s = outcome$mux0
    mu0ys = outcome$mu0y
    #
    us = matrix(outcome$u,nrow=tr$nbX,ncol=tr$nbY)
    vs = matrix(outcome$v,nrow=tr$nbX,ncol=tr$nbY,byrow=T)
    
    du_psis = du_Psi(tr,us,vs)
    dv_psis = 1 - du_psis
    #
    dtheta_psis = matrix(dtheta_Psi(tr,us,vs,dthetaPsi),nrow=tr$nbX*tr$nbY)
    mudthetapsi = array(c(mu)*c(dtheta_psis),dim=c(tr$nbX,tr$nbY,rangeParams))
    
    d_1 = apply(mudthetapsi, c(1,3), sum) / sigma
    d_2 = apply(mudthetapsi, c(2,3), sum) / sigma
    num = rbind(d_1,d_2)
    #
    Delta11 = diag(mux0s + apply(mu*du_psis,1,sum),nrow=tr$nbX)
    Delta22 = diag(mu0ys + apply(mu*dv_psis,2,sum),nrow=tr$nbY)
    Delta12 = mu * dv_psis
    Delta21 = t(mu * du_psis)
    Delta = rbind(cbind(Delta11,Delta12),cbind(Delta21,Delta22))
    #
    dlogmusingles = solve(Delta,num)
    dlogmux0 = dlogmusingles[1:tr$nbX,,drop=FALSE]
    dlogmu0y = dlogmusingles[(tr$nbX+1):(tr$nbX+tr$nbY),,drop=FALSE]
    dlogmux0full = array(0,dim=c(tr$nbX,tr$nbY,rangeParams))
    dlogmu0yfull = array(0,dim=c(tr$nbX,tr$nbY,rangeParams))
    #
    for(y in 1:tr$nbY){
        dlogmux0full[,y,] = dlogmux0
    }
    for(x in 1:tr$nbX){
        dlogmu0yfull[x,,] = dlogmu0y
    }
    #
    dlogmu = c(du_psis)*matrix(dlogmux0full, ncol=rangeParams) + 
             c(dv_psis)*matrix(dlogmu0yfull, ncol=rangeParams) - 
             matrix(dtheta_psis,ncol=rangeParams) / sigma    
    dmu    = c(mu) * dlogmu
    #
    return(list(mu = c(mu), mux0s = mux0s, mu0ys = mu0ys, dmu = dmu))
}

dtheta_mu_numeric <- function (model, market, theta, dtheta=diag(length(theta)))
{ 
    thef <- function(theparams){
        ret = solveEquilibrium(parametricMarket(model,theparams),notifications=FALSE)$mu
        return(ret)
    }
    #
    outcome = solveEquilibrium( parametricMarket(model,theta),notifications=FALSE)
    #
    mu=outcome$mu
    mux0s = outcome$mux0
    mu0ys = outcome$mu0y
    #
    dmu = jacobian(thef,theta) %*% dtheta
    #
    return(list(mu = c(mu), mux0s = mux0s, mu0ys = mu0ys, dmu = dmu))
}

dtheta_mu <- function(model, theta, dtheta=diag(length(theta)))
{
    market = parametricMarket(model,theta)
    #
    ret <- 0
    check_1 = (class(market$hetG)=="logit")
    check_2 = (class(market$hetH)=="logit")
    check_3 = (market$hetG$sigma == market$hetH$sigma)
    
    if(check_1 && check_2 && check_3){
        ret = dtheta_mu_logit(model,market,theta,dtheta)
    }else{
        ret = dtheta_mu_default(model,market,theta,dtheta)
    }
    #
    return(ret)
}

mLogLikelihood <- function(theta,model,muhat,muhatx0,muhat0y, scale=1) UseMethod("mLogLikelihood", model)

mLogLikelihood.default <- function(theta, model, muhat, muhatx0, muhat0y,scale=1)
{
    mudmu = try( dtheta_mu(model,theta),silent=T)
    #
    ret <- 0
    if(class(mudmu)!="try-error"){ 
        mLL = - sum(2* muhat * log(mudmu$mu)) - sum(muhatx0 * log(mudmu$mux0s)) - sum(muhat0y * log(mudmu$mu0ys))
        #
        term_1 = t(2*muhat/matrix(mudmu$mu,nrow=model$nbX) - muhatx0/mudmu$mux0s)
        term_2 = muhat0y / mudmu$mu0ys
        term_grad = c(t(term_1 - term_2))*mudmu$dmu
        #
        mGradLL = - apply(term_grad,2,sum) 
        #
        ret = list(objective=mLL / scale ,gradient=mGradLL / scale)
    }else{
        ret = list(objective=NaN, gradient=rep(NaN,model$nbParams))
    }
    #
    return(ret)
}

mle <- function(model, muhat, theta0=NULL, xtol_rel=1e-8, maxeval=1e5, print_level=0)
{
  nbX = length(model$n)
  nbY = length(model$m)
  scale = max(sum(model$n),sum(model$n))
  nbParams = length(model$nbParams)
    if(print_level > 0){
        message(paste0("Maximum Likelihood Estimation of ",class(model)," model."))
    }
    #
    muhatx0  = model$n-apply(muhat,1,sum)
    muhat0y  = model$m-apply(muhat,2,sum)
    #
    if(is.null(theta0)){
        theta0 = initparam(model)$param
    }
    #
    lb     = initparam(model)$lb
    ub     = initparam(model)$ub
    #
    res = nloptr(x0=theta0, 
                 eval_f=mLogLikelihood,
                 lb=lb, ub=ub,
                 opt = list(algorithm='NLOPT_LD_LBFGS',
                            xtol_rel=xtol_rel,
                            maxeval=maxeval,
                            "print_level"=print_level), 
                 model=model,
                 muhat=muhat,
                 muhatx0=muhatx0,
                 muhat0y=muhat0y,
                 scale = scale)
    #
    if(print_level > 0){
        print(res, show.controls=((1+nbX*nbY):(nbParams+nbX*nbY)))
    }
    #
    return(list(thetahat=res$solution))
}

MomentMatchingTUSmooth <- function(n, m, hetG, hetH, kron, Chat, theta0, xtol_rel=1e-4, maxeval=1e5, print_level=0)
{
    if(print_level>0){
        message("BFGS optimization used.")
    }
    #
    nbX = length(n)
    nbY = length(m)
    
    nbParams = dim(kron)[2]
    #
    eval_f <- function(thearg){
        theU = matrix(thearg[1:(nbX*nbY)],nbX,nbY)
        thetheta = thearg[(1+nbX*nbY):(nbParams+nbX*nbY)]
        
        phi = kron %*% thetheta
        phimat = matrix(phi,nbX,nbY)
        #
        resG = G(hetG,theU,n)
        resH = G(hetH,t(phimat-theU),m)
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

MomentMatchingTUEmpirical <- function(n, m, hetG, hetH, kron, Chat, theta0, xtol_rel=1e-4, maxeval=1e5, print_level=0)
{
    if (print_level>0){
        message("LP optimization used.")
    }
    #
    nbX = length (n)
    nbY = length (m)
    nbParams = length(Chat)
    #
    res1 = build_disaggregate_epsilon(n,nbX,nbY,hetG)
    res2 = build_disaggregate_epsilon(m,nbY,nbX,hetH)
    #
    epsilon_iy = res1$epsilon_iy
    epsilon0_i = c(res1$epsilon0_i)
    I_ix = res1$I_ix
    
    eta_xj = t(res2$epsilon_iy)
    eta0_j = c(res2$epsilon0_i)  
    I_yj = t(res2$I_ix)
    #
    ni = c(I_ix %*% n)/res1$nbDraws
    mj = c( m %*% I_yj)/res2$nbDraws
    
    nbI = length(ni)
    nbJ = length(mj)
    #
    # based on this, can compute aggregated equilibrium in LP 
    #
    A_11 = kronecker(matrix(1,nbY,1),sparseMatrix(1:nbI,1:nbI,x=1))
    A_12 = sparseMatrix(i=NULL,j=NULL,dims=c(nbI*nbY,nbJ),x=0)
    A_13 = kronecker(sparseMatrix(1:nbY,1:nbY,x=-1),I_ix)
    A_14 = sparseMatrix(i=NULL,j=NULL,dims=c(nbI*nbY,nbParams),x=0)
    
    A_21 = sparseMatrix(i=NULL,j=NULL,dims=c(nbX*nbJ,nbI),x=0)
    A_22 = kronecker(sparseMatrix(1:nbJ,1:nbJ,x=1),matrix(1,nbX,1))
    A_23 = kronecker(t(I_yj),sparseMatrix(1:nbX,1:nbX,x=1))
    A_24 = -t(matrix(matrix(t(kron),nbParams*nbX,nbY) %*% I_yj, nbParams, nbX*nbJ))
    
    A_1  = cbind(A_11,A_12,A_13, A_14)
    A_2  = cbind(A_21,A_22,A_23, A_24)
    
    A    = rbind(A_1,A_2)
    #
    nbconstr = dim(A)[1]
    nbvar = dim(A)[2]
    #
    lb  = c(epsilon0_i,t(eta0_j), rep(-Inf,nbX*nbY+nbParams))
    rhs = c(epsilon_iy, eta_xj)
    obj = c(ni,mj,rep(0,nbX*nbY),c(-Chat))
    #
    result = genericLP(obj=obj,A=A,modelsense="min",rhs=rhs,sense=rep(">=",nbconstr),lb=lb)
    #
    U = matrix(result$solution[(nbI+nbJ+1):(nbI+nbJ+nbX*nbY)],nrow=nbX)
    thetahat = result$solution[(nbI+nbJ+nbX*nbY+1):(nbI+nbJ+nbX*nbY+nbParams)]
    V = matrix(kron %*% thetahat,nbX,nbY) - U
    
    muiy = matrix(result$pi[1:(nbI*nbY)],nrow=nbI)
    mu = t(I_ix) %*% muiy
    
    val = result$objval
    #
    ret = list(thetahat=thetahat,
               U=U, V=V,
               val=val)
    #
    return(ret)
    
}

MomentMatchingTUNone <- function(n, m, kron, Chat, print_level=0)
{
    if(print_level > 0){
        message("LP optimization used.")
    }
    #
    nbX = length (n)
    nbY = length (m)
    nbParams = length(Chat)
    #
    A_1 = kronecker(matrix(1,nbY,1),sparseMatrix(1:nbX,1:nbX))
    A_2 = kronecker(sparseMatrix(1:nbY,1:nbY),matrix(1,nbX,1))
    A_3 = -kron
    
    A   = cbind(A_1,A_2,A_3)
    #
    nbconstr = dim(A)[1]
    nbvar = dim(A)[2]
    #
    rhs = rep(0,nbX*nbY)
    obj = c(n,m,c(-Chat))
    lb =c(rep(0,nbX+nbY),rep(-Inf,nbParams))
    #
    result = genericLP(obj=obj,A=A,modelsense="min",rhs=rhs,sense=rep(">=",nbconstr),lb=lb)
    #
    u = result$solution[1:nbX]
    v = result$solution[(nbX+1):(nbX+nbY)]
    thetahat = result$solution[(1+nbX+nbY):(nbParams+nbX+nbY)]
    mu = matrix(result$pi,nbX,nbY)
    val = result$objval
    #
    ret = list(thetahat=thetahat,
               u=u, v=v,
               val=val)
    #
    return(ret)
}

mme <- function(model, muhat, print_level=0)
{
    if(print_level>0){
        message(paste0("Moment Matching Estimation of ",class(model)," model."))
    }
    #
    theta0 = initparam(model)$param
    #  lb    =initparam(model)$lb
    #  ub    =initparam(model)$ub
    dtheta = dparam(model)
    market = parametricMarket(model,theta0)
    #
    if(class(market$transfers)!="TU"){
        stop("Moment Matching Estimation only applies for TU models.")
    }
    if(length(dtheta$dparamsG) + length(dtheta$dparamsH) > 0){
        stop("Moment Matching Estimation does not support parameterization of heterogeneity.")
    }
    if(!model$isLinear){
        stop("Moment Matching Estimation does not support nonlinear parameterizations.")
    }
    #
    kron = dtheta$dparamsPsi
    Chat = c(c(muhat) %*% kron)
    #
    if((class(market$hetG)=="empirical") & (class(market$hetH)=="empirical")){
        outcome = MomentMatchingTUEmpirical(market$n,market$m,market$hetG,market$hetH,
                                            kron,Chat,theta0,print_level=print_level)
    }else if((class(market$hetG)=="none") & (class(market$hetH)=="none")){
        outcome = MomentMatchingTUNone(market$n,market$m,kron,Chat,print_level=print_level)
    }else{
        outcome = MomentMatchingTUSmooth(market$n,market$m,market$hetG,market$hetH,
                                         kron,Chat,theta0,print_level=print_level)
    }
    #
    return(outcome)
}
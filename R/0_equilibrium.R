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

# References:
# A. Galichon, B. Salanie: " Cupid's Invisible Hand: Social Surplus and Identification in Matching Models"
# A. Galichon, Y.-W. Hsieh: "Love and Chance: Equilibrium and Identification in a Large NTU matching markets with stochastic choice"
# A. Galichon, S.D. Kominers, and S. Weber: "An Empirical Framework for Matching with Imperfectly Transferable Utility"  
# O. Bonnet, A. Galichon, and M. Shum: "Yoghurt Chooses Man: The Matching Approach to Identification of Nonadditive Random Utility Models".
#

ipfp <- function(market, xFirst=T, notifications=TRUE, debugmode=FALSE, tol=1e-12, mu0ystart=market$m)
    # Computes equilibrium in the logit case via IPFP in the all-logit case
{
    #
  if (!("mmfs" %in% market$types)) {stop("ipfp only defined for market whose types contains mmfs.")}
  mmfs = market$mmfs
    noSingles = !is.null(mmfs$neededNorm)
    if(noSingles){
        warning("There are known issues with the current implementation of the IPFP in the case without unassigned agents.")
        H = mmfs$neededNorm$H_edge_logit
        if(is.null(H)){
            stop("Function H_edge not included in market$neededNorm")
        }
    }
    #
    if(notifications){
        message('Solving for equilibrium in ITU_logit problem using IPFP.') 
    }
    n = mmfs$n
    m = mmfs$m
    nbX = length(n)
    nbY = length(m)
    #
    # Algorithm: Loop
    #
    mu0y = mu0ystart
    mux0 = rep(NA,length=nbX)
    #
    error = 2*tol
    iter = 0
    tm = proc.time()  
    while(max(error,na.rm=TRUE)>tol){
        iter = iter+1
        val = c(mux0,mu0y)
        
        #Solve for mux0 and then mu0y
        mux0 = margxInv(1:nbX,mmfs=mmfs,Mu0y=mu0y)
        mu0y = margyInv(1:nbY,mmfs=mmfs,Mux0=mux0)
        
        if(noSingles){
            rescale = H(mux0,mu0y)
            mux0 = mux0 * rescale
            mu0y = mu0y / rescale
        }
        
        error = abs(c(mux0,mu0y)-val)
        if(debugmode & notifications){
            message(paste0("Iter: ", iter, ". Error: ", max(error)))
        }
    }
    #
    time = proc.time()-tm  
    time = time["elapsed"] 
    if(notifications){
        message(paste0("IPFP converged in ", iter," iterations and ", round(time,digits=2), " seconds.\n"))
    }
    # Construct the equilibrium outcome based on mux0 and mu0y obtained from above
    mu = M(mmfs,mux0,mu0y)  
    #
    outcome = list(mu = mu,
                   mux0 = mux0, mu0y = mu0y,
                   iter=iter, time=time)
    #
    return(outcome)
}


nodalNewton <- function(market, xFirst=TRUE, notifications=FALSE, sigma = 1E-6, maxiter = 100, tol = 1e-3, xtol = 1e-3)
{
  if(!is.null(market$neededNorm)){
    stop("nodalNewton does not yet allow for the case without unmatched agents.")
  }
  if(notifications){
    message('Solving for equilibrium in MFE problem using nodalNewton.')
  }
  #
  if((class(market$arumsG)!="logit") || (class(market$arumsH)!="logit")){
    stop("Error: Heterogeneities are not all logit.")
  }
  if(market$arumsG$sigma != market$arumsH$sigma){
    stop("Error: scaling parameters of the logit differs on both sides of the market.")
  }
  #
  n = market$n
  m = market$m
  sigma = market$arumsG$sigma
  tr = market$transfers
  #
  theus = rep(0,tr$nbX)
  thevs = rep(0,tr$nbY)
  themux0s = rep(0,tr$nbX)
  themu0ys = rep(0,tr$nbY)
  themu = matrix(0,tr$nbX,tr$nbY)
  #
  Z <- function(uv)
  {
    theus <<- uv[1:tr$nbX]
    thevs <<- uv[(tr$nbX+1):(tr$nbX+tr$nbY)]
    themux0s <<- exp(- theus / sigma)
    themu0ys <<- exp(- thevs / sigma)
    themu <<- MMF(tr,themux0s,themu0ys,sigma=sigma)
    return( c(themux0s + apply(themu,1,sum) - n,
              themu0ys + apply(themu,2,sum) - m ))
  }
  #
  JZ <- function(uv)
  {
    #     theus = uv[1:tr$nbX]
    #     thevs = uv[(tr$nbX+1):(tr$nbX+tr$nbY)] ## note -- this is redundant with the computation of Z
    #     themux0s = exp(- theus / sigma) ## note -- this is redundant with the computation of Z
    #     themu0ys = exp(- thevs / sigma)  ## note -- this is redundant with the computation of Z
    #     themu = MMF(tr,themux0s,themu0ys,sigma=sigma) ## note -- this is redundant with the computation of Z
    du_psis = du_Psi(tr,theus,thevs)
    dv_psis = 1 - du_psis
    Delta11 = diag(themux0s + apply(themu*du_psis,1,sum),nrow=tr$nbX)
    Delta22 = diag(themu0ys + apply(themu*dv_psis,2,sum),nrow=tr$nbY)
    Delta12 = themu * dv_psis
    Delta21 = t(themu * du_psis)
    Delta = rbind(cbind(Delta11,Delta12),cbind(Delta21,Delta22))
    hess = (-1/sigma) * Delta
    return(hess)
  }
  #
  xinit = -sigma* c(log(n/2),log(m/2))
  tm = proc.time()  
  sol = nleqslv(x = xinit,
                fn = Z, jac = JZ,
                method = "Broyden", # "Newton"
                control = list(xtol=xtol,maxit = maxiter)
  )
  time = proc.time()-tm  
  time = time["elapsed"] 
  iter = sol$iter
  if(notifications){
    message(paste0("nodalNewton converged in ", iter," iterations and ", round(time,digits=2), " seconds.\n"))
  }
  #
  theus = sol$x[1:tr$nbX]
  thevs = sol$x[(tr$nbX+1):(tr$nbX+tr$nbY)]
  themux0s = exp(- theus / sigma)
  themu0ys = exp(- thevs / sigma)
  themu = MMF(tr,themux0s,themu0ys,sigma=sigma)
  error = max(abs(c(themux0s + apply(themu,1,sum) - n,themu0ys + apply(themu,2,sum) - m ))) # calling this gives the right value to themu, themux0s and themu0ys
  U = sigma * log(themu/themux0s)
  V = sigma * t(log(t(themu) / themu0ys))
  #
  outcome = list(mu = themu,
                 mux0 = themux0s, mu0y = themu0ys,
                 U = U, V = V,
                 u = theus,
                 v = thevs,
                 iter=iter, time=time,
                 error=error)
  #
  return(outcome)
}



arcNewton <- function(market, xFirst=TRUE, notifications=TRUE, wup=NULL, xtol=1e-5, method="Broyden")
    # method is "Broyden" or "Newton"
{
    if(method!="Broyden" && method!="Newton"){
        stop("invalid input passed to function 'newton': unknown method selected.\nMethod must be either \"Broyden\" or \"Newton\".")
    }
    if(!is.null(market$neededNorm)){
        stop("Newton does not allow for the case without unmatched agents")
    }
    if(notifications){
        message('Solving for equilibrium in ITU_general problem using Newton descent.\n')
    }
    #
    n = market$n
    m = market$m
    
    arumsG = market$arumsG
    arumsH = market$arumsH
    tr = market$transfers
    
    nbX = length(n)
    nbY = length(m)
    #
    if(is.null(wup)){
        winit =  wupperbound(market)
    }else{
        winit = wup
        if(min(G(arumsG,UW(tr,wup),n)$mu - t(G(arumsH,t(VW(tr,wup)),m)$mu)) < 0){
            stop("wup provided is not an actual upper bound")
        }
    }
    #
    ED <- function(thew){
        term_1 = G(arumsG,UW(tr,matrix(thew,nbX,nbY)),n)$mu
        term_2 = t(G(arumsH,t(VW(tr,matrix(thew,nbX,nbY))),m)$mu)
        #
        ret = c(term_1 - term_2)
        #
        return(ret)
    }
    
    jacED <- function(thew){
        term_1 = c(dw_UW(tr,matrix(thew,nbX,nbY)))*D2G(arumsG,UW(tr,matrix(thew,nbX,nbY)),n)
        term_2 = c(dw_VW(tr,matrix(thew,nbX,nbY)))*D2G(arumsH,t(VW(tr,matrix(thew,nbX,nbY))),m,xFirst=F)
        #
        ret = term_1 - term_2
        #
        return(ret)
    }
    #
    sol = nleqslv(x = c(winit),
                  fn = ED, jac = jacED,
                  method = "Broyden", #"Newton",
                  control = list(xtol=xtol,maxit = 200)
    )
    #
    w = sol$x
    #
    U = UW(tr,w)
    V = VW(tr,w)
    
    mu = G(arumsG,U,n)$mu
    
    mux0 = n - apply(mu,1,sum)
    mu0y = m - apply(mu,2,sum)
    #
    ret = list(mu=mu,
               mux0=mux0, mu0y=mu0y,
               U=U, V=V, sol=sol)
    #
    return(ret)
}

wupperbound <- function(market)
{
    n = market$n
    m = market$m
    
    arumsG = market$arumsG
    arumsH = market$arumsH
    tr = market$transfers
    
    nbX = length(n)
    nbY = length(m)
    
    w = matrix(0,nbX,nbY)
    U = matrix(0,nbX,nbY)
    V = matrix(0,nbX,nbY)
    #
    k=1
    cont = TRUE
    while(cont==TRUE){
        for(x in 1:nbX){
            typeTransfersx = determineType(tr,x)
            #
            if(typeTransfersx==1){
                mu_condx = rep( 1 / (2^(-k)+nbY),nbY)
                
                U[x,] = Gstarx(arumsG,mu_condx,x )$Ux
                w[x,] = WU(tr,U[x,],xs=x)
                V[x,] = VW(tr,w[x,],xs=x)
            }else{ # if transfers = type2
                if(typeTransfersx==2){
                    w[x,] = rep(2^k,nbY)
                    U[x,] = UW(tr,w[x,],xs=x)
                    V[x,] = VW(tr,w[x,],xs=x)
                }else{
                    stop("Multiple transfer types is not allowed.")
                }
            }
        }
        #
        Z = G(arumsG,U,n)$mu - t(G(arumsH,t(V),m)$mu)
        if(min(Z) < 0){
            k = 2*k
        }else{
            cont = FALSE
        }
    } #end while
    #
    return(w)
}

jacobi <- function(market, xFirst=TRUE, notifications=TRUE, wlow=NULL, wup=NULL, tol=1e-10)
{
    if(!is.null(market$neededNorm)){
        stop("Jacobi does not yet allow for the case without unmatched agents.")
    }
    if(notifications){
        message('Solving for equilibrium in ITU_general problem using Jacobi iterations.')
    }
    #
    n = market$n
    m = market$m
    
    arumsG = market$arumsG
    arumsH = market$arumsH
    tr = market$transfers
    
    nbX = length(n)
    nbY = length(m)
    #
    if(is.null(wup)){
        w = wupperbound(market)
    }else{
        w = wup
        Z = G(arumsG,UW(tr,wup),n)$mu - t(G(arumsH,t(VW(tr,wup)),m)$mu)
        if(min(Z) < 0 ){
            stop("wup provided not an actual upper bound")
        }
    }
    #
    if(is.null(wlow)){
        wlow = -t(wupperbound(marketTranspose(market)))
    }
    else{
        wlow = wlow
        Z = G(arumsG,UW(tr,wlow),n)$mu - t(G(arumsH,t(VW(tr,wlow)),m)$mu)
        
        if(max(Z)>0){
            stop("wlow provided not an actual lower bound")
        }
    }
    #
    U = UW(tr,w)
    V = VW(tr,w)
    Z = G(arumsG,U,n)$mu - t(G(arumsH,t(V),m)$mu)
    #
    iter=0
    while(max(abs(Z / (n %*% t(m))))>1e-4){
        #
        iter=iter+1
        #
        for(x in 1:nbX){
            for(y in 1:nbY){
                tofit <- function(thew)
                {
                    U[x,y] = UW(tr,thew,x,y)
                    V[x,y] = VW(tr,thew,x,y)
                    ed = G(arumsG,U,n)$mu - t(G(arumsH,t(V),m)$mu)
                    return(ed[x,y])
                }
                #
                w[x,y] = uniroot(tofit,c(wlow[x,y],w[x,y]),tol = tol, extendInt="upX")$root
                U[x,y] = UW(tr,w[x,y],x,y)
                V[x,y] = VW(tr,w[x,y],x,y)
            }
        }
        Z = G(arumsG,U,n)$mu - t(G(arumsH,t(V),m)$mu)
    }
    #
    if(notifications){
        message(paste0("Jacobi iteration converged in ", iter," iterations.\n"))
    }
    #
    mu = G(arumsG,U,n)$mu  
    mux0 = n - apply(mu,1,sum)
    mu0y = m - apply(mu,2,sum)
    #
    ret = list(mu=mu,
               mux0=mux0, mu0y=mu0y,
               U=U, V=V)
    #
    return(ret)
}

maxWelfare <- function(market, xFirst=TRUE, notifications=FALSE, tol_rel=1e-8)
{
    if(!is.null(market$neededNorm)){
        stop("maxWelfare does not yet allow for the case without unmatched agents.")
    }
    if(class(market$transfers)!="TU"){
        stop("maxWelfare only works for TU transfers.")
    }
    if(notifications){
        message('Solving for equilibrium in TU_general problem using convex optimization.\n')
    }
    #
    nbX = length(market$n)
    nbY = length(market$m)
    
    phi = market$transfers$phi
    #
    eval_f <- function(theU)
    {
        theU = matrix(theU,nbX,nbY)
        resG = G(market$arumsG,theU,market$n)
        resH = G(market$arumsH,t(phi-theU),market$m)
        val  = resG$val + resH$val
        #
        ret = list(objective = val,
                   gradient = c(resG$mu - t(resH$mu)))
        #
        return(ret)
    }
    #
    U_init = phi / 2
    
    resopt = nloptr(x0 = U_init, eval_f = eval_f,
                    opt = list("algorithm" = "NLOPT_LD_LBFGS",
                               "xtol_rel"=tol_rel,
                               "ftol_rel"=1e-15))
    #
    U = matrix(resopt$solution,nbX,nbY)
    resG = G(market$arumsG,U,market$n)
    resH = G(market$arumsH,t(phi-U),market$m)
    mu = resG$mu
    V = phi - U
    #
    mux0 = market$n - apply(mu,1,sum)
    mu0y = market$m - apply(mu,2,sum)
    #
    ret = list(mu=mu,
               mux0=mux0, mu0y=mu0y,
               U=U, V=V,
               val = resG$val + resH$val)
    #
    return(ret)
}

darum <- function(market, xFirst=TRUE, notifications=FALSE, tol=1e-8)
{
    if(!is.null(market$neededNorm)){
        stop("darum does not yet allow for the case without unmatched agents.")
    }
    if(class(market$transfers)!="NTU"){
        stop("darum only works for NTU transfers.")
    }
    if(notifications){
        message('Solving for equilibrium in NTU_general problem using DARUM.')
    }
    #
    alpha = market$transfers$alpha
    gamma = market$transfers$gamma
    
    n = market$n
    m = market$m
    
    arumsG = market$arumsG
    arumsH = market$arumsH
    
    nbX = length(n)
    nbY = length(m)
    #
    muNR = pmax(n %*% matrix(1,1,nbY), matrix(1,nbX,1) %*% t(m))
    #
    error = 2*tol
    iter = 0
    #
    while((max(error,na.rm=TRUE)>tol) & (iter<10000)){
        iter  = iter + 1
        #
        resP  = Gbar(arumsG,alpha,n,muNR)
        muP   = resP$mu
        #
        resD  = Gbar(arumsH,t(gamma),m,t(muP))
        muD   = t(resD$mu)
        #
        muNR  = muNR - (muP - muD)
        error = abs(muP- muD)
    }
    #
    if(notifications){
        message(paste0("Darum iteration converged in ", iter," iterations.\n"))
    }
    #
    mux0 = n-apply(muD,1,sum)
    mu0y = m-apply(muD,2,sum)
    #
    outcome = list(mu=muD,
                   mux0=mux0, mu0y=mu0y,
                   U = resP$U, V = t(resD$U))
    #
    return(outcome)
}

build_disaggregate_epsilon <- function(n, nbX, nbY, arums) 
    # takes U_xy and arums as input; returns U_iy as output
{
    nbDraws = arums$aux_nbDraws
    nbI = nbX*nbDraws
    
    epsilons = matrix(0,nbI,nbY+1)
    I_ix = matrix(0,nbI,nbX)
    #
    for(x in 1:nbX){
        if(arums$xHomogenous){
            epsilon = arums$atoms
        }else{
            epsilon = matrix(arums$atoms[,,x],nrow=arums$aux_nbDraw)
        }
        #
        epsilons[((x-1)*nbDraws+1):(x*nbDraws),] = epsilon
        #
        I_01 = ifelse((1:nbX)==x,1,0)
        I_ix[((x-1)*nbDraws+1):(x*nbDraws),] = matrix(rep(t(I_01),nbDraws),nrow=nbDraws,byrow=T)
    }
    epsilon_iy = epsilons[,1:nbY]
    epsilon0_i = epsilons[,nbY+1]
    #
    ret = list(epsilon_iy=epsilon_iy,
               epsilon0_i=epsilon0_i,
               I_ix=I_ix,
               nbDraws=nbDraws)
    #
    return(ret)  
}

CupidsLP <- function(market, xFirst=TRUE, notifications=FALSE)
{
    if(!is.null(market$neededNorm)){
        stop("CupidsLP does not yet allow for the case without unmatched agents.")
    }
    if((class(market$transfers)!="TU")||(class(market$arumsG)!="empirical")||(class(market$arumsH)!="empirical")){
        stop("cupidsLP only works for TU-empirical markets.")
    }
    if(notifications){
        message('Solving for equilibrium in TU_empirical problem using LP.\n')
    }
    #
    nbX = length (market$n)
    nbY = length (market$m)
    #
    phi = market$transfers$phi
    res1 = build_disaggregate_epsilon(market$n,nbX,nbY,market$arumsG)
    res2 = build_disaggregate_epsilon(market$m,nbY,nbX,market$arumsH)
    #
    epsilon_iy = res1$epsilon_iy
    epsilon0_i = c(res1$epsilon0_i)
    
    I_ix = res1$I_ix
    
    eta_xj = t(res2$epsilon_iy)
    eta0_j = c(res2$epsilon0_i)  
    
    I_yj = t(res2$I_ix)
    
    ni = c(I_ix %*% market$n)/res1$nbDraws
    mj = c(market$m %*% I_yj)/res2$nbDraws
    
    nbI = length(ni)
    nbJ = length(mj)
    #
    # based on this, can compute aggregated equilibrium in LP 
    #
    A_11 = suppressMessages( Matrix::kronecker(matrix(1,nbY,1),sparseMatrix(1:nbI,1:nbI,x=1)) )
    A_12 = sparseMatrix(i=NULL,j=NULL,dims=c(nbI*nbY,nbJ),x=0)
    A_13 = suppressMessages( Matrix::kronecker(sparseMatrix(1:nbY,1:nbY,x=-1),I_ix) )
    
    A_21 = sparseMatrix(i=NULL,j=NULL,dims=c(nbX*nbJ,nbI),x=0)
    A_22 = suppressMessages( Matrix::kronecker(sparseMatrix(1:nbJ,1:nbJ,x=1),matrix(1,nbX,1)) )
    A_23 = suppressMessages( Matrix::kronecker(t(I_yj),sparseMatrix(1:nbX,1:nbX,x=1)) )
    
    A_1  = cbind(A_11,A_12,A_13)
    A_2  = cbind(A_21,A_22,A_23)
    
    A    = rbind(A_1,A_2)
    # 
    nbconstr = dim(A)[1]
    nbvar = dim(A)[2]
    #
    lb  = c(epsilon0_i,t(eta0_j), rep(-Inf,nbX*nbY))
    rhs = c(epsilon_iy, eta_xj+phi %*% I_yj)
    obj = c(ni,mj,rep(0,nbX*nbY))
    sense = rep(">=",nbconstr)
    modelsense = "min"
    #
    result = genericLP(obj=obj,A=A,modelsense=modelsense,rhs=rhs,sense=sense,lb=lb)
    #
    U = matrix(result$solution[(nbI+nbJ+1):(nbI+nbJ+nbX*nbY)],nrow=nbX)
    V = phi - U
    
    muiy = matrix(result$pi[1:(nbI*nbY)],nrow=nbI)
    mu = t(I_ix) %*% muiy
    
    val = sum(ni*result$solution[1:nbI]) + sum(mj*result$solution[(nbI+1):(nbI+nbJ)])
    #
    ret = list(success=TRUE,
               mu=mu,
               mux0 = market$n - apply(mu,1,sum),
               mu0y = market$m - apply(mu,2,sum),
               U=U, V=V,
               val=result$objval)
    #
    return(ret)
}

oapLP <- function(market, xFirst=TRUE, notifications=FALSE)
{
    if(!is.null(market$neededNorm)){
        stop("oapLP does not yet allow for the case without unmatched agents.")
    }
    if(class(market$transfers)!="TU"){
        stop("oapLP only works for TU transfers.")
    }
    if(notifications){
        message('Solving for equilibrium in TU_none problem using Linear Programming.\n')
    }
    #
    phi = market$transfers$phi
    nbX = dim(phi)[1]
    nbY = dim(phi)[2]
    #
    obj = c(phi)
    
    A1 = Matrix::kronecker(matrix(1,1,nbY),sparseMatrix(1:nbX,1:nbX))
    A2 = Matrix::kronecker(sparseMatrix(1:nbY,1:nbY),matrix(1,1,nbX))
    A = rbind2(A1,A2)
    
    d = c(market$n,market$m)
    pi_init = c(market$n %*% t(market$m))
    #
    result = genericLP(obj=obj,A=A,modelsense="max",rhs=d,sense="<",start=pi_init)
    #
    mu = matrix(result$solution,nrow=nbX)
    u0 = result$pi[1:nbX] 
    v0 = result$pi[(nbX+1):(nbX+nbY)]
    val = result$objval
    #
    objBis = ifelse(xFirst==TRUE,1,-1)*c(market$n,-market$m)
    ABis = rbind2(Matrix::t(A),c(market$n,market$m))
    dBis = c(phi,val)
    uvinit = c(u0,v0)

    resultBis = genericLP(obj=objBis,A=ABis,modelsense="max",rhs=dBis,sense=c(rep(">",nbX*nbY),"="),start=uvinit)
    #
    u = resultBis$solution[1:nbX] 
    v = resultBis$solution[(nbX+1):(nbX+nbY)]
    
    residuals = Psi(market$transfers,matrix(u,nrow=nbX,ncol=nbY),matrix(v,nrow=nbX,ncol=nbY,byrow=T))
    #
    outcome = list(success=TRUE,
                   mu=mu,
                   mux0 = market$n - apply(mu,1,sum),
                   mu0y = market$m - apply(mu,2,sum),
                   u=u, v=v,
                   val=val,
                   residuals=residuals)
    #
    return(outcome)
}

updatev <- function(market, v, xFirst)
{
    tr = market$transfers
    
    n=market$n
    m=market$m
    
    nbX = length(n)
    nbY = length(m)
    #
    themat=matrix(0,nbX,nbY)
    vupdated = rep(0,nbY)
    #
    for(y in 1:nbY){    
        for(x in 1:nbX){
            for (yp in 1:nbY){
                themat[x,yp] = Vcal(tr,ifelse(yp==y,0,Ucal(tr,v[yp],x,yp)),x,y)
            }
        }
        #
        d = rep(0,nbX)
        A =  cbind2(sparseMatrix(1:nbX,1:nbX),rep(1,nbX))
        lb = c(-apply(themat,1,min),0)
        obj = c(market$n,market$m[y])
        #
        result = genericLP(obj=obj,A=A,modelsense="min",rhs=d,sense=">",lb=lb)
        #
        u0 = result$solution[1:nbX] 
        v0y = result$solution[nbX+1]
        val = result$objval
        #
        objBis = c(rep(0,nbX),1)
        ABis = rbind2(A,c(market$n,market$m[y]))
        modelsenseBis = ifelse(xFirst,"min","max")
        dBis = c(d,val)
        senseBis = c(rep(">",nbX),"=")
        uvinit = c(u0,v0y)
        
        resultBis = genericLP(obj=objBis,A=ABis,modelsense=modelsenseBis,rhs=dBis,sense=senseBis,lb=lb,start=uvinit)
        #
        u = resultBis$solution[1:nbX] 
        vupdated[y] = resultBis$solution[nbX+1]
    }
    #
    return(vupdated)  
}

eapNash <- function(market, xFirst=TRUE, notifications=FALSE, tol=1e-8, debugmode=FALSE)
{
    #
    if(!is.null(market$neededNorm)){
        stop("eapNash does not yet allow for the case without unmatched agents.")
    }
    nb_digits = 3
    #
    tr = market$transfers
    
    n = market$n
    m = market$m
    
    nbX = length(n)
    nbY = length(m)
    
    OutputFlag = 0  #ifelse(notifications,1,0)
    #
    if(notifications){
        message('Solving for equilibrium in ITU_none problem using Nash-ITU.')
    }
    if(xFirst){
        vcur = vfromus(tr,rep(0,nbX))$v
    }else{
        vcur = rep(0,nbY)
    }
    if(notifications & debugmode){
        message(c("vcur = ",round(vcur,nb_digits)))
    }
    #
    iter = 0
    cont = TRUE
    
    while(cont==TRUE){
        vnext = updatev(market,vcur,xFirst)
        if(max(abs(vcur-vnext)) < tol){
            cont = FALSE
        }
        vcur = vnext
        if(notifications & debugmode){
            message(c("vcur = ",round(vcur,nb_digits)))
        }
        #
        iter = iter + 1
    }
    #
    v = vcur
    #
    if(notifications){
        message(paste0("Nash-ITU converged in ",iter," iterations."))
    }
    #
    res = ufromvs(tr,v,tol)
    u = res$u
    
    A1 = Matrix::kronecker(matrix(1,1,nbY),sparseMatrix(1:nbX,1:nbX))
    A2 = Matrix::kronecker(sparseMatrix(1:nbY,1:nbY),matrix(1,1,nbX))
    A = rbind2(A1,A2)
    #
    sense = ifelse(abs(c(u,v) - 0) < tol, "<", "=")
    
    result = genericLP(obj=c(res$subdiff),A=A,modelsense="max",rhs=c(n,m),sense=sense)
    #
    mu = matrix(result$solution,nrow=nbX)
    #
    mux0 = n - apply(mu,1,sum)
    mu0y = m - apply(mu,2,sum)
    #
    outcome = list(mu=mu,
                   mux0=mux0, mu0y=mu0y,
                   u=u, v=v)
    #
    if(min( c(res$subdiff,(abs(u)<tol),(abs(v)<tol)) >= c((mu>0),(mux0>0),(mu0y>0)) )==0){
        if(notifications){
            message("Equilibrium not found.")
            message(paste0("mu = ",mu))
            message(paste0("Psi_xy(u_x,v_y) = ", round((residuals),2)))
            message(paste0("mux0 = ",mux0))
            message(paste0("u(x) = ",u))
            message(paste0("mu0y = ",mu0y))
            message(paste0("v(y) = ",v,"\n"))
        }
        outcome$success = FALSE
    }else{
        if(notifications){
            message("Equilibrium found.\n")
        }
    }
    #
    return(outcome)
}

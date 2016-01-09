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

marketTranspose <- function(market)
{
    thelist=list(n=market$m, m=market$n, 
                 hetG=market$hetH, hetH=market$hetG,
                 transfers=transfersTranspose(market$transfers),
                 neededNorm=normalizationTranspose(market$neededNorm))
    #
    add = 0
    names = names(market)
    # here, transpose additional elements; add 1 to add each time
    if(length(names) > 5 + add){
        message("Warning: in bipartite market transposition, 
         some elements have not been copied.")
    }
    #
    return(structure(thelist,class=class(market)))
}

normalizationTranspose <- function(neededNorm)
{
    if(!is.null(neededNorm)){
        names = names(neededNorm)
        thelist = list()
        nbElts = 0
        #
        if("H_edge_logit" %in% neededNorm){
            nbElts = nbElts + 1
            thelist$H_edge_logit = function(a,b) (neededNorm$H_edge_logit(1/b,1/a))
        }
        #
        return(thelist)  
    }else{
        return(NULL)
    }
}

outcomeTranspose <- function(outcome)
{
    names = names(outcome)
    thelist = list()
    nbElts = 0
    #
    if("mu" %in% names){
        nbElts = nbElts + 1
        thelist$mu = t(outcome$mu)
    }
    if("mux0" %in% names){
        nbElts = nbElts + 1
        thelist$mu0y = outcome$mux0
    }
    if("mu0y" %in% names){
        nbElts = nbElts + 1
        thelist$mux0 = outcome$mu0y
    }
    if("U" %in% names){
        nbElts = nbElts + 1
        thelist$V = t(outcome$U)
    }
    if("V" %in% names){
        nbElts = nbElts + 1
        thelist$U = t(outcome$V)
    }
    if("val" %in% names){
        nbElts = nbElts + 1
        thelist$val = outcome$val
    }
    if("u" %in% names){
        nbElts = nbElts + 1
        thelist$v = outcome$u
    }
    if("v" %in% names){
        nbElts = nbElts + 1
        thelist$u = outcome$v
    }
    if("success" %in% names){
        nbElts = nbElts + 1
        thelist$success = outcome$success
    }
    if("residuals" %in% names){
        nbElts = nbElts + 1
        thelist$residuals = outcome$residuals
    }  
    if(length(outcome) > nbElts){
        message("Warning: in outcome transposition, some elements have not been copied.")
    }
    #
    return(structure(thelist))
}

checkNorm <- function(neededNorm, n, m, hetG, hetH)
{
    if(sum(n)!=sum(m)){
        stop("Normalization asked but sum(n) does not coincide with sum(m)")
    }
    if(hetG$outsideOption==TRUE){
        stop("Normalization asked but hetG should not allow an outside option")
    }
    if(hetH$outsideOption==TRUE){
        stop("Normalization asked but hetH should not allow an outside option")
    }
}

#################################################
########     Methods for  markets       #########
#################################################
#################################################

margxInv <- function(xs, mkt, ...) UseMethod("margxInv",mkt)

margxInv.default <- function(xs, mkt, Mu0ys, sigma=1) 
{
    message("hello default") # Keith: should this be here?
    coeff = ifelse(is.null(mkt$neededNorm),1,0)
    #
    if(is.null(xs)){
        xs = 1:mkt$transfers$nbX
    }
    themux0s = rep(0,length(xs))
    #
    for(x in xs){
        root_fn <- function(z) (coeff*z - mkt$n[x] + sum(MMF(mkt$transfers,z,Mu0ys,xs=x,sigma=sigma)))
        themux0s[x] = uniroot(root_fn, c(0,mkt$n[x]), tol = 1e-300)$root # Keith: fix tolerence   
    }
    #
    return(themux0s)
}

margyInv <- function(ys, mkt, ...) UseMethod("margyInv",mkt)

margyInv.default <- function(ys, mkt, Mux0s, sigma=1)
{
    coeff = ifelse(is.null(mkt$neededNorm),1,0)
    #
    if(is.null(ys)){
        ys = 1:mkt$transfers$nbY
    }
    themu0ys = rep(0,mkt$transfers$nbY)
    #
    for(y in ys){
        root_fn <- function(z) (coeff*z - mkt$m[y] + sum(MMF(mkt$transfers,Mux0s,z,ys=y,sigma=sigma)))
        themu0ys[y] = uniroot(root_fn, c(0,mkt$m[y]), tol=1e-300)$root
    }
    #
    return(themu0ys)
}

solveEquilibrium <- function(market, ...) UseMethod("solveEquilibrium")

build_market_TU_none <- function(n, m, phi, neededNorm=NULL)
{
    if(!is.null(neededNorm) && (sum(n) != sum(m))){
        stop("Normalization asked but sum(n) does not coincide with sum(m)")
    }
    #
    nbX = length(n)
    nbY = length(m)
    #
    TUs = build_TUs(phi)
    noneM = build_none(nbX,nbY)
    noneW = build_none(nbY,nbX)
    #
    ret = list(n=n, m=m,
               hetG=noneM, hetH=noneW,
               transfers=TUs,
               neededNorm=neededNorm)
    class(ret) = "TU_none"
    #
    return(ret)  
}

solveEquilibrium.TU_none = oapLP

build_market_TU_general <- function(n, m, phi, hetG, hetH, neededNorm=NULL)
{
    if(!is.null(neededNorm)){
        checkNorm(neededNorm,n,m,hetG,hetH)
    }
    #
    nbX = length(n)
    nbY = length(m)
    #
    TUs = build_TUs(phi)
    #
    ret = list(n=n, m=m,
               hetG=hetG, hetH=hetH,
               transfers=TUs,
               neededNorm=neededNorm)
    class(ret) = "TU_general"
    #
    return(ret)
}

solveEquilibrium.TU_general = maxWelfare

build_market_TU_empirical <- function(n, m, phi, hetG, hetH, nbDraws, seed=NULL, neededNorm=NULL)
{
    if(!is.null(neededNorm)){
        checkNorm(neededNorm,n,m,hetG,hetH)
    }
    #
    hetGsim = simul(hetG,nbDraws,seed)
    hetHsim = simul(hetH,nbDraws,seed)
    #
    nbX = length(n)
    nbY = length(m)
    #
    TUs = build_TUs(phi)
    #
    ret = list(n=n, m=m,
               hetG=hetGsim, hetH=hetHsim,
               transfers=TUs,
               neededNorm=neededNorm)
    class(ret) = "TU_empirical"
    #
    return(ret)
}

solveEquilibrium.TU_empirical = CupidsLP

build_market_NTU_none <- function(n, m, alpha, gamma, neededNorm=NULL)
{
    if(!is.null(neededNorm) && (sum(n) != sum(m))){
        stop("Normalization asked but sum(n) does not coincide with sum(m)")
    }
    #
    nbX = length(n)
    nbY = length(m)
    #
    NTUs = build_NTUs(alpha,gamma)
    noneM = build_none(nbX,nbY)
    noneW = build_none(nbY,nbX)
    #
    ret = list(n=n, m=m,
               hetG=noneM, hetH=noneW,
               transfers=NTUs,
               neededNorm=neededNorm)
    class(ret) = "NTU_none"
    #
    return(ret)
}

solveEquilibrium.NTU_none = darum

build_market_NTU_general <- function(n, m, alpha, gamma, hetG, hetH, neededNorm=NULL)
{
    if(!is.null(neededNorm)){
        checkNorm(neededNorm,n,m,hetG,hetH)
    }
    #
    NTUs = build_NTUs(alpha,gamma)
    #
    ret = list(alpha=alpha, gamma=gamma,
               n=n,m=m,
               hetG=hetG, hetH=hetH,
               transfers=NTUs,
               neededNorm=neededNorm)
    class(ret) = "NTU_general"
    #
    return(ret)
}

solveEquilibrium.NTU_general = darum

build_market_LTU_none <- function(n, m, lambda, phi, neededNorm=NULL)
{
    if(!is.null(neededNorm) && (sum(n) != sum(m))){
        stop("Normalization asked but sum(n) does not coincide with sum(m)")
    }
    #
    nbX = length(n)
    nbY = length(m)
    #
    LTUs = build_LTUs(lambda,phi)
    noneM = build_none(nbX,nbY)
    noneW = build_none(nbY,nbX)
    #
    ret = list(n=n, m=m,
               hetG=noneM, hetH=noneW,
               transfers=LTUs,
               neededNorm=neededNorm)
    class(ret) = "LTU_none"
    #
    return(ret)
}

build_market_LTU_logit <- function(n, m, lambda, phi, sigma=1, neededNorm=NULL)
{
    if(!is.null(neededNorm) && (sum(n) != sum(m))){
        stop("Normalization asked but sum(n) does not coincide with sum(m)")
    }
    #
    nbX = length(n)
    nbY = length(m)
    #
    LTUs = build_LTUs(lambda,phi)
    logitM = build_logits(nbX,nbY,sigma,outsideOption=is.null(neededNorm))
    logitW = build_logits(nbY,nbX,sigma,outsideOption=is.null(neededNorm))
    #
    ret = list(n=n, m=m,
               hetG=logitM, hetH=logitW,
               transfers=LTUs,
               neededNorm=neededNorm)
    class(ret) = "LTU_logit"
    #
    return(ret)
}

build_market_ITU_logit <- function(n, m, transfers, sigma=1, neededNorm=NULL)
{
    if(!is.null(neededNorm) && (sum(n) != sum(m))){
        stop("Normalization asked but sum(n) does not coincide with sum(m)")
    }
    #
    nbX = length(n)
    nbY = length(m)
    #
    logitM = build_logits(nbX,nbY,sigma,outsideOption=is.null(neededNorm))
    logitW = build_logits(nbY,nbX,sigma,outsideOption=is.null(neededNorm))
    #
    ret = list(n=n, m=m,
               hetG=logitM, hetH=logitW,
               transfers=transfers) # Keith: should this include neededNorm?
    class(ret) = "ITU_logit"
    #
    return(ret)
}

solveEquilibrium.ITU_logit = ipfp

build_market_ITU_general <- function(n, m, hetG, hetH, transfers, neededNorm=NULL)
{
    ret = list(n=n, m=m,
               hetG=hetG, hetH=hetH,
               transfers=transfers,
               neededNorm=neededNorm)
    class(ret) = "ITU_general"
    #
    return(ret)
}

solveEquilibrium.ITU_general = jacobi

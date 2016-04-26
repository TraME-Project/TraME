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
  thelist=list(types = market$types,
               n=market$m, m=market$n, 
               neededNorm=normalizationTranspose(market$neededNorm),
               mmfs = mmfsTranspose(market$mmfs),
               arumsG=market$arumsH, arumsH=market$arumsG,
               transfers=transfersTranspose(market$transfers))
  #
  add = 0
  names = names(market)
  # here, transpose additional elements; add 1 to add each time
  if(length(names) > length(thelist) + add){
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

defaultNorm <- function(noSingles=FALSE)
{
  if(noSingles){
    return(list(H_edge_logit = function(mux0,mu0y) (mu0y[1])))
  }else{
    return(NULL)
  }
} 

checkNorm <- function(neededNorm, n, m, arumsG, arumsH)
{
  if(sum(n)!=sum(m)){
    stop("Normalization asked but sum(n) does not coincide with sum(m)")
  }
  if(arumsG$outsideOption==TRUE){
    stop("Normalization asked but arumsG should not allow an outside option")
  }
  if(arumsH$outsideOption==TRUE){
    stop("Normalization asked but arumsH should not allow an outside option")
  }
}

#################################################
########     Methods for  markets       #########
#################################################

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
  ret = list(types = "itu-rum",
             n=n, m=m,
             arumsG=noneM, arumsH=noneW,
             transfers=TUs,
             neededNorm=neededNorm)
  class(ret) = "TU_none"
  #
  return(ret)  
}

solveEquilibrium.TU_none = oapLP

build_market_TU_general <- function(n, m, phi, arumsG, arumsH, neededNorm=NULL)
{
  if(!is.null(neededNorm)){
    checkNorm(neededNorm,n,m,arumsG,arumsH)
  }
  #
  nbX = length(n)
  nbY = length(m)
  #
  TUs = build_TUs(phi)
  #
  ret = list(types = "itu-rum",
             n=n, m=m,
             arumsG=arumsG, arumsH=arumsH,
             transfers=TUs,
             neededNorm=neededNorm)
  class(ret) = "TU_general"
  #
  return(ret)
}

solveEquilibrium.TU_general = maxWelfare

build_market_TU_empirical <- function(n, m, phi, arumsG, arumsH, nbDraws, seed=NULL, neededNorm=NULL)
{
  if(!is.null(neededNorm)){
    checkNorm(neededNorm,n,m,arumsG,arumsH)
  }
  #
  arumsGsim = simul(arumsG,nbDraws,seed)
  arumsHsim = simul(arumsH,nbDraws,seed)
  #
  nbX = length(n)
  nbY = length(m)
  #
  TUs = build_TUs(phi)
  #
  ret = list(types = "itu-rum",
             n=n, m=m,
             arumsG=arumsGsim, arumsH=arumsHsim,
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
  ret = list(types = "itu-rum",
             n=n, m=m,
             arumsG=noneM, arumsH=noneW,
             transfers=NTUs,
             neededNorm=neededNorm)
  class(ret) = "NTU_none"
  #
  return(ret)
}

solveEquilibrium.NTU_none = darum

build_market_NTU_general <- function(n, m, alpha, gamma, arumsG, arumsH, neededNorm=NULL)
{
  if(!is.null(neededNorm)){
    checkNorm(neededNorm,n,m,arumsG,arumsH)
  }
  #
  NTUs = build_NTUs(alpha,gamma)
  #
  ret = list(types = "itu-rum",
             alpha=alpha, gamma=gamma,
             n=n,m=m,
             arumsG=arumsG, arumsH=arumsH,
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
  ret = list(types = "itu-rum",
             n=n, m=m,
             arumsG=noneM, arumsH=noneW,
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
  if(is.null(neededNorm)){
    outsideOption = TRUE
  }else{
    outsideOption = FALSE
  }
  #
  nbX = length(n)
  nbY = length(m)
  #
  LTUs = build_LTUs(lambda,phi)
  logitM = build_logits(nbX,nbY,sigma=sigma,outsideOption=outsideOption)
  logitW = build_logits(nbY,nbX,sigma=sigma,outsideOption=outsideOption)
  # 
  LTUmmfs = build_LTUmmfs(n,m,lambda,exp(phi/sigma),neededNorm)
  #
  ret = list(types = c("itu-rum","mfe"),
             n=n, m=m,
             neededNorm=neededNorm,
             #
             arumsG=logitM, arumsH=logitW,
             transfers=LTUs,
             #
             mmfs = LTUmmfs
  )
  class(ret) = "LTU_logit"
  #
  return(ret)
}

build_market_ETU_logit <- function(n, m, alpha,gamma,tau,sigma=1, neededNorm=NULL)
{
  if(!is.null(neededNorm) && (sum(n) != sum(m))){
    stop("Normalization asked but sum(n) does not coincide with sum(m)")
  }
  #
  nbX = length(n)
  nbY = length(m)
  #
  #
  ETUs = build_ETUs(alpha, gamma, tau)
  logitM = build_logits(nbX,nbY,sigma,outsideOption=is.null(neededNorm))
  logitW = build_logits(nbY,nbX,sigma,outsideOption=is.null(neededNorm))
  #
  ETUmmfs = build_ETUmmfs(n,m,exp(-alpha/tau),exp(-gamma/tau),-1/tau,neededNorm)
  #
  ret = list(types = c("itu-rum","mfe"),
             n=n, m=m,
             neededNorm=neededNorm,
             #
             arumsG=logitM, arumsH=logitW,
             transfers=ETUs,
             #
             mmfs = ETUmmfs
  ) 
  class(ret) = "ETU_logit"
  #
  return(ret)
}
#
solveEquilibrium.ETU_logit = ipfp
#
build_market_ITU_general <- function(n, m, arumsG, arumsH, transfers, neededNorm=NULL)
{
  ret = list(types = "itu-rum",
             n=n, m=m,
             arumsG=arumsG, arumsH=arumsH,
             transfers=transfers,
             neededNorm=neededNorm)
  class(ret) = "ITU_general"
  #
  return(ret)
}
#
solveEquilibrium.ITU_general = jacobi
#
build_market_NTU_logit <- function(n, m, alpha, gamma, sigma=1, neededNorm=NULL)
{
  if(!is.null(neededNorm) && (sum(n) != sum(m))){
    stop("Normalization asked but sum(n) does not coincide with sum(m)")
  }
  #
  nbX = length(n)
  nbY = length(m)
  #
  NTUs = build_NTUs(alpha,gamma)
  logitM = build_logits(nbX,nbY,sigma=sigma,outsideOption=is.null(neededNorm))
  logitW = build_logits(nbY,nbX,sigma=sigma,outsideOption=is.null(neededNorm))
  #
  NTUmmfs = build_NTUmmfs(n,m,exp(alpha/sigma),exp(gamma/sigma),neededNorm)
  #
  ret = list(types = c("itu-rum","mfe"),
             n=n,m=m,
             neededNorm=neededNorm,
             #
             arumsG=logitM,arumsH=logitW,
             transfers=NTUs,
             #
             mmfs = NTUmmfs
  )
  class(ret) = "NTU_logit"
  #
  return(ret)
}

solveEquilibrium.NTU_logit = ipfp


build_market_TU_logit <- function(n, m, phi, sigma=1, neededNorm=NULL)
{
  if(!is.null(neededNorm) && (sum(n) != sum(m))){
    stop("Normalization asked but sum(n) does not coincide with sum(m)")
  }
  #
  nbX = length(n)
  nbY = length(m)
  #
  TUs = build_TUs(phi)
  logitM = build_logits(nbX,nbY,sigma=sigma,outsideOption=is.null(neededNorm))
  logitW = build_logits(nbY,nbX,sigma=sigma,outsideOption=is.null(neededNorm))
  #
  TUmmfs = build_TUmmfs(n,m,exp(phi/(2*sigma)),neededNorm)
  #
  ret = list(types = c("itu-rum","mfe"),
             n=n,m=m,
             neededNorm=neededNorm,
             #
             arumsG=logitM, arumsH=logitW,
             transfers=TUs,
             #
             mmfs = TUmmfs
  )
  class(ret) = "TU_logit"
  #
  return(ret)
}

solveEquilibrium.TU_logit = ipfp


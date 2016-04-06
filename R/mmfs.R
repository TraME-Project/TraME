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
# References:
#
# A. Galichon, S. Weber: "Estimation of Matching Function Equilibria"
#           
#           
# mmfs has n, and m and neededNorm
# M.mmfs <- function(mmfs, mux0, mu0y)



margxInv.default <- function(xs, mmfs, Mu0ys) 
{
  nbX = length(mmfs$n)
  if (is.null(mmfs$neededNorm))
  {
    coeff = 1
    ubs = mmfs$n
  }
  else
  {
    coeff = 0
    ubs = rep(1e10,nbX)
  }
  #
  if(is.null(xs)){
    xs = 1:nbX
  }
  themux0s = rep(0,length(xs))
  #
  for(x in xs){
    root_fn <- function(z) (coeff*z - mmfs$n[x] + sum(M(mmfs,z,Mu0ys,xs=x)))
    themux0s[x] = uniroot(root_fn, c(0,ubs[x]), tol = 1e-300)$root # Keith: fix tolerence   
  }
  #
  return(themux0s)
}

margyInv.default <- function(ys, mmfs, Mux0s)
{
  nbY = length(mmfs$m)
  if (is.null(mmfs$neededNorm))
  {
    coeff = 1
    ubs = mmfs$m
  }
  else
  {
    coeff = 0
    ubs = rep(1e10,nbY)
  }
  #
  if(is.null(ys)){
    ys = 1:nbY
  }
  themu0ys = rep(0,nbY)
  #
  for(y in ys){
    root_fn <- function(z) (coeff*z - mmfs$m[y] + sum(M(mmfs,Mux0s,z,ys=y)))
    themu0ys[y] = uniroot(root_fn, c(0,ubs[y]), tol=1e-300)$root
  }
  #
  return(themu0ys)
}
#########################  TU #########################################
build_TUmmfs <- function(n,m,phi,neededNorm)
{
  ret = list(n=n,
             m=m,
             neededNorm=neededNorm,
             expphiover2 = exp(phi/2)
  )
  class(ret)="TUmmfs"
  return(ret)
}
#
M.TUmmfs <- function(mmfs, mux0s, mu0ys, xs=length(mmfs$n), ys=length(mmfs$m))
{
  term_1 = mmfs$expphiover2[xs,ys]
  term_2 = sqrt(mux0s %*% t(mu0ys))
  #
  ret = term_1 * term_2
  #
  return(ret)
}
#
margxInv.TUmmfs <- function(xs, mmfs, Mu0ys)
{
  if (is.null(xs)) {xs = 1:length(mmfs$n)}
  #
  sqrtMu0ys = sqrt(Mu0ys)
  if(is.null(mmfs$neededNorm)){
    b = (mmfs$expphiover2[xs,] %*% sqrtMu0ys)/2
    sqrtMux0s = sqrt(mmfs$n[xs]+ b*b) - b
  }else{
    sqrtMux0s = mmfs$n / c(mmfs$expphiover2[xs,] %*% sqrtMu0ys)
  }
  #
  ret = c(sqrtMux0s*sqrtMux0s)
  #
  return(ret)
}
#
margyInv.TUmmfs <- function(ys, mmfs, Mux0s)
{
  if (is.null(ys)) {xs = 1:length(mmfs$m)}
  #
  sqrtMux0s = sqrt(Mux0s)
  if(is.null(mmfs$neededNorm)){
    b = c(sqrtMux0s %*% mmfs$expphiover2[,ys])/2
    sqrtMu0ys = sqrt(mmfs$m[ys] + b*b) - b
  }else{
    sqrtMu0ys = mkt$m / c(sqrtMux0s %*% mmfs$expphiover2[,ys])
  }
  #
  ret = sqrtMu0ys*sqrtMu0ys
  #
  return(ret)
}
######################### LTU #########################################
build_LTUmmfs <- function(n,m,lambda,phi,neededNorm)
{
  ret = list(n=n,
             m=m,
             neededNorm=neededNorm,
             lambda = lambda,
             expphi=exp(phi),
             aux_zeta = 1-lambda)
  class(ret)="LTUmmfs"
  return(ret)
}
#
M.LTUmmfs <- function(mmfs, mux0s, mu0ys, xs=length(mmfs$n), ys=length(mmfs$m))
{
  term_1 = mux0s^mmfs$lambda[xs,ys]
  term_2 = t( mu0ys^t(mmfs$aux_zeta[xs,ys]) )
  term_3 = mmfs$expphi[xs,ys]
  #
  ret = term_1 * term_2 * term_3
  #
  return(ret)
}
#
######################### NTU #########################################
build_NTUmmfs <- function(n,m,alpha,gamma,neededNorm)
{
  ret = list(n=n,
             m=m,
             neededNorm=neededNorm,
             expalpha = exp(alpha),
             expgamma=exp(gamma))
  class(ret)="NTUmmfs"
  return(ret)
}
#
M.NTUmmfs <- function(mmfs, mux0s, mu0ys, xs=1:length(mmfs$n), ys=length(mmfs$m))
{
  term_1 = mux0s * mmfs$expalpha[xs,ys]
  term_2 = t( mu0ys * t(mmfs$expgamma[xs,ys] ))
  #
  ret = pmin(term_1, term_2)
  #
  return(ret)
}
#
margxInv.NTUmmfs = function(xs,mmfs,Mu0ys)
{
  if (is.null(xs)) {xs = 1:length(mmfs$n)}
  if (!is.null(mmfs$neededNorm)) {stop('not supported yet')}
  themux0s = inversePWA(mmfs$n[xs], t ( t(mmfs$expgamma[xs,] / mmfs$expalpha[xs,])* Mu0ys ),mmfs$expalpha[xs,]) 
  return(themux0s)
}
#
margyInv.NTUmmfs = function(ys,mmfs,Mux0s)
{
  if (is.null(ys)) {ys = 1:length(mmfs$m)}
  if (!is.null(mmfs$neededNorm)) {stop('not supported yet')}
  themu0ys = inversePWA(mmfs$m[ys], t(( mmfs$expalpha[,ys] / mmfs$expgamma[,ys]) * Mux0s), t(mmfs$expgamma[,ys] ) ) 
  return(themu0ys)
}
#
######################### ITU #########################################
#
build_ETUmmfs <- function(n,m,alpha,gamma,tau,neededNorm)
{
  ret = list(n=n,
             m=m,
             neededNorm=neededNorm,
             expminusalphaovertau = exp(-alpha/tau),
             expminusgammaovertau = exp(-gamma/tau),
             tauinv = 1/tau
  )
  class(ret)="ETUmmfs"
}
#
M.ETUmmfs <- function(mmfs, mux0s, mu0ys, xs=1:length(mmfs$n), ys=length(mmfs$m))
{
  
  term_1 = mmfs$expminusalphaovertau[xs,ys] / (mux0s^mmfs$tauinv)
  term_2 = mmfs$expminusgammaovertau[xs,ys] / t(mu0ys^t(mmfs$tauinv))
  term_exp = mmfs$tau[xs,ys]
  #
  ret = (2/(term_1 + term_2))^(1/mmfs$tauinv)
  #
  return(ret)
}


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
#########################  TU #########################################
build_LTUmmfs <- function(n,m,phi,neededNorm)
{
  ret = list(n=n,
             m=m,
             neededNorm=neededNorm,
             expphiover2 = exp(phi/2)
             )
  class(ret)="TUmmfs"
}
#
M.TUmmfs <- function(mmfs, mux0s, mu0ys, xs=1:tr$nbX, ys=1:tr$nbY)
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
           aux_zeta = aux_zeta)
class(ret)="LTUmmfs"
}
#
M.LTUmmfs <- function(mmfs, mux0s, mu0ys, xs=1:tr$nbX, ys=1:tr$nbY)
{
  term_1 = mux0s^mmfs$lambda[xs,ys]
  term_2 = t( mu0ys^t(mmfs$aux_zeta[xs,ys]) )
  term_3 = tr$expphi[xs,ys]
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
}
#
M.NTUmmfs <- function(mmfs, mux0s, mu0ys, xs=1:tr$nbX, ys=1:tr$nbY)
{
  term_1 = mux0s * mmfs$expalpha[xs,ys]
  term_2 = t( mu0ys * t(mmfs$expgamma[xs,ys] )
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
  themu0ys = inversePWA(mmfs$m[ys], t(( mmfs$expalpha[,ys] / mmfs$expgamma[,ys]) * Mux0s), t(mmfs$expgamma[,ys] ) 
  return(themu0ys)
}
#
######################### ITU #########################################
#

build_ITUmmfs <- function(n,m,transfers,neededNorm)
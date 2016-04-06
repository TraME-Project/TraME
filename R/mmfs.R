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


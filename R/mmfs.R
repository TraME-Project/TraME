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
# M.mmfs <- function(mmfs, ax, by)
#
################################################################################
########################    Default and generic methods ########################
################################################################################
#
margxInv.default <- function(xs, mmfs, Bys)
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
  theaxs = rep(0,length(xs))
  #
  i = 0
  for(x in xs){
    i = i+1
    root_fn <- function(z) (ifelse(coeff,Mx0(mmfs,z),0) - mmfs$n[x] + sum(M(mmfs,z,Bys,xs=x)))
    theaxs[i] = uniroot(root_fn, c(0,ubs[x]), tol = 1e-300)$root # Keith: fix tolerence   
  }
  #
  return(theaxs)
}
#
margyInv.default <- function(ys, mmfs, Axs)
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
  thebys = rep(0,length(ys))
  #
  j = 0
  for(y in ys){
    j = j+1
    root_fn <- function(z) (ifelse(coeff,M0y(mmfs,z),0) - mmfs$m[y] + sum(M(mmfs,Axs,z,ys=y)))
    thebys[j] = uniroot(root_fn, c(0,ubs[y]), tol=1e-300)$root
  }
  #
  return(thebys)
}
#
Mx0.default <- function(mmfs, ax, ...)
{
  return(ax)
}
#
M0y.default <- function(mmfs, by, ...)
{
  return(by)
}
#
################################################################################
########################             TU MMfs            ########################
################################################################################
build_TUmmfs <- function(n,m,K,neededNorm)
{
  ret = list(n=n,
             m=m,
             neededNorm=neededNorm,
             K = K
  )
  class(ret)="TUmmfs"
  return(ret)
}
#
mmfsTranspose.TUmmfs <- function(mmfs)
{
  ret = list(n=mmfs$m,
             m=mmfs$n,
             neededNorm=normalizationTranspose(mmfs$neededNorm),
             K = t(mmfs$K)
  )
  class(ret)="TUmmfs"
  return(ret)
}
#
M.TUmmfs <- function(mmfs, axs, bys, xs=1:length(mmfs$n), ys=1:length(mmfs$m))
{
  term_1 = mmfs$K[xs,ys]
  term_2 = sqrt(axs %*% t(bys))
  #
  ret = term_1 * term_2
  #
  return(ret)
}
#
margxInv.TUmmfs <- function(xs, mmfs, Bys)
{
  if (is.null(xs)) {xs = 1:length(mmfs$n)}
  #
  sqrtBys = sqrt(Bys)
  if(is.null(mmfs$neededNorm)){
    b = (mmfs$K[xs,] %*% sqrtBys)/2
    sqrtAxs = sqrt(mmfs$n[xs]+ b*b) - b
  }else{
    sqrtAxs = mmfs$n / c(mmfs$K[xs,] %*% sqrtBys) # Keith: should this be n[xs] ?
  }
  #
  ret = c(sqrtAxs*sqrtAxs)
  #
  return(ret)
}
#
margyInv.TUmmfs <- function(ys, mmfs, Axs)
{
  if (is.null(ys)) {xs = 1:length(mmfs$m)}
  #
  sqrtAxs = sqrt(Axs)
  if(is.null(mmfs$neededNorm)){
    b = c(sqrtAxs %*% mmfs$K[,ys])/2 # Keith: will these be conformable? maybe need t(sqrtAxs)
    sqrtBys = sqrt(mmfs$m[ys] + b*b) - b
  }else{
    sqrtBys = mmfs$m / c(sqrtAxs %*% mmfs$K[,ys])
  }
  #
  ret = sqrtBys*sqrtBys
  #
  return(ret)
}

################################################################################
########################            NTU MMfs            ########################
################################################################################
build_NTUmmfs <- function(n,m,A,B,neededNorm)
{
  ret = list(n=n,
             m=m,
             neededNorm=neededNorm,
             A = A,
             B = B)
  class(ret)="NTUmmfs"
  return(ret)
}
#
mmfsTranspose.NTUmmfs <- function(mmfs)
{
  ret = list(n = mmfs$m,
             m = mmfs$n,
             neededNorm = normalizationTranspose(mmfs$neededNorm),
             A = t(mmfs$B),
             B = t(mmfs$A)
  )
  class(ret)="NTUmmfs"
  return(ret)
}
#
M.NTUmmfs <- function(mmfs, axs, bys, xs=1:length(mmfs$n), ys=1:length(mmfs$m))
{
  term_1 = axs * mmfs$A[xs,ys]
  term_2 = t( bys * t(mmfs$B[xs,ys] ))
  #
  ret = pmin(term_1, term_2)
  #
  return(ret)
}
#
margxInv.NTUmmfs <- function(xs,mmfs,Bys) # Keith: shouldn't this have mmfs first?
{
  if (is.null(xs)) {xs = 1:length(mmfs$n)}
  if (!is.null(mmfs$neededNorm)) {stop('not supported yet')}
  theaxs = inversePWA(mmfs$n[xs], t ( t(mmfs$B[xs,] / mmfs$A[xs,])* Bys ),mmfs$A[xs,]) 
  return(theaxs)
}
#
margyInv.NTUmmfs <- function(ys,mmfs,Axs) # Keith: shouldn't this have mmfs first?
{
  if (is.null(ys)) {ys = 1:length(mmfs$m)}
  if (!is.null(mmfs$neededNorm)) {stop('not supported yet')}
  thebys = inversePWA(mmfs$m[ys], t(( mmfs$A[,ys] / mmfs$B[,ys]) * Axs), t(mmfs$B[,ys] ) ) 
  return(thebys)
}
################################################################################
########################            LTU MMfs            ########################
################################################################################
build_LTUmmfs <- function(n,m,lambda,K,neededNorm)
{
  ret = list(n=n,
             m=m,
             neededNorm=neededNorm,
             lambda = lambda,
             K=K,
             aux_zeta = 1-lambda)
  class(ret)="LTUmmfs"
  return(ret)
}
#
mmfsTranspose.LTUmmfs <- function(mmfs)
{
  ret = list(n=mmfs$m,
             m=mmfs$n,
             neededNorm=normalizationTranspose(mmfs$neededNorm),
             lambda = t(mmfs$aux_zeta),
             K=t(mmfs$K),
             aux_zeta = t(mmfs$lambda)
             )
  class(ret)="LTUmmfs"
  return(ret)
}
#
M.LTUmmfs <- function(mmfs, axs, bys, xs=1:length(mmfs$n), ys=1:length(mmfs$m))
{
  term_1 = axs^mmfs$lambda[xs,ys]
  term_2 = t( bys^t(mmfs$aux_zeta[xs,ys]) )
  term_3 = mmfs$K[xs,ys]
  #
  ret = term_1 * term_2 * term_3
  #
  return(ret)
}
################################################################################
########################            ETU MMfs            ########################
################################################################################
build_ETUmmfs <- function(n,m,C,D,kappa,neededNorm)
{
  ret = list(n=n,
             m=m,
             neededNorm=neededNorm,
             C = C,
             D = D,
             kappa = kappa
  )
  class(ret)="ETUmmfs"
  return(ret)
}
#
mmfsTranspose.ETUmmfs <- function(mmfs)
{
  ret = list(n = mmfs$m,
             m = mmfs$n,
             neededNorm = normalizationTranspose(mmfs$neededNorm),
             C = t(mmfs$D),
             D =  t(mmfs$C),
             kappa = t(mmfs$kappa)
             )
  class(ret)="ETUmmfs"
  return(ret)
}
#
M.ETUmmfs <- function(mmfs, axs, bys, xs=1:length(mmfs$n), ys=1:length(mmfs$m))
{
  
  term_1 = mmfs$C[xs,ys] * (axs^mmfs$kappa[xs,ys])
  term_2 = mmfs$D[xs,ys] * t(bys^t(mmfs$kappa[xs,ys]))
  #
  ret = ((term_1 + term_2)/2)^(1/mmfs$kappa[xs,ys])
  #
  return(ret)
}


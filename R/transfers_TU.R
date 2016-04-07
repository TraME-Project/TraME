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
# A. Galichon, S.D. Kominers, and S. Weber: "An Empirical Framework for Matching with Imperfectly Transferable Utility"  
#           

build_TUs <- function(phi)
{
    nbX=dim(phi)[1]
    nbY=dim(phi)[2]
    #
    ret = list(nbX = nbX, nbY = nbY,
               nbParams = nbX*nbY,
               phi = phi)
    class(ret) = "TU"
    #
    return(ret)
}
#
transfersTranspose.TU <- function(tr)
{
  ret = list(nbX = tr$nbY, nbY = tr$nbX,
             nbParams = tr$nbParams,
             phi = t(tr$phi)
  )
  class(ret) = "TU"
  #
  return(ret)
}
#
Psi.TU <- function(tr, U, V)
{
    return((U + V - tr$phi)/2)
}

Psi_sub.TU <- function(tr,U,V,xs,ys) # here for vectorization purposes only 
{
  return((U + V - tr$phi[xs,ys])/2)
}
  

du_Psi.TU <- function(tr, U, V)
{
    ret = matrix(1/2, nrow=tr$nbX, ncol=tr$nbY)
    return(ret)
}

du_Psi_sub.TU <- function(tr, U, V, xs, ys)
{
  ret = matrix(1/2, nrow=length(xs), ncol=length(ys))
  return(ret)
}

dtheta_Psi.TU <- function(tr, U, V, dtheta=NULL) 
{
    ret <- 0
    if(is.null(dtheta)){
        ret = Diagonal(tr$nbX*tr$nbY,-1/2)
    }else{
        ret = -dtheta/2
    }
    return(ret)
}

determineType.TU <- function(tr, ...)
{
    return(1)
}


Ucal.TU <- function(tr, vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
    ret = tr$phi[xs,ys] - matrix(vs,length(xs),length(ys),byrow=TRUE)
    return(ret)
}

Vcal.TU <- function(tr, us, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(tr$phi[xs,ys] - us)
}


WU.TU <- function(tr, Us, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(2*Us - tr$phi[xs,ys])
}

WV.TU <- function(tr, Vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(tr$phi[xs,ys] - 2*Vs)
}

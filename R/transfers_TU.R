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
               phi = phi,
               aux_expphiover2 = exp(phi/2))
    class(ret) = "TU"
    #
    return(ret)
}

Psi.TU <- function(tr, U, V)
{
    return((U + V - tr$phi)/2)
}

du_Psi.TU <- function(tr, U, V)
{
    ret = matrix(1/2, nrow=tr$nbX, ncol=tr$nbY)
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

determineType.TU <- function(tr, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(1)
}

MMF.TU <- function(tr, mux0s, mu0ys, xs=1:tr$nbX, ys=1:tr$nbY, sigma=1)
{
    term_1 = tr$aux_expphiover2[xs,ys]^(1/sigma)
    term_2 = sqrt(mux0s %*% t(mu0ys))
    #
    ret = term_1 * term_2
    #
    return(ret)
}

Ucal.TU <- function(tr, vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
    ret = tr$phi[xs,ys] - matrix(vs,tr$nbX,tr$nbY,byrow=TRUE)
    return(ret)
}

Vcal.TU <- function(tr, us, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(tr$phi[xs,ys] - us)
}

UW.TU <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return((tr$phi[xs,ys] + Ws)/2)
}

dw_UW.TU <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    ret = matrix(1/2,length(xs),length(ys))
    return(ret)
}

VW.TU <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return((tr$phi[xs,ys]-Ws)/2)
}

dw_VW.TU <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    ret = matrix(-1/2,length(xs),length(ys))
    return(ret)
}

WU.TU <- function(tr, Us, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(2*Us - tr$phi[xs,ys])
}

WV.TU <- function(tr, Vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(tr$phi[xs,ys] - 2*Vs)
}

transfersTranspose.TU <- function(tr)
{
    ret = list(nbX = tr$nbY, nbY = tr$nbX,
               nbParams = tr$nbParams,
               phi = t(tr$phi),
               aux_expphiover2 = t(tr$aux_expphiover2))
    class(ret) = "TU"
    #
    return(ret)
}
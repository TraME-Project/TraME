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

build_NTUs <- function(alpha, gamma) 
{
    nbX = dim(alpha)[1]
    nbY = dim(gamma)[2]
    #
    ret = list(nbX = nbX, nbY = nbY,
               nbParams = 2 * nbX * nbY,
               alpha = alpha, gamma = gamma,
               aux_expalpha = exp(alpha),
               aux_expgamma = exp(gamma))
    class(ret) = "NTU"
    #
    return(ret)
}

Psi.NTU <- function(tr, U, V)
{
    return(pmax(U - tr$alpha, V - tr$gamma))
}

du_Psi.NTU <- function(tr, U, V)
{
    ret = ifelse(U-tr$alpha >= V - tr$gamma,1,0)
    return(ret)
}

dtheta_Psi.NTU <- function(tr, U, V, dtheta=NULL) 
{
    dupsi = c(du_Psi(tr,U,V))
    if(is.null(dtheta)){
        ret = -cbind(Diagonal(x=dupsi),
                     Diagonal(x=1- dupsi)) 
        return(ret)
    }else{
        dtheta1 = dtheta[1:(tr$nbX*tr$nbY)]
        dtheta2 = dtheta[(1+tr$nbX*tr$nbY):(2*tr$nbX*tr$nbY)]
        #
        ret = -c(dupsi*dtheta1 + (1-dupsi)*dtheta2)
        return(ret)
    }
}

determineType.NTU <- function(tr, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(2)
}

MMF.NTU <- function(tr, mux0s, mu0ys, xs=1:tr$nbX, ys=1:tr$nbY, sigma=1)
{
    term_1 = mux0s * tr$aux_expalpha[xs,ys]^(1/sigma)
    term_2 = t( mu0ys * t(tr$aux_expgamma[xs,ys]^(1/sigma)) )
    #
    ret = pmin(term_1, term_2)
    #
    return(ret)
}

UW.NTU <- function(tr, ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(pmin(tr$alpha[xs,ys], ws+tr$gamma[xs,ys]))
}

dw_UW.NTU <- function(tr, ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    ret = ifelse(tr$alpha[xs,ys]>ws+tr$gamma[xs,ys],1,0)
    return(ret)
}

VW.NTU <- function(tr, ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(pmin(tr$alpha[xs,ys]-ws,tr$gamma[xs,ys]))
}

dw_VW.NTU <- function(tr, ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    ret = ifelse(tr$alpha[xs,ys]-ws<tr$gamma[xs,ys],-1,0)
    return(ret)
}

WU.NTU <- function(tr, Us, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(Us - tr$alpha[xs,ys])
}

WV.NTU <- function(tr, Vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
    return(tr$gamma[xs,ys] - Vs)
}

transfersTranspose.NTU <- function(tr)
{
    ret = list(nbX = tr$nbY, nbY = tr$nbX,
               nbParams = tr$nbParams,
               alpha = t(tr$gamma), gamma = t(tr$alpha),
               aux_expalpha = t(tr$aux_expgamma),
               aux_expgamma = t(tr$aux_expalpha))
    class(ret) = "NTU"
    #
    return(ret)
}

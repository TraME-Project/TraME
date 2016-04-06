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

build_LTUs <- function(lambda, phi) 
{
    nbX = dim(lambda)[1]
    nbY = dim(lambda)[2]
    #
    aux_zeta = 1 - lambda
    if(min(c(lambda,aux_zeta)) <= 0){
        stop ("lambda not strictly between 0 and 1")
    }
    #
    ret = list(nbX = nbX, nbY = nbY,
               nbParams = 2* nbX * nbY,
               lambda = lambda, phi = phi, 
               aux_zeta = aux_zeta)
    class(ret) = "LTU"
    #
    return(ret)
}

Psi.LTU <- function(tr, U, V)
{
    ret = tr$lambda * U + tr$aux_zeta * V - tr$phi
    return(ret)
}

du_Psi.LTU <- function(tr, U, V) ( tr$lambda )

dtheta_Psi.LTU <- function(tr, U, V, dtheta=NULL) 
{
    UminusV = c(U - V)
    dtheta1 <- dtheta2 <- ret <- 0
    #
    if(is.null(dtheta)){
        ret = cbind(Diagonal(x=UminusV), Diagonal(tr$nbX*tr$nbY,-1))
        return(ret)
    }else{
        dtheta1 = dtheta[1:(tr$nbX*tr$nbY)]
        dtheta2 = dtheta[(1+tr$nbX*tr$nbY):(2*tr$nbX*tr$nbY)]
        #
        ret = c(UminusV*dtheta1 - dtheta2)
        return(ret)
    }
}

determineType.LTU <- function(tr, xs=1:tr$nbX, ys=1:tr$nbY) (1)


Ucal.LTU <- function(tr, vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
    term_1 = tr$phi[xs,ys]
    term_2 = tr$aux_zeta[xs,ys] * matrix(vs,length(xs),length(ys),byrow=TRUE)
    ret = (term_1 - term_2) / tr$lambda[xs,ys]
    #
    return(ret)
}

Vcal.LTU <- function(tr, us, xs=1:tr$nbX, ys=1:tr$nbY)
{
    term_1 = tr$phi[xs,ys]
    term_2 = tr$lambda[xs,ys] * matrix(us,length(xs),length(ys))
    ret = (term_1 - term_2) / tr$aux_zeta[xs,ys]
    #
    return(ret)
}

UW.LTU <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    ret = tr$phi[xs,ys] + tr$aux_zeta[xs,ys]*Ws
    return(ret)
}

dw_UW.LTU <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    ret = tr$aux_zeta[xs,ys]
    return(ret)
}

VW.LTU <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    ret = tr$phi[xs,ys] - tr$lambda[xs,ys]*Ws
    return(ret)
}

dw_VW.LTU <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    ret = -tr$lambda[xs,ys]
    return(ret)
}

WU.LTU <- function(tr, Us, xs=1:tr$nbX, ys=1:tr$nbY)
{
    numer = Us - tr$phi[xs,ys]
    denom = tr$aux_zeta[xs,ys]
    #
    ret = numer / denom
    #
    return(ret)
}

WV.LTU <- function(tr, Vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
    numer = tr$phi[xs,ys] - Vs
    denom = tr$lambda[xs,ys] 
    #
    ret = numer / denom
    #
    return(ret)
}

transfersTranspose.LTU <- function(tr)
{
    ret = list(nbX = tr$nbY, nbY = tr$nbX,
               nbParams = tr$nbParams,
               lambda = t(tr$aux_zeta), phi = t(tr$phi),
               aux_zeta = t(tr$lambda))
    class(ret) = "LTU"
    #
    return(ret)
}
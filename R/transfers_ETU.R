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
# A. Galichon, S.D. Kominers, and S. Weber: "An Empirical Framework for Matching with Imperfectly Transferable Utility". 
#

build_ETUs <- function(alpha, gamma, tau)
{
    nbX = dim(alpha)[1]
    nbY = dim(alpha)[2]
    #
    ret = list(nbX = nbX, nbY = nbY,
               nbParams = 3*nbX*nbY,
               alpha=alpha, gamma=gamma, tau=tau, 
               aux_expminusalphaovertau = exp(-alpha/tau),
               aux_expminusgammaovertau = exp(-gamma/tau))
    class(ret) = "ETU"
    #
    return(ret)
}

Psi.ETU <- function(tr, U, V)
{
    term_1 = exp(U/tr$tau)*tr$aux_expminusalphaovertau
    term_2 = exp(V/tr$tau)*tr$aux_expminusgammaovertau
    #
    ret = tr$tau * log((term_1 + term_2)/2)
    #
    return(ret)
}

du_Psi.ETU <- function(tr, U, V)
{
    term_1 = V - U + tr$alpha - tr$gamma
    term_2 = tr$tau
    #
    ret = 1/(1 + exp(term_1/term_2))
    #
    return(ret)
}

dtheta_Psi.ETU <- function(tr, U, V, dtheta=NULL) 
{
    dupsimat = du_Psi(tr,U,V)
    dupsi = c(dupsimat)
    #
    term_1 <- term_2 <- ret <- 0
    #
    if(is.null(dtheta)){
        term_1 = (U - tr$alpha )*dupsi
        term_2 = (V - tr$gamma)*(1-dupsi)
        #
        dsigmapsi = c((Psi(tr,U,V) - term_1 - term_2) / tr$tau)
        #
        ret = cbind(Diagonal(x = -dupsi),
                    Diagonal(x = dupsi-1),
                    Diagonal(x = dsigmapsi))
        #
        return(ret)
    }else{
        dtheta1 = dtheta[1:(tr$nbX*tr$nbY),]
        dtheta2 = dtheta[(1+tr$nbX*tr$nbY):(2*tr$nbX*tr$nbY),]
        dtheta3 = dtheta[(1+2*tr$nbX*tr$nbY):(3*tr$nbX*tr$nbY),]
        #
        if(min(dtheta3==0)){
            dsigmapsidtheta = 0
        }else{
            term_1 = (U - tr$alpha )*dupsimat
            term_2 = (V - tr$gamma)*(1-dupsimat)
            #
            dsigmapsidtheta = dtheta3*c((Psi(tr,U,V) - term_1 - term_2) / tr$tau)
        }
        #
        ret = c(-dupsi*dtheta1 - (1-dupsi)*dtheta2 + dsigmapsidtheta)
        #
        return(ret)
    }
}

determineType.ETU <- function(tr, xs=1:tr$nbX, ys=1:tr$nbY) (2)


Ucal.ETU <- function(tr, vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
    term_1 = matrix(vs,tr$nbX,tr$nbY,byrow=TRUE) - tr$gamma[xs,ys]
    term_log = 2 - exp(term_1/tr$tau[xs,ys])
    #
    ret = tr$alpha[xs,ys] + tr$tau[xs,ys] * log(term_log)
    #
    return(ret)
}

Vcal.ETU <- function(tr, us, xs=1:tr$nbX, ys=1:tr$nbY)
{
    term_1 = us - tr$alpha[xs,ys]
    term_log = 2 - exp(term_1/tr$tau[xs,ys])
    #
    ret = tr$gamma[xs,ys] + tr$tau[xs,ys] * log(term_log)
    #
    return(ret)
}

UW.ETU <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    term_1 = tr$aux_expminusalphaovertau[xs,ys]
    term_2 = exp(-Ws/tr$tau[xs,ys]) * tr$aux_expminusgammaovertau[xs,ys]
    term_log = (term_1 + term_2)/2
    #
    ret = -tr$tau[xs,ys] * log(term_log)
    #
    return(ret)
    
}

VW.ETU <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
    term_1 = tr$aux_expminusgammaovertau[xs,ys]
    term_2 = exp(Ws/tr$tau[xs,ys]) * tr$aux_expminusalphaovertau[xs,ys]
    term_log = (term_2 + term_1)/2
    #
    ret = -tr$tau[xs,ys] * log(term_log)
    #
    return(ret)
}

WU.ETU <- function(tr, Us, xs=1:tr$nbX, ys=1:tr$nbY)
{
    term_1 = 2*exp( (tr$gamma[xs,ys] - Us)/tr$tau[xs,ys] )
    term_2 = exp( (tr$gamma[xs,ys] - tr$alpha[xs,ys])/tr$tau[xs,ys] )
    term_log = term_1 - term_2
    #
    ret = -tr$tau[xs,ys] * log(term_log)
    #
    return(ret)
}

WV.ETU <- function(tr, Vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
    term_1 = 2*exp( (tr$alpha[xs,ys] - Vs)/tr$tau[xs,ys] )
    term_2 = exp( (tr$alpha[xs,ys] - tr$gamma[xs,ys])/tr$tau[xs,ys] )
    term_log = term_1 - term_2
    #
    ret = tr$tau[xs,ys] * log(term_log)
}

transfersTranspose.ETU <- function(tr)
{
    ret = list(nbX = tr$nbY, nbY = tr$nbX,
               nbParams = tr$nbParams,
               alpha=t(tr$gamma), gamma=t(tr$alpha), tau=t(tr$tau),
               aux_expminusalphaovertau = t(tr$aux_expminusgammaovertau),
               aux_expminusgammaovertau = t(tr$aux_expminusalphaovertau))
    class(ret) = "ETU"
    #
    return(ret)
}
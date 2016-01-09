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

build_market_TU_logit <- function(n, m, phi, sigma=1, neededNorm=NULL)
{
    if(is.null(neededNorm)){
        outsideOption = TRUE
    }else{
        outsideOption = FALSE
    }
    #
    nbX = length(n)
    nbY = length(m)
    #
    TUs = build_TUs(phi)
    logitM = build_logits(nbX,nbY,sigma=sigma,outsideOption=outsideOption)
    logitW = build_logits(nbY,nbX,sigma=sigma,outsideOption=outsideOption)
    #
    ret = list(n=n,m=m,
               hetG=logitM,hetH=logitW,
               transfers=TUs,
               neededNorm=neededNorm)
    class(ret) = "TU_logit"
    #
    return(ret)
}

solveEquilibrium.TU_logit = ipfp

margxInv.TU_logit <- function(xs, mkt, Mu0ys, sigma=1)
{
    if(is.null(xs)){
        xs = 1:mkt$transfers$nbX
    }
    #
    sqrtMu0ys = sqrt(Mu0ys)
    if(is.null(mkt$neededNorm)){
        b = (mkt$transfers$aux_expphiover2[xs,]^(1/sigma) %*% sqrtMu0ys)/2
        sqrtMux0s = sqrt(mkt$n[xs]+ b*b) - b
    }else{
        sqrtMux0s = mkt$n / c(mkt$transfers$aux_expphiover2[xs,]^(1/sigma) %*% sqrtMu0ys)
    }
    #
    ret = c(sqrtMux0s*sqrtMux0s)
    #
    return(ret)
}

margyInv.TU_logit <- function(ys, mkt, Mux0s, sigma=1)
{
    if(is.null(ys)){
        ys = 1:mkt$transfers$nbY
    }
    #
    sqrtMux0s = sqrt(Mux0s)
    if(is.null(mkt$neededNorm)){
        b = c(sqrtMux0s %*% mkt$transfers$aux_expphiover2[,ys]^(1/sigma))/2
        sqrtMu0ys = sqrt(mkt$m[ys] + b*b) - b
    }else{
        sqrtMu0ys = mkt$m / c(sqrtMux0s %*% mkt$transfers$aux_expphiover2[,ys]^(1/sigma))
    }
    #
    ret = sqrtMu0ys*sqrtMu0ys
    #
    return(ret)
}
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

buildModel_affinity <- function(Xvals, Yvals, n=NULL, m=NULL, noSingles=FALSE)
{
    nbX = dim(Xvals)[1]
    nbY = dim(Yvals)[1]
    #
    dX = dim(Xvals)[2]
    dY = dim(Yvals)[2]
    #
    if(is.null(n)){
        n = rep(1,nbX)
    }
    if(is.null(m)){
        m = rep(1,nbY)
    }
    #
    if(noSingles){
        neededNorm = list(H_edge_logit = function(mux0,mu0y) (mu0y[1]))
    }else{
        neededNorm = NULL
    }
    #
    ret = list(kron=kronecker(Yvals,Xvals),
               nbParams=dX*dY,
               nbX=nbX, nbY=nbY,
               dX=dX, dY=dY,
               n=n, m=m,
               isLinear=TRUE,
               neededNorm=neededNorm)
    class(ret) = "affinity"
    #
    return(ret)
}

parametricMarket.affinity <- function(model, params)
    # params is simply the affinity matrix
{
    phi = matrix(model$kron %*% c(params),nrow=model$nbX)
    #
    ret = build_market_TU_logit(model$n,model$m,phi,
                                neededNorm=model$neededNorm)
    #
    return(ret)
}

dparam.affinity <- function(model, dparams=diag(model$nbParams))
    # params is simply the affinity matrix
{
    dparamsPsi = model$kron %*% dparams
    dparamsG = matrix(0,nrow=0,ncol=dim(dparams)[2])
    dparamsH = matrix(0,nrow=0,ncol=dim(dparams)[2])
    #
    ret = list(dparamsPsi = dparamsPsi,
               dparamsG   = dparamsG,
               dparamsH   = dparamsH)
    #
    return(ret)
}
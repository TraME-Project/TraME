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

buildModel_etu <- function(Xvals, Yvals, n=NULL, m=NULL)
{
    nbX = dim(t(t(Xvals)))[1]
    nbY = dim(t(t(Yvals)))[1]
    #
    dX = dim(t(t(Xvals)))[2]
    dY = dim(t(t(Yvals)))[2]
    #
    eX = matrix(rep(1,nbX),ncol=1)
    eY = matrix(rep(1,nbY),ncol=1)
    #
    diff = abs(kronecker(eY,Xvals)-kronecker(Yvals,eX))
    #
    if(is.null(n)){
        n=rep(1,nbX)
    }
    if(is.null(m)){
        m=rep(1,nbY)
    }
    #
    ret = list(diff=diff,
               nbParams=2*dim(t(t(diff)))[2]+1,
               nbX=nbX, nbY=nbY,
               dX=dX, dY=dY,
               n=n, m=m,
               isLinear=TRUE)
    class(ret) = "etu"
    #
    return(ret)
}

parametricMarket.etu <- function(model, theta)
    # the theta are the parameters for alpha, gamma and tau
{
    theta1 = theta[1:model$dX]
    theta2 = theta[(model$dX+1):(model$dX+model$dY)]
    theta3 = theta[length(theta)]
    #
    alpha = matrix(model$diff %*% theta1,nrow=model$nbX)
    gamma = matrix(model$diff %*% theta2,nrow=model$nbX)
    
    tau = matrix(theta3, model$nbX, model$nbY)
    #
    tr = build_ETUs(alpha, gamma, tau)
    #
    ret = build_market_ITU_logit(model$n,model$m,tr,sigma=1)
    #
    return(ret)
}

dparam.etu <- function(model, dparams=diag(model$nbParams))
    # params is simply the affinity matrix
{
    zero1 = matrix(0,model$nbX*model$nbY,model$dX)
    zero2 = matrix(0,model$nbX*model$nbY,model$dY)
    zero3 = matrix(0,model$nbX*model$nbY,1)
    #
    dparamsPsi = rbind(cbind(model$diff,zero2,zero3),
                       cbind(zero1,model$diff,zero3),
                       cbind(zero1,zero2,rep(1,model$nbX*model$nbY)))
    dparamsG = matrix(0,nrow=1,ncol=dim(dparams)[2])
    dparamsH = matrix(0,nrow=1,ncol=dim(dparams)[2])
    #
    ret = list(dparamsPsi = dparamsPsi,
               dparamsG   = dparamsG,
               dparamsH   = dparamsH)
    #
    return(ret)
}

initparam.etu <- function(model)
{
    ret = list(param=c(rep(0,model$nbParams-1),1),
               lb=NULL,ub=NULL)
    #
    return(ret)
}
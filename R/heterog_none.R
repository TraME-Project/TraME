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

build_none <- function(nbX, nbY)
    # Extreme value type I errors; sigma is the scaling parameter
{
    ret = list(nbX=nbX,nbY=nbY,nbParams=0)
    class(ret) = "none"
    #
    return(ret)
}

Gx.none <- function(heterog, Ux, x)
{
    nbY = length(Ux)
    y = which.max(c(Ux,0))
    #
    mux = rep(0,nbY+1)
    mux[y] = 1
    #
    return(list(valx = max(c(Ux,0),mux = mux)))
}

Gstar.none <- function(heterog, mu)
{
    stop("Gstar not yet defined for no heterogeneity case.")
}

dtheta_NablaGstar.none <- function(het, mu, n, dtheta=NULL, xFirst=TRUE)
{
    return(rep(0,het$nbX*het$nbY))
}

Gbarx.none <- function(heterog, Ubarx, mubarx, x)
{
    nbY0 = length(Ubarx)
    #
    srt = sort(Ubarx,decreasing=TRUE,index.return=TRUE)
    rk = 1
    #
    mux = rep(0,nbY0)
    cumul = mubarx[srt$ix[rk]]
    #
    while((rk<nbY0) & (cumul<1) & (Ubarx[srt$ix[rk]]>0)){
        mux[srt$ix[rk]] = mubarx[srt$ix[rk]]
        rk = rk + 1
        cumul = cumul + mubarx[srt$ix[rk]]
    }
    #
    if(Ubarx[srt$ix[rk]] > 0){
        mux[srt$ix[rk]] = mubarx[srt$ix[rk]] + 1 - cumul
    }
    #
    return(list(valx=sum(mux*Ubarx),
                Ux=rep(NA,nbY0), 
                mux=mux))
}

simul.none <- function(heterog, nbDraws, seed=NULL)
{
    set.seed(seed)
    #
    ret = list(nbX=het$nbX, nbY=het$nbY,
               nbParams = nbDraws * (het$nbY+1) * het$nbX,
               atoms = array(0,dim=c(nbDraws,het$nbY+1, het$nbX)),
               aux_nbDraws=nbDraws,
               xHomogenous=FALSE)
    class(ret) = "empirical"
    #
    return(ret)
}
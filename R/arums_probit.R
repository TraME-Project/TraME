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
# D. McFadden: Modelling the choice of residential location. 
#   In A. Karlquist et. al., editor, Spatial Interaction Theory and Residential Location. 
#   North Holland Pub. Co., 1978. 
# Anderson, S., de Palma, A., and Thisse, J.-F. 
# K. Train: Discrete Choice Models with Simulation. Cambridge.
# A. Galichon, B. Salanie: " Cupid's Invisible Hand: Social Surplus and Identification in Matching Models"
# K. Chiong, A. Galichon, M. Shum: "Duality in dynamic discrete choice models"
#           

build_probit <- function(Covar, outsideOption=TRUE)
    # dim(Covar)=c(nbX,nbOptions,nbOptions), where nbOptions=nbY+1 if outsideOption==T, =nbY else
    # S[x,,] is the correlation matrix between utility shocks associated to nbOptions alternatives
    # there are nbOptions alternatives (when outsideOption==T,the last one is associated to the outside option)
{
    nbX = dim(Covar)[1]
    nbOptions = dim(Covar)[2]
    nbY = nbOptions - ifelse(outsideOption,1,0)
    #
    ret = list(nbX=nbX, nbY=nbY,
               nbParams = nbX*nbOptions*(nbOptions-1)/2,
               Covar=Covar,
               aux_nbOptions=nbOptions,
               outsideOption=outsideOption)
    class(ret) = "probit"
    #
    return(ret)
}

simul.probit <- function(arums, nbDraws, seed=NULL)
{
    set.seed(seed)
    #
    atoms = array(0,dim=c(nbDraws,arums$aux_nbOptions,arums$nbX))
    for(x in 1:arums$nbX){
        E = eigen(arums$Covar[x,,])
        V = E$values 
        Q = E$vectors
        
        SqrtCovar = Q%*%diag(1/sqrt(V))%*%t(Q) # Keith: benefit of this vs chol?
        
        atoms[,,x] = matrix(rnorm(nbDraws*arums$aux_nbOptions),ncol=arums$aux_nbOptions) %*% SqrtCovar 
    }
    #
    ret = list(nbX=arums$nbX, nbY=arums$nbY,
               nbParams=length(atoms),
               atoms=atoms,
               xHomogenous=FALSE,
               aux_nbDraws=nbDraws,
               outsideOption=arums$outsideOption)
    class(ret) = "empirical"
    #
    return(ret)
}

unifCorrelCovMatrices <- function(nbX, nbY, rho, outsideOption=TRUE)
{
    if(outsideOption){
        nbOptions = nbY + 1
        #
        Sig = rho * matrix(1,nbOptions,nbOptions) + (1-rho) * diag(1,nbOptions)
        Sig[,nbOptions] = 0
        Sig[nbOptions,] = 0
        Sig[nbOptions,nbOptions] = 1
        #
        Covar = array(0,c(nbX,nbOptions,nbOptions))
    }else{
        nbOptions = nbY
        Sig = rho * matrix(1,nbOptions,nbOptions) + (1-rho) * diag(1,nbOptions)
    }
    #
    Covar = array(0,c(nbX,nbOptions,nbOptions))
    for(x in 1:nbX){
        Covar[x,,] = Sig
    }
    #
    return(Covar)
}  

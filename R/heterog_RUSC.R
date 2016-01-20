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
# A. Galichon, B. Salanie: " Cupid's Invisible Hand: Social Surplus and Identification in Matching Models"
#           
### This function has the only purpose to launch a Tutorial for the Package TraME

build_RUSC <- function(zeta, outsideOption=TRUE)
    # RUSC errors; sigma is the scaling parameter; 
    # dim(zeta)=c(nbX,nbY+1) and the last column corresponds to singlehood
{
    if(!outsideOption){
        stop("outsideOption=F not implemented yet on RUSC heterogeneity")
    }
    #
    nbX = dim(zeta)[1]
    nbY = dim(zeta)[2]-1
    
    aux_ord = array(0,dim=c(nbX,nbY+1))
    aux_A = array(0,dim=c(nbX,nbY,nbY))
    aux_b = array(0,dim=c(nbX,nbY))
    aux_c = rep(0,nbX)
    #
    for(x in 1:nbX){
        z = zeta[x,1:nbY]
        z0 = zeta[x,nbY+1]
        #
        maxzz = pmax(z %*% matrix(1,1,nbY), matrix(1,nbY,1) %*% t(z))
        maxzz0 = pmax(z,z0)
        maxzz0Mat = matrix(maxzz0,ncol=1) %*% matrix(1,1,nbY)
        
        Ax =  maxzz0Mat + t(maxzz0Mat) - maxzz - z0
        #
        aux_A[x,,] = Ax
        #Ainv[x,,] = solve(Ax)
        aux_b[x,] = z0 - maxzz0
        aux_c[x]  = -z0/2
        aux_ord[x,] = order(zeta[x,])
    }
    #
    ret = list(nbX=nbX, nbY=nbY,
               nbParams=length(zeta), zeta=zeta,
               aux_A=aux_A, aux_b=aux_b,
               aux_c=aux_c, aux_ord=aux_ord,
               outsideOption=outsideOption)
    class(ret) = "RUSC"
    #
    return(ret)
}

Gx.RUSC <- function(het, Ux, x)
{
    M = length(Ux) + 1
    #
    muxtilde = rep(0,M)
    Uxtilde  = c(Ux,0)
    #
    for(i in 1:M){
        y = het$aux_ord[x,i]
        runmax = 0
        #
        j = 1
        while(j < i){
            z = het$aux_ord[x,j]
            runmax = max(runmax, (Uxtilde[y]-Uxtilde[z])/(het$zeta[x,z]-het$zeta[x,y]))
            j = j + 1
        }
        #
        runmin = 1
        j = M
        while(j > i){
            z = het$aux_ord[x,j]
            runmin = min(runmin, (Uxtilde[y]-Uxtilde[z])/(het$zeta[x,z]-het$zeta[x,y]))
            j = j - 1 
        }
        #
        muxtilde[y] = max(runmin-runmax,0)
    }
    #
    mux = muxtilde[1:M-1]
    #
    valx = sum(mux*(Ux-het$aux_b[x,])) - matrix(mux,nrow=1)%*%het$aux_A[x,,]%*%matrix(mux,ncol=1)/2 - het$aux_c[x]
    ret = list(valx = valx, mux  = mux)
    #
    return(ret)
}

Gstarx.RUSC <- function(het, mux, x)
{
    Amu = c(het$aux_A[x,,] %*% matrix(mux,ncol=1))
    #
    ret = list(valx = sum(mux*Amu)/2 + sum(mux*het$aux_b[x,]) + het$aux_c[x],
               Ux   = Amu + het$aux_b[x,])
    #
    return(ret)
}

Gbarx.RUSC <- function (het, Ubarx, mubarx, x)
{
    M = length(Ubarx) + 1
    #
    obj = c(het$aux_b[x,] - Ubarx,0)
    A = matrix(1,1,M)
    rhs = c(1)
    
    Q = matrix(0,M,M)
    Q[1:M-1,1:M-1] = het$aux_A[x,,] / 2
    
    lb = rep(0,M)
    ub = c(mubarx,1)
    #
    result = genericLP(obj=obj,A=A,modelsense="min",rhs=rhs,sense="=",Q=Q,lb=lb,ub=ub)
    #
    mux = result$solution[1:M-1]
    Amu = c(het$aux_A[x,,] %*% matrix(mux,ncol=1))
    Ux  = Amu + het$aux_b[x,]
    #
    ret = list(valx = -result$objval - het$aux_c[x],
               mux=mux, Ux=Ux)
    #
    return(ret)
}

simul.RUSC <- function(het, nbDraws, seed=NULL)
{  
    set.seed(seed)
    #
    atoms = array(0,dim=c(nbDraws,het$nbY+1,het$nbX))
    for(x in 1:het$nbX){
        atoms[,,x] =  matrix(runif(nbDraws),ncol=1) %*% matrix(het$zeta[x,],nrow=1)
    }
    #
    ret = list(nbX=het$nbX, nbY=het$nbY,
               nbParams=length(atoms),
               atoms=atoms,
               aux_nbDraws=nbDraws,
               xHomogenous=FALSE,
               outsideOption=het$outsideOption)
    class(ret) = "empirical"
    #
    return(ret)
}

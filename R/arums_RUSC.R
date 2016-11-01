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
        stop("outsideOption=F not implemented yet on RUSC arums")
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

Gx.RUSC <- function(arums, Ux, x)
{
    nbAlt = length(Ux) + 1
    #
    muxtilde = rep(0,nbAlt)
    Uxtilde  = c(Ux,0)
    #
    for(i in 1:nbAlt){
        y = arums$aux_ord[x,i]
        runmax = 0
        #
        j = 1
        while(j < i){
            z = arums$aux_ord[x,j]
            runmax = max(runmax, (Uxtilde[y]-Uxtilde[z])/(arums$zeta[x,z]-arums$zeta[x,y]))
            j = j + 1
        }
        #
        runmin = 1
        j = nbAlt
        while(j > i){
            z = arums$aux_ord[x,j]
            runmin = min(runmin, (Uxtilde[y]-Uxtilde[z])/(arums$zeta[x,z]-arums$zeta[x,y]))
            j = j - 1 
        }
        #
        muxtilde[y] = max(runmin-runmax,0)
    }
    #
    mux = muxtilde[1:nbAlt-1]
    #
    valx = sum(mux*(Ux-arums$aux_b[x,])) - matrix(mux,nrow=1)%*%arums$aux_A[x,,]%*%matrix(mux,ncol=1)/2 - arums$aux_c[x]
    ret = list(valx = valx, mux  = mux)
    #
    return(ret)
}

Gstarx.RUSC <- function(arums, mux, x)
{
    Amu = c(arums$aux_A[x,,] %*% matrix(mux,ncol=1))
    #
    ret = list(valx = sum(mux*Amu)/2 + sum(mux*arums$aux_b[x,]) + arums$aux_c[x],
               Ux   = Amu + arums$aux_b[x,])
    #
    return(ret)
}

Gbarx.RUSC <- function (arums, Ubarx, mubarx, x)
{
    nbAlt = length(Ubarx) + 1
    #
    obj = c(arums$aux_b[x,] - Ubarx,0)
    A = matrix(1,1,nbAlt)
    rhs = c(1)
    
    Q = matrix(0,nbAlt,nbAlt)
    Q[1:(nbAlt-1),1:(nbAlt-1)] = arums$aux_A[x,,] / 2
    
    lb = rep(0,nbAlt)
    ub = c(mubarx,1)
    #
    result = genericLP(obj=obj,A=A,modelsense="min",rhs=rhs,sense="=",Q=Q,lb=lb,ub=ub)
    #
    mux = result$solution[1:(nbAlt-1)]
    Amu = c(arums$aux_A[x,,] %*% matrix(mux,ncol=1))
    Ux  = Amu + arums$aux_b[x,]
    #
    ret = list(valx = -result$objval - arums$aux_c[x],
               mux=mux, Ux=Ux)
    #
    return(ret)
}

simul.RUSC <- function(arums, nbDraws, seed=NULL)
{  
    set.seed(seed)
    #
    atoms = array(0,dim=c(nbDraws,arums$nbY+1,arums$nbX))
    for(x in 1:arums$nbX){
        atoms[,,x] =  matrix(runif(nbDraws),ncol=1) %*% matrix(arums$zeta[x,],nrow=1)
    }
    #
    ret = list(nbX=arums$nbX, nbY=arums$nbY,
               nbParams=length(atoms),
               atoms=atoms,
               aux_nbDraws=nbDraws,
               xHomogenous=FALSE,
               outsideOption=arums$outsideOption)
    class(ret) = "empirical"
    #
    return(ret)
}

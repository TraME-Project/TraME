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

build_empirical <- function(nbX, nbY, atoms, outsideOption=TRUE)
    # Empirical distribution; 
    # xHomogenous = T if all Gx are the same, in which case the empirical counterpart will be the same
    # if xHomogenous=T, then dim(atoms)=(aux_nbDraws,nbY+1)  
    # dim(atoms) = c(aux_nbDraws,nbY+1,nbX)
{
    dims = dim(atoms)
    #
    if(length(dims)==2){
        xHomogenous = TRUE
    }else if(length(dims)==3){
        xHomogenous = FALSE
    }else{
        stop("array atoms does not have the right dimensions")
    }
    #
    ret = list(nbX=nbX,nbY=nbY,
               nbParams=length(atoms),
               atoms=atoms,
               aux_nbDraws=dim(atoms)[1],
               xHomogenous=xHomogenous,
               outsideOption=outsideOption)
    class(ret) = "empirical"
    #
    return(ret)
}

Gx.empirical <- function(arums, Ux, x)
{
    if(arums$outsideOption){
        Uxs = c(Ux,0)
    }else{
        Uxs = Ux
    }
    if(arums$xHomogenous){
        Utilde = matrix(1,nrow=arums$aux_nbDraws,ncol=1) %*% matrix(Uxs,nrow=1) + arums$atoms
    }else{
        Utilde = matrix(1,nrow=arums$aux_nbDraws,ncol=1) %*% matrix(Uxs,nrow=1) + arums$atoms[,,x]
    }
    #
    argmaxs = apply(Utilde, 1, which.max)
    #
    thesum = 0
    for(t in 1:arums$aux_nbDraws){
        thesum = thesum + Utilde[t,argmaxs[t]]
    }
    #
    ret = list(valx = thesum/arums$aux_nbDraws,
               mux = apply(matrix(c(1:arums$nbY),ncol=1),1,function(x) sum(argmaxs==x)) / arums$aux_nbDraws)
    #
    return(ret)
}

Gstarx.empirical <- function(arums, mux, x)
{
    nbOptions = ifelse(arums$outsideOption,arums$nbY+1,arums$nbY)
    if(arums$xHomogenous){
        Phi = arums$atoms
    }else{
        Phi = arums$atoms[,,x]
    }
    #
    p = rep(1,arums$aux_nbDraws) / arums$aux_nbDraws
    q = 0
    if(arums$outsideOption){
        q = c(mux, 1-sum(mux))
    }else{
        q = mux
    }
    #
    obj = c(Phi)
    
    A1 = Matrix::kronecker(matrix(1,1,nbOptions),sparseMatrix(1:arums$aux_nbDraws,1:arums$aux_nbDraws))
    A2 = Matrix::kronecker(sparseMatrix(1:nbOptions,1:nbOptions),matrix(1,1,arums$aux_nbDraws))
    
    A = rbind2(A1,A2)
    d = c(p,q)
    #
    result = genericLP(obj=obj,A=A,modelsense="max",rhs=d,sense="=")
    #
    # pi = matrix(result$solution,nrow=arums$aux_nbDraws)
    u = result$pi[1:arums$aux_nbDraws]
    
    if(arums$outsideOption){
        temp = result$pi[(arums$aux_nbDraws+1):(arums$aux_nbDraws+arums$nbY+1)]
        Uoptx = -temp[1:arums$nbY]+temp[arums$nbY+1]
    }else{
        Uoptx = -result$pi[(arums$aux_nbDraws+1):(arums$aux_nbDraws+arums$nbY)] - sum(p*u)
    }
    
    valx = -result$objval
    #
    return(list(valx=valx,Ux=Uoptx))
}

Gbarx.empirical <- function(arums, Ubarx, mubarx, x)
{
    if(!arums$outsideOption){
        stop("Gbarx not implemented for empirical with outsideOption=F")
    }
    if(arums$xHomogenous){
        Phi = arums$atoms
    }else{
        Phi = arums$atoms[,,x]
    }
    #
    A1 = diag(1,arums$nbY)
    A2 = sparseMatrix(i=1,j=1,x=0,dims=c(arums$nbY,arums$aux_nbDraws))
    A3 = Matrix::kronecker(sparseMatrix(i=1:arums$nbY,j=1:arums$nbY),matrix(1,arums$aux_nbDraws,1))
    A4 = - Matrix::kronecker(matrix(1,arums$nbY,1),sparseMatrix(i=1:arums$aux_nbDraws,j=1:arums$aux_nbDraws))
    A5 = sparseMatrix(i=1,j=1,x=0,dims=c(arums$aux_nbDraws,arums$nbY))
    A6 = - sparseMatrix(i=1:arums$aux_nbDraws,j=1:arums$aux_nbDraws)
    # should this be rbind or rbind2?
    A  = rbind2(cbind2(A1,A2),rbind2(cbind2(A3,A4),cbind2(A5,A6)))
    #
    d1 = matrix(t(Ubarx),ncol=1)
    d2 = matrix(-Phi, ncol=1)
    d  = rbind2(d1,d2)
    #
    c1 = matrix(t(mubarx),ncol=1)
    c2 = matrix(-1/arums$aux_nbDraws,arums$aux_nbDraws,1)
    c  = rbind(c1,c2)
    #
    z1 = matrix(0,arums$nbY,1)
    z2 = matrix(apply(Phi, 1, function(x) max(x)),ncol=1)
    z_init = rbind(z1,z2)
    #
    # ------------- this program really computes Glowerbarx but the result is adapted
    result = genericLP(obj=c,A=A,modelsense="max",rhs=d,sense="<",start=z_init)
    #
    Uoptx = c(matrix(result$solution[1:arums$nbY],nrow=1))
    deltamux = matrix(result$pi[1:arums$nbY],nrow=1)
    mux = c(mubarx - deltamux)
    valx = sum(mubarx*Ubarx) - result$objval
    #
    ret = list(valx=valx,
               mux=mux,
               Ux=Uoptx)
    #
    return(ret)
}

simul.empirical <- function(arums, nbDraws, seed=NULL) 
{
    set.seed(seed)
    #
    nbOptions = ifelse(arums$outsideOption,arums$nbY+1,arums$nbY)
    if(arums$xHomogenous){ 
        bootstraped_atoms = vector(0,dim=c(nbDraws,nbOptions))
        for(y in 1:nbOptions){
            bootstraped_atoms[,y] = sample(arums$atoms[,y],nbDraws,replace=TRUE)
        }
    }else{
        bootstraped_atoms = vector(0,dim=c(nbDraws,nbOptions,arums$nbX))
        for(x in 1:nbX){
            for(y in 1:nbOptions){
                bootstraped_atoms[,y,x] = sample(arums$atoms[,y,x],nbDraws,replace=TRUE)
            }
        }
    }
    #
    ret = list(nbX=arums$nbX, nbY=arums$nbY,
               nbParams=length(bootstrapped_atoms),
               atoms=bootstraped_atoms,
               aux_nbDraws=nbDraws,
               xHomogenous=arums$xHomogenous,
               outsideOption=arums$outsideOption)
    class(ret) = "empirical"
    #
    return(ret)  
}

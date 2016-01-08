################################################################################
##
##   Copyright (C) 2015 - 2016 Alfred Galichon
##
##   This file is part of the R package TraME.
##
##   The R package TraME free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 2 of the License, or
##   (at your option) any later version.
##
##   The R package BMR is distributed in the hope that it will be useful,
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

Gx.empirical <- function(het, Ux, x)
{
    if(het$outsideOption){
        Uxs=c(Ux,0)
    }else{
        Uxs = Ux
    }
    if(het$xHomogenous){
        Utilde = matrix(1,nrow=het$aux_nbDraws,ncol=1) %*% matrix(Uxs,nrow=1) + het$atoms
    }else{
        Utilde = matrix(1,nrow=het$aux_nbDraws,ncol=1) %*% matrix(Uxs,nrow=1) + het$atoms[,,x]
    }
    #
    argmaxs = apply(Utilde, 1, which.max)
    #
    thesum = 0
    for(t in 1:het$aux_nbDraws){
        thesum = thesum + Utilde[t,argmaxs[t]]
    }
    #
    ret = list(valx = thesum/het$aux_nbDraws,
               mux = apply(matrix(c(1:het$nbY),ncol=1),1,function(x) sum(argmaxs==x)) / het$aux_nbDraws)
    #
    return(ret)
}

Gstarx.empirical <- function(het, mux, x)
{
    nbOptions = ifelse(het$outsideOption,het$nbY+1,het$nbY)
    if(het$xHomogenous){
        Phi = het$atoms
    }else{
        Phi = het$atoms[,,x]
    }
    #
    p = rep(1,het$aux_nbDraws) / het$aux_nbDraws
    q = 0
    if(het$outsideOption){
        q = c(mux, 1-sum(mux))
    }else{
        q = mux
    }
    #
    obj = c(Phi)
    
    A1 = kronecker(matrix(1,1,nbOptions),sparseMatrix(1:het$aux_nbDraws,1:het$aux_nbDraws))
    A2 = kronecker(sparseMatrix(1:nbOptions,1:nbOptions),matrix(1,1,het$aux_nbDraws))
    
    A = rbind2(A1,A2)
    d = c(p,q)
    #
    gurobiModel = list(A=A,obj=obj,
                       modelsense="max",
                       rhs=d,
                       sense="=")
    result = gurobi(gurobiModel, params=list(OutputFlag=0))
    #
    if(result$status=="OPTIMAL"){
        # pi = matrix(result$x,nrow=het$aux_nbDraws)
        u = result$pi[1:het$aux_nbDraws]
        #
        if(het$outsideOption){
            temp = result$pi[(het$aux_nbDraws+1):(het$aux_nbDraws+het$nbY+1)]
            Uoptx = -temp[1:het$nbY]+temp[het$nbY+1]
        }else{ 
            Uoptx = -result$pi[(het$aux_nbDraws+1):(het$aux_nbDraws+het$nbY)] - sum(p*u)
        }
        #
        valx = -result$objval
    }else{
        stop("optimization problem with Gurobi")
    }
    #
    return(list(valx=valx,Ux=Uoptx))
}

Gbarx.empirical <- function(het, Ubarx, mubarx, x)
{
    if(!het$outsideOption){
        stop("Gbarx not implemented for empirical with outsideOption=F")
    }
    if(het$xHomogenous){
        Phi = het$atoms
    }else{
        Phi = het$atoms[,,x]
    }
    #
    A1 = diag(1,het$nbY)
    A2 = sparseMatrix(i=1,j=1,x=0,dims=c(het$nbY,het$aux_nbDraws))
    A3 = kronecker(sparseMatrix(i=1:het$nbY,j=1:het$nbY),matrix(1,het$aux_nbDraws,1))
    A4 = - kronecker(matrix(1,het$nbY,1),sparseMatrix(i=1:het$aux_nbDraws,j=1:het$aux_nbDraws))
    A5 = sparseMatrix(i=1,j=1,x=0,dims=c(het$aux_nbDraws,het$nbY))
    A6 = - sparseMatrix(i=1:het$aux_nbDraws,j=1:het$aux_nbDraws)
    # should this be rbind or rbind2?
    A  = rbind2(cbind2(A1,A2),rbind2(cbind2(A3,A4),cbind2(A5,A6)))
    #
    d1 = matrix(t(Ubarx),ncol=1)
    d2 = matrix(-Phi, ncol=1)
    d  = rbind2(d1,d2)
    #
    c1 = matrix(t(mubarx),ncol=1)
    c2 = matrix(-1/het$aux_nbDraws,het$aux_nbDraws,1)
    c  = rbind(c1,c2)
    #
    z1 = matrix(0,het$nbY,1)
    z2 = matrix(apply(Phi, 1, function(x) max(x)),ncol=1)
    z_init = rbind(z1,z2)
    #
    # ------------- this program really computes Glowerbarx but the result is adapted
    gurobiModel = list(A=A,obj=c,
                       modelsense="max",
                       rhs=d,
                       sense="<",
                       start=z_init)
    result = gurobi(gurobiModel, params=list(OutputFlag=0))
    #
    if(result$status=="OPTIMAL") {
        Uoptx = c(matrix(result$x[1:het$nbY],nrow=1))
        deltamux = matrix(result$pi[1:het$nbY],nrow=1)
        mux = c(mubarx - deltamux)
        valx = sum(mubarx*Ubarx) - result$objval
    }
    else{
        stop("optimization problem with Gurobi")
    }
    #
    ret = list(valx=valx,
               Ux=Uoptx,
               mux=mux)
    #
    return(ret)
}

simul.empirical <- function(het, nbDraws, seed=NULL) 
{
    set.seed(seed)
    #
    nbOptions = ifelse(het$outsideOption,het$nbY+1,het$nbY)
    if(het$xHomogenous){ 
        bootstraped_atoms = vector(0,dim=c(nbDraws,nbOptions))
        for(y in 1:nbOptions){
            bootstraped_atoms[,y] = sample(het$atoms[,y],nbDraws,replace=TRUE)
        }
    }else{
        bootstraped_atoms = vector(0,dim=c(nbDraws,nbOptions,het$nbX))
        for(x in 1:nbX){
            for(y in 1:nbOptions){
                bootstraped_atoms[,y,x] = sample(het$atoms[,y,x],nbDraws,replace=TRUE)
            }
        }
    }
    #
    ret = list(nbX=het$nbX, nbY=het$nbY,
               nbParams=length(bootstrapped_atoms),
               atoms=bootstraped_atoms,
               aux_nbDraws=nbDraws,
               xHomogenous=het$xHomogenous,
               outsideOption=het$outsideOption)
    class(ret) = "empirical"
    #
    return(ret)  
}

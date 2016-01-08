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

build_RSC <- function(zeta, aux_cdf_eps, aux_quant_eps, aux_pot_eps=NULL, aux_pdf_eps=NULL, outsideOption=TRUE)
    # RSC errors; 
    # dim(zeta)=c(nbX,nbY+1) and the last column corresponds to singlehood
    # aux_cdf_eps is the cdf ; aux_quant_eps is the quantile (option) ; 
    # aux_pot_eps is a primitive of the quantile
{
    if(!outsideOption){
        stop("outsideOption=F not implemented yet on RSC heterogeneity")
    }
    if(is.null(aux_pdf_eps)==TRUE){
        aux_pdf_eps <- function (xs) sapply( xs, function (x) (grad(aux_cdf_eps,x) ))
    }
    if(is.null(aux_pot_eps)==TRUE){
        aux_pot_eps <- function (xs) sapply( xs, function (x) (integrate(aux_quant_eps,0,x)$value))
    }
    #
    nbX = dim(zeta)[1]
    nbY = dim(zeta)[2] - 1
    #
    aux_ord = array(0,dim=c(nbX,nbY+1))
    
    D    = diag(nbY+1) - rbind(rep(0,nbY+1), cbind(diag(nbY),rep(0,nbY)))
    Dinv = solve(D)
    
    ts = array(0,c(nbX,nbY+1,nbY+1))
    N  = cbind(diag(nbY),rep(-1,nbY))
    
    aux_Influence_lhs = c()
    aux_Influence_rhs = c()
    aux_Psigma        = c()
    aux_DinvPsigma    = c()
    #
    for(x in 1:nbX){
        aux_ordx =  order(zeta[x,])
        aux_ord[x,] = aux_ordx
        
        Psigmax = as(aux_ordx,"pMatrix")
        
        aux_Influence_lhs[[x]] = N %*% Matrix::crossprod(Psigmax,Dinv) # Keith: Matrix:: is required to get around a weird bug
        aux_Influence_rhs[[x]] = D %*% Psigmax
        aux_Psigma[[x]]        = diag(nbY+1) %*% Psigmax
        aux_DinvPsigma[[x]]    = Dinv %*% Psigmax
    }
    #
    ret = list(nbX=nbX, nbY=nbY,
               nbParams=length(zeta), zeta=zeta, 
               aux_ord=aux_ord,
               aux_Influence_lhs=aux_Influence_lhs,
               aux_Influence_rhs=aux_Influence_rhs,
               aux_DinvPsigma=aux_DinvPsigma, aux_Psigma=aux_Psigma,  
               aux_cdf_eps=aux_cdf_eps, aux_quant_eps=aux_quant_eps,
               aux_pot_eps=aux_pot_eps, aux_pdf_eps=aux_pdf_eps,
               outsideOption=TRUE)
    class(ret) = "RSC"
    #
    return(ret)
}

build_RSCunif <- function(zeta, alpha, beta)
    # RSC beta errors; 
    # dim(zeta)=c(nbX,nbY+1) and the last column corresponds to singlehood
    # epsilon is a uniform distribution
{
    warning("only use RSCunif for testing purposes -- use RUSC otherwise")
    #
    ret = build_RSC(zeta,aux_cdf_eps=punif,aux_quant_eps=qunif,aux_pdf_eps=dunif)
    #
    return(ret)
}

build_RSCbeta <- function(zeta, alpha, beta)
    # RSC beta errors; 
    # dim(zeta)=c(nbX,nbY+1) and the last column corresponds to singlehood
    # epsilon is a beta(alpha,beta) distribution
{
    aux_cdf_eps   <- function(x) (pbeta(x,alpha,beta))
    aux_quant_eps <- function(x) (qbeta(x,alpha,beta))
    aux_pdf_eps   <- function(x) (dbeta(x,alpha,beta))
    #
    ret = build_RSC(zeta,aux_cdf_eps=aux_cdf_eps,
                    aux_quant_eps=aux_quant_eps,
                    aux_pdf_eps=aux_pdf_eps)
    #
    return(ret)
}

# build_RSCnorm <- function( zeta,mean = 0,sd = 1 )
#   # RSC normal heterogeneity; 
#   # dim(zeta)=c(nbX,nbY+1) and the last column corresponds to singlehood
#   # epsilon is a N(mean,sd) distribution
# {
#   
#   if ((mean==0) && (sd==1) )
#   {
#     aux_pot_eps <- function (xs) sapply( xs, function (x) (integrate(qnorm,1e-6,x)$value ))
#   return(build_RSC(zeta,aux_cdf_eps=pnorm,aux_quant_eps=qnorm,aux_pdf_eps=dnorm,aux_pot_eps=aux_pot_eps))
#   }
#   else
#   {
#     aux_cdf_eps <- function(x) (pnorm(x,mean,sd))
#     aux_quant_eps <- function(x) (qnorm(x,mean,sd))
#     aux_pdf_eps <- function(x) (dnorm(x,mean,sd))
#     aux_pot_eps <- function (xs) sapply( xs, function (x) (integrate(aux_quant_eps,1e-6,x)$value ))
#     return(build_RSC(zeta,aux_cdf_eps=aux_cdf_eps,aux_quant_eps=aux_quant_eps,aux_pdf_eps=aux_pdf_eps,aux_pot_eps=aux_pot_eps))
#   }
# }

Gx.RSC <- function(het, Ux, x)
{
    M =  het$nbY + 1
    #
    muxtilde = rep(0,M)
    Uxtilde = c(Ux,0)
    #
    valx = 0
    Eepssofar = 0
    cumulmusofar = 0
    #
    for(i in 1:M){
        y = het$aux_ord[x,i]
        runmax = het$aux_quant(0)
        #
        j = 1
        while(j < i){
            z = het$aux_ord[x,j]
            
            if(het$zeta[x,z] != het$zeta[x,y]){
                runmax = max(runmax, (Uxtilde[y]-Uxtilde[z])/(het$zeta[x,z]-het$zeta[x,y]))
            }else if(Uxtilde[y] < Uxtilde[z]){
                runmax=Inf
            }
            
            j = j + 1
        }
        #
        runmin = het$aux_quant(1)
        j = M
        while(j > i){
            z = het$aux_ord[x,j]
            runmin = min(runmin, (Uxtilde[y]-Uxtilde[z])/(het$zeta[x,z]-het$zeta[x,y]))
            j = j - 1 
        }
        #print(c("y",y,"low",round(runmax,2),"high",round(runmin,2)))
        if(runmin > runmax){
            muxtildey = max(het$aux_cdf_eps(runmin)-het$aux_cdf_eps(runmax),0)
        }else{
            muxtildey = 0
        }
        #
        muxtilde[y] = muxtildey
        #
        if(muxtildey > 0){ 
            cumulmusofar = cumulmusofar + muxtildey
            ey = het$aux_quant_eps(cumulmusofar)
            EepssofarNext = ey*cumulmusofar - het$aux_pot_eps(cumulmusofar)
            valx = valx + muxtildey*Uxtilde[y] + het$zeta[x,y]*(EepssofarNext - Eepssofar) 
        }
    }  
    #
    mux = muxtilde[1:M-1]
    #
    ret = list(valx = sum(mux*Ux) - Gstarx.RSC(het,mux,x)$valx,
               mux  = mux)
    #
    return(ret)
}

Gstarx.RSC <- function(het, mux, x)
{
    tsfull = het$aux_DinvPsigma[[x]] %*% c(mux,1-sum(mux))
    ts = tsfull[1:het$nbY]
    
    pots = het$aux_pot_eps(c(0,tsfull))
    diffpots = pots[2:(het$nbY+2)] - pots[1:(het$nbY+1)]
    #
    valx = -sum( (het$aux_Psigma[[x]] %*% het$zeta[x,]) * diffpots )
    
    e = diag( c(0,het$aux_quant_eps(ts)) )
    Ux = -c( het$aux_Influence_lhs[[x]] %*% e %*% het$aux_Influence_rhs[[x]] %*% het$zeta[x,] )
    #
    ret = list(valx = valx, Ux = Ux)
    #
    return(ret)
    
}

D2Gstar.RSC <- function(het, mu, n, xFirst=TRUE)
{ 
    mux0 = n - apply(mu,1,sum)
    if(xFirst){
        ders = array(0,dim=c(het$nbX,het$nbY,het$nbX,het$nbY))
    }else{
        ders = array(0,dim=c(het$nbY,het$nbX,het$nbY,het$nbX))
    }
    #
    for(x in 1:het$nbX){
        C = -t( (c(het$aux_Influence_rhs[[x]] %*% het$zeta[x,]) ) * t(het$aux_Influence_lhs[[x]]))
        
        tsfull = het$aux_DinvPsigma[[x]] %*% c(mu[x,],mux0[x])/n[x] 
        erestr = het$aux_quant_eps(tsfull)[1:het$nbY]
        
        d_mue_temp = rbind(0,diag( c(1/het$aux_pdf_eps(erestr)) ) %*% (het$aux_DinvPsigma[[x]][1:het$nbY,]))
        d_mue = d_mue_temp[,1:het$nbY] - d_mue_temp[,het$nbY+1]
        #
        if(xFirst){
            ders[x,,x,] = (C %*% d_mue)/n[x]
        }else{
            ders[,x,,x] = (C %*% d_mue)/n[x]
        }
    }
    #
    ret = array(ders,dim=c(het$nbX*het$nbY,het$nbX*het$nbY))
    return(ret)
}

dtheta_NablaGstar.RSC <- function(het, mu, n, dtheta=diag(het$nbParams), xFirst=TRUE)
{
    if(length(dtheta)==0){
        return(matrix(0,nrow=het$nbX*het$nbY,ncol=0))
    }
    #
    nbDirs = length(dtheta) %/% (het$nbX * het$nbX * (het$nbY+1))
    dthetamat = array(dtheta,dim=c(het$nbX,het$nbY+1,het$nbX,nbDirs))
    #
    mux0 = n - apply(mu,1,sum)
    if(xFirst){
        ders = array(0,dim=c(het$nbX,het$nbY,het$nbX,nbDirs))
    }else{
        ders = array(0,dim=c(het$nbY,het$nbX,het$nbX,nbDirs))
    }
    #
    for (x in 1:het$nbX){
        tsfull = het$aux_DinvPsigma[[x]] %*% c(mu[x,],n[x]-sum(mu[x,]))/n[x] 
        e = diag( c(0,het$aux_quant_eps(tsfull[1:het$nbY])) )
        
        if(xFirst){
            ders[x,,x,] = -het$aux_Influence_lhs[[x]] %*% e %*% het$aux_Influence_rhs[[x]] %*% dthetamat[x,,x,]
        }else{
            ders[,x,x,] = -het$aux_Influence_lhs[[x]] %*% e %*% het$aux_Influence_rhs[[x]] %*% dthetamat[x,,x,]
        }
    }
    #
    return(array(ders,dim=c(het$nbX*het$nbY,het$nbX*nbDirs)))
}

simul.RSC <- function(het, nbDraws, seed=NULL)
{  
    set.seed(seed)
    #
    atoms = array(0,dim=c(nbDraws,het$nbY+1,het$nbX))
    for(x in 1:het$nbX){
        atoms[,,x] = matrix(het$aux_quant_eps(runif(nbDraws)),ncol=1) %*% matrix(het$zeta[x,],nrow=1)
    }
    #
    ret = list(nbX=het$nbX, nbY = het$nbY,
               nbParams=length(atoms),
               atoms=atoms,
               xHomogenous=FALSE,
               aux_nbDraws=nbDraws,
               outsideOption=het$outsideOption)
    class(ret) = "empirical"
    #
    return(ret)
}
############################################
############################################
############################################
############################################
#####  Alternative Gstar.RSC function ######
############################################
# Gstarx.RSC <- function(het,mux,x )
# {
#   q = c(mux, (1-sum(mux)))
#   aux_ord = het$aux_ord[x,]
#   zeta = het$zeta[x,]
#   M = het$nbY+1
#   v = rep(0,M) 
#   t = q[aux_ord[1]]
#   eps = het$aux_quant_eps(t)
#   v[aux_ord[1]] = zeta[aux_ord[1]]*eps
#   valx = zeta[aux_ord[1]] * het$aux_pot_eps(q[aux_ord[1]])
#   for (j in 2:M)
#   {
#     valx = valx + zeta[aux_ord[j]] * (het$aux_pot_eps(t+q[aux_ord[j]])-het$aux_pot_eps(t)) 
#     v[aux_ord[j]] = v[aux_ord[j-1]] + het$aux_quant_eps(t)* (zeta[aux_ord[j]]-zeta[aux_ord[j-1]])
#     t = t + q[aux_ord[j]]
#   }
#   
#   
#   Ux = v[M] - v[1:(M-1)]  # CHANGE HERE TO BE VERIFIED
#   
#   return(list(valx= -valx,
#               Ux=  Ux ))
# }


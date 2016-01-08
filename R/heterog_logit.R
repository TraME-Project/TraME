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
# D. McFadden: Modelling the choice of residential location. 
#   In A. Karlquist et. al., editor, Spatial Interaction Theory and Residential Location. 
#   North Holland Pub. Co., 1978. 
# Anderson, S., de Palma, A., and Thisse, J.-F. 
# K. Train: Discrete Choice Models with Simulation. Cambridge.
# A. Galichon, B. Salanie: " Cupid's Invisible Hand: Social Surplus and Identification in Matching Models"
# K. Chiong, A. Galichon, M. Shum: "Duality in dynamic discrete choice models"
#           
#

build_logits <- function(nbX, nbY, sigma=1, outsideOption=TRUE)
    # Extreme value type I errors; sigma is the scaling parameter
{
    ret = list(nbX=nbX,nbY=nbY,
               nbParams=1,
               sigma=sigma,
               outsideOption=outsideOption)
    class(ret) = "logit"
    #
    return(ret)
}

G.logit <- function(het, U, n)
{
    expU = exp(U/het$sigma)
    #
    if(het$outsideOption){
        denom = 1 + apply(expU,1,sum)
    }else{
        denom = apply(expU,1,sum)
    }
    #
    ret = list(val = het$sigma*sum(n*log(denom)),
               mu  = (n /denom)*expU)
    #
    return(ret)
}

Gstar.logit <- function(het, mu, n)
{
    mux0 <- ret <- 0
    #
    if(het$outsideOption){
        mux0 = n - apply(mu,1,sum)  
        ret = list(val = het$sigma*(sum(mu*log(mu/n)) + sum(mux0*log(mux0/n))),
                   U   = het$sigma*log(mu/mux0))
    }else{
        ret = list(val = het$sigma*(sum(mu*log(mu/n))),
                   U   = het$sigma*log(mu/n))
    }
    #
    return(ret)
}

Gstarx.logit <- function(het, mux, x)
{
    mu0 <- ret <- 0
    #
    if(het$outsideOption){
        mu0 = 1 - sum(mux)
        ret = list(valx = het$sigma*(mu0*log(mu0) + sum(mux*log(mux))),
                   Ux   = het$sigma*log(mux/mu0))
    }else{
        ret = list(valx = het$sigma*(sum(mux*log(mux))),
                   Ux   = het$sigma*log(mux))
    }
    #
    return(ret)
}

D2G.logit <- function(het, U, n, xFirst=TRUE)
{ 
    # NOTE; the formula is the same regardless of whether outsideOption == TRUE or FALSE
    muxy = G(het,U,n)$mu
    H = matrix(0,het$nbX*het$nbY,het$nbX*het$nbY)
    #
    for(x in 1:het$nbX){
        for(y in 1:het$nbY){
            for(yprime in 1:het$nbY){
                if(xFirst){
                    H[x+het$nbX*(y-1),x+het$nbX*(yprime-1)] = ifelse(y==yprime,
                                                                     muxy[x,y]*(1-muxy[x,y]/n[x])/het$sigma,
                                                                     -muxy[x,y]*muxy[x,yprime]/(n[x]*het$sigma))
                }else{
                    H[(x-1)*het$nbY+y,(x-1)*het$nbY+yprime] = ifelse(y==yprime,
                                                                     muxy[x,y]*(1-muxy[x,y]/n[x])/het$sigma,
                                                                     -muxy[x,y]*muxy[x,yprime]/(n[x]*het$sigma))
                }
            }
        }
    }
    #
    return(H)
}

D2Gstar.logit <- function(het, mu, n, xFirst=TRUE) # Keith: changed layout a lot
{ 
    mux0 = n - apply(mu,1,sum)
    #
    oneovermux0 = 1/mux0 
    oneovermuxy = 1/mu
    #
    #H = Diagonal(het$nbX*het$nbY,0)
    H = matrix(0,het$nbX*het$nbY,het$nbX*het$nbY)
    #
    for(x in 1:het$nbX){
        for(y in 1:het$nbY){
            for(yprime in 1:het$nbY){
                if(xFirst){
                    H[x+het$nbX*(y-1),x+het$nbX*(yprime-1)] = oneovermux0[x] + ifelse(y==yprime,oneovermuxy[x,y],0)
                }else{
                    H[(x-1)*het$nbY+y,(x-1)*het$nbY+yprime] = oneovermux0[x] + ifelse(y==yprime,oneovermuxy[x,y],0)
                }
            }
        }
    }
    #
    return(het$sigma*H)
}

dtheta_NablaGstar.logit <- function(het, mu, n, dtheta=diag(1), xFirst=TRUE)
{
    mux0 <- logmu <- ret <- 0
    #
    if(het$outsideOption){
        if(length(dtheta)==0){
            return(matrix(0,nrow=het$nbX*het$nbY,ncol=0))
        }
        #
        mux0 = n - apply(mu,1,sum)
        #
        if(xFirst){
            logmuovermux0 = c(log(mu/mux0))
        }else{
            logmuovermux0 = c(t(log(mu/mux0)))
        }
        #
        ret = matrix(c(dtheta)*logmuovermux0,ncol=1)
    }else{
        if(length(dtheta)==0){
            return(matrix(0,nrow=het$nbX*het$nbY,ncol=0))
        }
        #
        if(xFirst){
            logmu = c(log(mu))
        }else{
            logmu = c(t(log(mu)))
        }
        #
        ret = matrix(c(dtheta)*logmu,ncol=1)
    }
    #
    return(ret)
}

Gbarx.logit <- function(het, Ubarx, mubarx, x)
{
    if(het$outsideOption){
        TOL   = 1e-100
        #
        sigma = het$sigma
        expUbarx = exp(Ubarx/sigma)
        #
        differMargx <- function(z){z + sum(pmin(z*expUbarx,mubarx)) - 1}
        #
        mux0 = uniroot(differMargx,c(0,1),tol=TOL)$root      
        mux = pmin(mux0*expUbarx,mubarx)
        #
        ret = list(valx = sum(mux*Ubarx) - sigma*(mux0*log(mux0)+sum(mux*log(mux))),
                   Ux   = sigma*log(mux/mux0),
                   mux  = mux)
        return(ret)
    }else{
        stop("Method Gbarx not implemented yet when outsideOption==F")
    }
}

simul.logit <- function(het, nbDraws, seed=NULL)
{
    set.seed(seed)
    #
    epsilon_biy <- 0
    if(het$outsideOption){
        epsilon_biy = array(digamma(1) - het$sigma*log(-log(runif(nbDraws*het$nbX*(het$nbY+1)))), dim=c(nbDraws,het$nbY+1,het$nbX))
    }else{
        epsilon_biy = array(digamma(1) - het$sigma*log(-log(runif(nbDraws*het$nbX*het$nbY))), dim=c(nbDraws,het$nbY,het$nbX))
    }
    #
    ret = list(nbX=het$nbX, nbY=het$nbY,
               nbParams=length(epsilon_biy),
               atoms=epsilon_biy,
               aux_nbDraws=nbDraws, 
               xHomogenous=FALSE,
               outsideOption=het$outsideOption)
    class(ret) = "empirical"
    #
    return(ret)
}
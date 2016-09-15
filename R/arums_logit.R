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

build_logit <- function(nbX, nbY, sigma=1, outsideOption=TRUE)
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

G.logit <- function(arums, U, n)
{
    expU = exp(U/arums$sigma)
    #
    if(arums$outsideOption){
        denom = 1 + apply(expU,1,sum)
    }else{
        denom = apply(expU,1,sum)
    }
    #
    ret = list(val = arums$sigma*sum(n*log(denom)),
               mu  = (n/denom)*expU)
    #
    return(ret)
}

Gstar.logit <- function(arums, mu, n)
{
    mux0 <- ret <- 0
    #
    if(arums$outsideOption){
        mux0 = n - apply(mu,1,sum)  
        ret = list(val = arums$sigma*(sum(mu*log(mu/n)) + sum(mux0*log(mux0/n))),
                   U   = arums$sigma*log(mu/mux0))
    }else{
        ret = list(val = arums$sigma*(sum(mu*log(mu/n))),
                   U   = arums$sigma*log(mu/n))
    }
    #
    return(ret)
}

Gstarx.logit <- function(arums, mux, x)
{   # Keith: where is input 'x' used?
    mu0 <- ret <- 0
    #
    if(arums$outsideOption){
        mu0 = 1 - sum(mux)
        ret = list(valx = arums$sigma*(mu0*log(mu0) + sum(mux*log(mux))),
                   Ux   = arums$sigma*log(mux/mu0))
    }else{
        ret = list(valx = arums$sigma*(sum(mux*log(mux))),
                   Ux   = arums$sigma*log(mux))
    }
    #
    return(ret)
}

D2G.logit <- function(arums, U, n, xFirst=TRUE)
{
    # NOTE: the formula is the same regardless of whether outsideOption == TRUE or FALSE
    muxy = G(arums,U,n)$mu
    H = matrix(0,arums$nbX*arums$nbY,arums$nbX*arums$nbY)
    #
    for(x in 1:arums$nbX){
        for(y in 1:arums$nbY){
            for(yprime in 1:arums$nbY){
                if(xFirst){
                    H[x+arums$nbX*(y-1),x+arums$nbX*(yprime-1)] = ifelse(y==yprime,
                                                                     muxy[x,y]*(1-muxy[x,y]/n[x])/arums$sigma,
                                                                     -muxy[x,y]*muxy[x,yprime]/(n[x]*arums$sigma))
                }else{
                    H[(x-1)*arums$nbY+y,(x-1)*arums$nbY+yprime] = ifelse(y==yprime,
                                                                     muxy[x,y]*(1-muxy[x,y]/n[x])/arums$sigma,
                                                                     -muxy[x,y]*muxy[x,yprime]/(n[x]*arums$sigma))
                }
            }
        }
    }
    #
    return(H)
}

D2Gstar.logit <- function(arums, mu, n, xFirst=TRUE)
{ 
    mux0 = n - apply(mu,1,sum)
    #
    oneovermux0 = 1/mux0 
    oneovermuxy = 1/mu
    #
    #H = Diagonal(arums$nbX*arums$nbY,0)
    H = matrix(0,arums$nbX*arums$nbY,arums$nbX*arums$nbY)
    #
    for(x in 1:arums$nbX){
        for(y in 1:arums$nbY){
            for(yprime in 1:arums$nbY){
                if(xFirst){
                    H[x+arums$nbX*(y-1),x+arums$nbX*(yprime-1)] = oneovermux0[x] + ifelse(y==yprime,oneovermuxy[x,y],0)
                }else{
                    H[(x-1)*arums$nbY+y,(x-1)*arums$nbY+yprime] = oneovermux0[x] + ifelse(y==yprime,oneovermuxy[x,y],0)
                }
            }
        }
    }
    #
    return(arums$sigma*H)
}

dtheta_NablaGstar.logit <- function(arums, mu, n, dtheta=diag(1), xFirst=TRUE)
{
    mux0 <- logmu <- ret <- 0
    #
    if(arums$outsideOption){
        if(length(dtheta)==0){
            return(matrix(0,nrow=arums$nbX*arums$nbY,ncol=0))
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
            return(matrix(0,nrow=arums$nbX*arums$nbY,ncol=0))
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

Gbarx.logit <- function(arums, Ubarx, mubarx, x)
{
    #
    if(arums$outsideOption){
        TOL   = 1e-100 # Keith: fix later
        #
        sigma = arums$sigma
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

simul.logit <- function(arums, nbDraws, seed=NULL)
{
    set.seed(seed)
    #
    epsilon_biy <- 0
    if(arums$outsideOption){
        epsilon_biy = array(digamma(1) - arums$sigma*log(-log(runif(nbDraws*arums$nbX*(arums$nbY+1)))), dim=c(nbDraws,arums$nbY+1,arums$nbX))
    }else{
        epsilon_biy = array(digamma(1) - arums$sigma*log(-log(runif(nbDraws*arums$nbX*arums$nbY))), dim=c(nbDraws,arums$nbY,arums$nbX))
    }
    #
    ret = list(nbX=arums$nbX, nbY=arums$nbY,
               nbParams=length(epsilon_biy),
               atoms=epsilon_biy,
               aux_nbDraws=nbDraws, 
               xHomogenous=FALSE,
               outsideOption=arums$outsideOption)
    class(ret) = "empirical"
    #
    return(ret)
}
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
#
# A. Galichon, B. Salanie: " Cupid's Invisible Hand: Social Surplus and Identification in Matching Models"
# A. Galichon, Yu-Wei Hsieh: "Love and Chance: Equilibrium and Identification in a Large NTU matching markets with stochastic choice"
# K. Chiong, A. Galichon, M. Shum: "Duality in dynamic discrete choice models"
#           

G.default <- function(arums, U, n)
{
    val = 0
    mu = matrix(0,arums$nbX,arums$nbY)
    #
    for(x in 1:arums$nbX){
        resGx = Gx(arums,U[x,],x)
        #
        val = val + n[x]*resGx$valx
        mu[x,] = n[x]*resGx$mux
    }
    #
    return(list(val=val,mu=mu))
}

Gstar.default <- function(arums, mu, n)
{
    val = 0
    Uopt = matrix(0,arums$nbX,arums$nbY)
    #
    for(x in 1:arums$nbX){
        resx = Gstarx(arums,mu[x,]/n[x],x)
        #
        val = val + n[x]*resx$valx
        Uopt[x,] = resx$Ux
    }
    #
    return(list(val=val,U=Uopt))
}

Gstarx.default <- function(arums, mux, x)
{ 
    if(!arums$outsideOption){
        stop("Gstarx.default not supported for outsideOption==F")
    }
    #
    thef <- function(Ux)
    {
        res = Gx(arums,Ux,x)
        return(list("objective" =  res$valx - sum(mux*Ux), "gradient" = res$mux - mux))
    }
    #
    resopt = nloptr(x0= rep(0,arums$nbY), eval_f = thef,
                    opt = list(algorithm = 'NLOPT_LD_LBFGS', xtol_rel = xtol_rel, "ftol_rel"=1e-15))
    #
    return(list(valx = -resopt$objective, Ux = resopt$solution))
}

Gbar.default <- function(arums, Ubar, n, mubar)
{
    if(!arums$outsideOption){
        stop("Method Gbar.default not implemented yet when outsideOption==F")
    }
    #
    val = 0
    Uopt = matrix(0,arums$nbX,arums$nbY)
    muopt = matrix(0,arums$nbX,arums$nbY)
    #
    for(x in 1:arums$nbX){
        resx = Gbarx(arums, Ubar[x,],mubar[x,]/n[x], x)
        val = val + n[x]*resx$valx
        Uopt[x,] = resx$Ux
        muopt[x,] = n[x]*resx$mux
    }
    #
    return(list(val=val,U=Uopt,mu=muopt))
}

Gbarx.default <- function(arums, Ubarx, mubarx, x)
{ 
    if(!arums$outsideOption){
        stop("Method Gbarx.default not implemented yet when outsideOption==F")
    }
    #
    thef <- function(mux){
        res = Gstarx(arums,mux,x)
        #
        return(list("objective" = res$valx - sum(mux*Ubarx),
                    "gradient" = res$Ux - Ubarx))
    }
    #
    theg <- function(mux){
        return(list("constraints" = sum(mux) - 1, 
                    "jacobian" = rep(1,arums$nbY)))
    }
    #
    lb = rep(0,arums$nbY)
    ub = mubarx
    #
    resopt = nloptr(x0 = mubarx/2, eval_f = thef, eval_g_ineq = theg,
                    lb=lb, ub=ub,
                    opt = list("algorithm" = "NLOPT_LD_AUGLAG", 
                               local_opts = list("algorithm" = "NLOPT_LD_MMA","xtol_rel"=1e-7),
                               "xtol_rel"=1e-7))
    #
    return(list(valx = -resopt$objective, Ux = Gstarx(arums,resopt$solution,x)$Ux, mux = resopt$solution))
}

D2Gstarx.default <- function (arums, mux, x)
{
    nablaGstar <- function(themux) (Gstarx(arums,themux,x)$Ux)
    return(jacobian(nablaGstar,mux))
}

D2Gstar.default <- function (arums, mu,n, xFirst=T)
{
    if(!arums$outsideOption){
        stop("Method D2Gstar.default not implemented yet when outsideOption==F")
    }
    #
    Hess = Diagonal(arums$nbX*arums$nbY,0)
    #
    for(x in 1:arums$nbX){
        if(xFirst){
            inds=x+arums$nbX*(0:(arums$nbY-1))
            Hess[inds,inds] = D2Gstarx(arums,mu[x,]/n[x],x) / n[x]
        }else{
            inds= (1:arums$nbY) + (x-1)*arums$nbY
            Hess[inds,inds] = D2Gstarx(arums,mu[x,]/n[x],x) / n[x]
        }
    }
    #
    return(Hess)
}

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

build_market_NTU_logit <- function(n, m, alpha, gamma, sigma=1)
{
    nbX = length(n)
    nbY = length(m)
    #
    NTUs = build_NTUs(alpha,gamma)
    logitM = build_logits(nbX,nbY,sigma)
    logitW = build_logits(nbY,nbX,sigma)
    #
    ret = list(n=n,m=m,
               hetG=logitM,hetH=logitW,
               transfers=NTUs)
    class(ret) = "NTU_logit"
    #
    return(ret)
}

solveEquilibrium.NTU_logit = ipfp

# margxInv.NTU_logit = function(xs,mkt,Mu0ys,sigma=1)
# {
#   if (is.null(xs)) {xs = 1:mkt$transfers$nbX}
#   if (!is.null(mkt$neededNorm)) {stop('not supported yet')}
#   themux0s = inversePWA(mkt$n[xs],t ( t(exp(mkt$transfers$gamma[xs,]/sigma - mkt$transfers$alpha[xs,]/sigma) )* Mu0ys ),mkt$transfers$aux_expalpha[xs,]^(1/sigma)) 
#   return(themux0s)
# }

# margyInv.NTU_logit = function(ys,mkt,Mux0s,sigma=1)
# {
#   if (is.null(ys)) {ys = 1:mkt$transfers$nbY}
#   if (!is.null(mkt$neededNorm)) {stop('not supported yet')}
#   themu0ys = inversePWA(mkt$m[ys],  t(exp( mkt$transfers$alpha[,ys]/sigma -mkt$transfers$gamma[,ys]/sigma) * Mux0s) ,t(mkt$transfers$aux_expgamma[,ys]^(1/sigma)) ) 
#   return(themu0ys)
# }

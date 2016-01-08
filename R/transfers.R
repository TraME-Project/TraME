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
# I. Singer: Abstract Convex Analysis. Wiley.
# L. Samuelson, and G. Noldeke: "Implementation Duality".
# A. Galichon, S.D. Kominers, and S. Weber: "An Empirical Framework for Matching with Imperfectly Transferable Utility". 
# O. Bonnet, A. Galichon, and M. Shum: "Yoghurt Chooses Man: The Matching Approach to Identification of Nonadditive Random Utility Models".
#

Psi <- function(tr, ...) UseMethod("Psi")

du_Psi <- function(tr, ...) UseMethod("du_Psi")

dtheta_Psi <- function(tr, ...) UseMethod("dtheta_Psi")

determineType <- function(tr, ...) UseMethod("determineType")

transfersTranspose <- function(tr, ...) UseMethod("transfersTranspose")

Ucal <- function(tr, ...) UseMethod("Ucal") 

Vcal <- function(tr, ...) UseMethod("Vcal")

UW <- function(tr, ...) UseMethod("UW")

dw_UW = function(tr, ...) UseMethod("dw_UW")

VW <- function(tr, ...) UseMethod("VW")

dw_VW = function(tr, ...) UseMethod("dw_VW")

WU <- function(tr, ...) UseMethod("WU")

WV <- function(tr, ...)  UseMethod("WV")

MMF <- function(tr, ...) UseMethod("MMF")

MMF.default <- function(tr, mux0s, mu0ys, xs, ys, sigma=1) # Keith: I broke this expression up into pieces
{
    term_2 = matrix(-sigma*log(mux0s),nrow=tr$nbX,ncol=tr$nbY)
    term_3 = matrix(-sigma*log(mu0ys),nrow=tr$nbX,ncol=tr$nbY,byrow=TRUE)
    #
    term_exp = Psi(tr,term_2,term_3)[xs,ys]
    #
    ret = exp(-term_exp/sigma)
    #
    return(ret)
} 

ufromvs <- function(tr,...) UseMethod("ufromvs")

ufromvs.default <- function(tr, v, tol=0)
{  
    us = Ucal(tr,v,1:tr$nbX,1:tr$nbY)
    u = pmax(apply(us,1,max),0)
    #
    subdiff = matrix(0,tr$nbX,tr$nbY)
    subdiff[which(abs(u-us) <= tol)] = 1
    #
    return(list(u=u,subdiff=subdiff))
}

vfromus <- function(tr,...) UseMethod("vfromus")

vfromus.default <- function(tr, u, tol=0)
{
    vs = Vcal(tr,u,1:tr$nbX,1:tr$nbY)
    v = pmax(apply(vs,2,max),0)
    #
    tsubdiff = matrix(0,tr$nbY,tr$nbX)
    tsubdiff[which(abs(v-t(vs)) <= tol)] = 1
    #
    return(list(v=v,subdiff=t(tsubdiff) ))
}

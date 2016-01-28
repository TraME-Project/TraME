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

# heterogeneity-related classes
G <- function(het, ...) UseMethod("G")

Gx <- function (het, ...) UseMethod("Gx")

Gstar <- function (het, ...) UseMethod("Gstar")

Gstarx <- function (het, ...) UseMethod("Gstarx")

D2Gx <- function (het, ...) UseMethod("D2Gx")

D2Gstarx <- function (het, ...) UseMethod("D2Gstarx")

D2G <- function (het, ...) UseMethod("D2G")

D2Gstar <- function (het, ...) UseMethod("D2Gstar")

dtheta_NablaGstar <- function (het, ...) UseMethod("dtheta_NablaGstar")

Gbar <- function(het, ...) UseMethod("Gbar")

Gbarx <- function(het, ...) UseMethod("Gbarx")

simul <- function (heterog, ...) UseMethod("simul")

# market classes

margxInv <- function(xs, mkt, ...) UseMethod("margxInv",mkt)

margyInv <- function(ys, mkt, ...) UseMethod("margyInv",mkt)

solveEquilibrium <- function(market, ...) UseMethod("solveEquilibrium")

# models

parametricMarket <- function(model, ...) UseMethod("parametricMarket")

initparam <- function(model, ...) UseMethod("initparam")

dparam <- function(model, ...) UseMethod("dparam")

estimate <- function(model, ...) UseMethod("estimate")

# transfer classes

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

ufromvs <- function(tr,...) UseMethod("ufromvs")

vfromus <- function(tr,...) UseMethod("vfromus")

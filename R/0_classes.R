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

# arums-related classes
G <- function(arums, ...) UseMethod("G")

Gx <- function (arums, ...) UseMethod("Gx")

Gstar <- function (arums, ...) UseMethod("Gstar")

Gstarx <- function (arums, ...) UseMethod("Gstarx")

D2Gx <- function (arums, ...) UseMethod("D2Gx")

D2Gstarx <- function (arums, ...) UseMethod("D2Gstarx")

D2G <- function (arums, ...) UseMethod("D2G")

D2Gstar <- function (arums, ...) UseMethod("D2Gstar")

dtheta_NablaGstar <- function (arums, ...) UseMethod("dtheta_NablaGstar")

Gbar <- function(arums, ...) UseMethod("Gbar")

Gbarx <- function(arums, ...) UseMethod("Gbarx")

simul <- function (arums, ...) UseMethod("simul")

# mmfs class

mmfsTranspose <- function(mmfs) UseMethod("mmfsTranspose")

M <- function(mmfs,...) UseMethod("M")

margxInv <- function(xs, mmfs, ...) UseMethod("margxInv",mmfs)

margyInv <- function(ys, mmfs, ...) UseMethod("margyInv",mmfs)

# market classes

solveEquilibrium <- function(market, ...) UseMethod("solveEquilibrium")

# models

parametricMarket <- function(model, ...) UseMethod("parametricMarket")

initparam <- function(model, ...) UseMethod("initparam")

dparam <- function(model, ...) UseMethod("dparam")

estimate <- function(model, ...) UseMethod("estimate")

# transfer classes

Psi <- function(tr, ...) UseMethod("Psi")

Psi_sub <- function(tr, ...) UseMethod("Psi_sub")

du_Psi <- function(tr, ...) UseMethod("du_Psi")

du_Psi_sub <- function(tr, ...) UseMethod("du_Psi")

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

M <- function(mmfs, ...) UseMethod("M")


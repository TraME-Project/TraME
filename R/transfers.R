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
# I. Singer: Abstract Convex Analysis. Wiley.
# L. Samuelson, and G. Noldeke: "Implementation Duality".
# A. Galichon, S.D. Kominers, and S. Weber: "An Empirical Framework for Matching with Imperfectly Transferable Utility". 
# O. Bonnet, A. Galichon, and M. Shum: "Yoghurt Chooses Man: The Matching Approach to Identification of Nonadditive Random Utility Models".
#
UW<- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
  return(-Psi_sub(tr,0,-Ws,xs,ys))
}
#
VW <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
  return(-Psi_sub(tr,Ws,0,xs,ys))
}
#
dw_UW <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
  return(1-du_Psi_sub.TU(tr,0,-Ws,xs,ys) )
}
#
dw_VW <- function(tr, Ws, xs=1:tr$nbX, ys=1:tr$nbY)
{
  return(-du_Psi_sub.TU(tr,Ws,0,xs,ys) )
}


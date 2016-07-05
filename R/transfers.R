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
################################################################################
########################    Default and generic methods ########################
################################################################################
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
################################################################################
########################            TU transfers         #######################
################################################################################
build_TUs <- function(phi)
{
  nbX=dim(phi)[1]
  nbY=dim(phi)[2]
  #
  ret = list(nbX = nbX, nbY = nbY,
             nbParams = nbX*nbY,
             phi = phi)
  class(ret) = "TU"
  return(ret)
}
#
transfersTranspose.TU <- function(tr)
{
  ret = list(nbX = tr$nbY, nbY = tr$nbX,
             nbParams = tr$nbParams,
             phi = t(tr$phi)
  )
  class(ret) = "TU"
  return(ret)
}
#
Psi.TU <- function(tr, U, V) ((U + V - tr$phi)/2)
#
Psi_sub.TU <- function(tr,U,V,xs,ys) ((U + V - tr$phi[xs,ys])/2) # here for vectorization purposes only 
#
du_Psi.TU <- function(tr, U, V) (matrix(1/2, nrow=tr$nbX, ncol=tr$nbY))
#
du_Psi_sub.TU <- function(tr, U, V, xs, ys) (matrix(1/2, nrow=length(xs), ncol=length(ys)))
#
dtheta_Psi.TU <- function(tr, U, V, dtheta=NULL) 
{
  ret <- 0
  if(is.null(dtheta)){
    ret = Diagonal(tr$nbX*tr$nbY,-1/2)
  }else{
    ret = -dtheta/2
  }
  return(ret)
}
#
determineType.TU <- function(tr, ...) (1)
#
Ucal.TU <- function(tr, vs, xs=1:tr$nbX, ys=1:tr$nbY) ( tr$phi[xs,ys] - matrix(vs,length(xs),length(ys),byrow=TRUE))
#
Vcal.TU <- function(tr, us, xs=1:tr$nbX, ys=1:tr$nbY) (tr$phi[xs,ys] - us)
#
WU.TU <- function(tr, Us, xs=1:tr$nbX, ys=1:tr$nbY) (2*Us - tr$phi[xs,ys])
#
WV.TU <- function(tr, Vs, xs=1:tr$nbX, ys=1:tr$nbY) (tr$phi[xs,ys] - 2*Vs)
#
################################################################################
########################         End of TU transfers         ###################
################################################################################
#
################################################################################
########################            NTU transfers         ######################
################################################################################
build_NTUs <- function(alpha, gamma) 
{
  nbX = dim(alpha)[1]
  nbY = dim(gamma)[2]
  #
  ret = list(nbX = nbX, nbY = nbY,
             nbParams = 2 * nbX * nbY,
             alpha = alpha, gamma = gamma)
  class(ret) = "NTU"
  return(ret)
}
#
transfersTranspose.NTU <- function(tr)
{
  ret = list(nbX = tr$nbY, nbY = tr$nbX,
             nbParams = tr$nbParams,
             alpha = t(tr$gamma), gamma = t(tr$alpha))
  class(ret) = "NTU"
  return(ret)
}
#
Psi.NTU <- function(tr, U, V) (pmax(U - tr$alpha, V - tr$gamma))
#
Psi_sub.NTU <- function(tr, U, V, xs, ys) (pmax(U - tr$alpha[xs,ys], V - tr$gamma[xs,ys]))
#
du_Psi.NTU <- function(tr, U, V) ( ifelse(U-tr$alpha >= V - tr$gamma,1,0) ) #keith: should this be a vector of ones or zeros?
#
du_Psi_sub.NTU <- function(tr, U, V, xs, ys) (ifelse(U-tr$alpha[xs,ys] >= V - tr$gamma[xs,ys],1,0))
#
dtheta_Psi.NTU <- function(tr, U, V, dtheta=NULL) 
{
  dupsi = c(du_Psi(tr,U,V))
  if(is.null(dtheta)){
    ret = -cbind(Diagonal(x=dupsi),
                 Diagonal(x=1- dupsi)) 
    return(ret)
  }else{
    dtheta1 = dtheta[1:(tr$nbX*tr$nbY)]
    dtheta2 = dtheta[(1+tr$nbX*tr$nbY):(2*tr$nbX*tr$nbY)]
    #
    ret = -c(dupsi*dtheta1 + (1-dupsi)*dtheta2)
    return(ret)
  }
}
#
determineType.NTU <- function(tr, ...) (2)
#
WU.NTU <- function(tr, Us, xs=1:tr$nbX, ys=1:tr$nbY) (Us - tr$alpha[xs,ys])
#
WV.NTU <- function(tr, Vs, xs=1:tr$nbX, ys=1:tr$nbY) (tr$gamma[xs,ys] - Vs)
################################################################################
########################        End of NTU transfers         ###################
################################################################################
#
################################################################################
########################            LTU transfers         ######################
################################################################################
build_LTUs <- function(lambda, phi) 
{
  nbX = dim(lambda)[1]
  nbY = dim(lambda)[2]
  #
  aux_zeta = 1 - lambda
  if(min(c(lambda,aux_zeta)) <= 0){
    stop ("lambda not strictly between 0 and 1")
  }
  #
  ret = list(nbX = nbX, nbY = nbY,
             nbParams = 2* nbX * nbY,
             lambda = lambda, phi = phi, 
             aux_zeta = aux_zeta)
  class(ret) = "LTU"
  return(ret)
}
#
transfersTranspose.LTU <- function(tr)
{
  ret = list(nbX = tr$nbY, nbY = tr$nbX,
             nbParams = tr$nbParams,
             lambda = t(tr$aux_zeta), phi = t(tr$phi),
             aux_zeta = t(tr$lambda))
  class(ret) = "LTU"
  return(ret)
}
#
Psi.LTU <- function(tr, U, V) (tr$lambda * U + tr$aux_zeta * V - tr$phi)
#
Psi_sub.LTU <- function(tr, U, V,xs,ys) (tr$lambda[xs,ys] * U 
                                         + tr$aux_zeta[xs,ys] * V - tr$phi[xs,ys])
#
du_Psi.LTU <- function(tr, ...) ( tr$lambda )
#
du_Psi_sub.LTU <- function(tr, U, V,xs,ys) ( tr$lambda[xs,ys] )
#
dtheta_Psi.LTU <- function(tr, U, V, dtheta=NULL) 
{
  UminusV = c(U - V)
  dtheta1 <- dtheta2 <- ret <- 0
  #
  if(is.null(dtheta)){
    ret = cbind(Diagonal(x=UminusV), Diagonal(tr$nbX*tr$nbY,-1))
    return(ret)
  }else{
    dtheta1 = dtheta[1:(tr$nbX*tr$nbY)]
    dtheta2 = dtheta[(1+tr$nbX*tr$nbY):(2*tr$nbX*tr$nbY)]
    #
    ret = c(UminusV*dtheta1 - dtheta2)
    return(ret)
  }
}
#
determineType.LTU <- function(tr, ...) (1)
#
Ucal.LTU <- function(tr, vs, xs=1:tr$nbX, ys=1:tr$nbY)((tr$phi[xs,ys] - tr$aux_zeta[xs,ys] 
                                                        * matrix(vs,length(xs),length(ys),byrow=TRUE)) 
                                                       / tr$lambda[xs,ys])
#
Vcal.LTU <- function(tr, us, xs=1:tr$nbX, ys=1:tr$nbY) ((tr$phi[xs,ys] 
                                                         - tr$lambda[xs,ys] 
                                                         * matrix(us,length(xs),length(ys))) 
                                                        / tr$aux_zeta[xs,ys])
#
WU.LTU <- function(tr, Us, xs=1:tr$nbX, ys=1:tr$nbY) ((Us - tr$phi[xs,ys]) 
                                                      /  tr$aux_zeta[xs,ys])
#
WV.LTU <- function(tr, Vs, xs=1:tr$nbX, ys=1:tr$nbY) ( (tr$phi[xs,ys] - Vs) 
                                                       /  tr$lambda[xs,ys]  )
################################################################################
########################        End of LTU transfers         ###################
################################################################################
#
################################################################################
########################            ETU transfers         ######################
################################################################################
build_ETUs <- function(alpha, gamma, tau)
{
  nbX = dim(alpha)[1]
  nbY = dim(alpha)[2]
  #
  ret = list(nbX = nbX, nbY = nbY,
             nbParams = 3*nbX*nbY,
             alpha=alpha, gamma=gamma, tau=tau, 
             aux_expminusalphaovertau = exp(-alpha/tau),
             aux_expminusgammaovertau = exp(-gamma/tau))
  class(ret) = "ETU"
  return(ret)
}
#
transfersTranspose.ETU <- function(tr)
{
  ret = list(nbX = tr$nbY, nbY = tr$nbX,
             nbParams = tr$nbParams,
             alpha=t(tr$gamma), gamma=t(tr$alpha), tau=t(tr$tau),
             aux_expminusalphaovertau = t(tr$aux_expminusgammaovertau),
             aux_expminusgammaovertau = t(tr$aux_expminusalphaovertau))
  class(ret) = "ETU"
  return(ret)
}
#
Psi.ETU <- function(tr, U, V) (tr$tau * log((exp(U/tr$tau) *tr$aux_expminusalphaovertau 
                                             + exp(V/tr$tau)*tr$aux_expminusgammaovertau
                                             )/2))
#
Psi_sub.ETU <- function(tr, U, V,xs,ys)
{
  tauxsys = tr$tau[xs,ys]
  term_1 = exp(U/ tauxsys)*tr$aux_expminusalphaovertau[xs,ys]
  term_2 = exp(V/ tauxsys)*tr$aux_expminusgammaovertau[xs,ys]
  ret = tauxsys * log((term_1 + term_2)/2)
  return(ret)
}
#
du_Psi.ETU <- function(tr, U, V) ( 1/(1 + exp((V - U + tr$alpha - tr$gamma)/(tr$tau))) )
#
du_Psi_sub.ETU <- function(tr, U, V, xs, ys) (1/(1 + exp((V - U + tr$alpha[xs,ys] - tr$gamma[xs,ys])
                                                         /(tr$tau[xs,ys]))))
#
dtheta_Psi.ETU <- function(tr, U, V, dtheta=NULL) 
{
  dupsimat = du_Psi(tr,U,V)
  dupsi = c(dupsimat)
  #
  term_1 <- term_2 <- ret <- 0
  #
  if(is.null(dtheta)){
    term_1 = (U - tr$alpha )*dupsi # keith: should this be dupsi_mat?
    term_2 = (V - tr$gamma)*(1-dupsi)
    #
    dsigmapsi = c((Psi(tr,U,V) - term_1 - term_2) / tr$tau)
    #
    ret = cbind(Diagonal(x = -dupsi),
                Diagonal(x = dupsi-1),
                Diagonal(x = dsigmapsi))
    #
    return(ret)
  }else{
    dtheta1 = dtheta[1:(tr$nbX*tr$nbY),]
    dtheta2 = dtheta[(1+tr$nbX*tr$nbY):(2*tr$nbX*tr$nbY),]
    dtheta3 = dtheta[(1+2*tr$nbX*tr$nbY):(3*tr$nbX*tr$nbY),]
    #
    if(min(dtheta3==0)){
      dsigmapsidtheta = 0
    }else{
      term_1 = (U - tr$alpha )*dupsimat
      term_2 = (V - tr$gamma)*(1-dupsimat)
      #
      dsigmapsidtheta = dtheta3*c((Psi(tr,U,V) - term_1 - term_2) / tr$tau)
    }
    #
    ret = c(-dupsi*dtheta1 - (1-dupsi)*dtheta2 + dsigmapsidtheta)
    #
    return(ret)
  }
}
#
determineType.ETU <- function(tr, ...) (2)
#
Ucal.ETU <- function(tr, vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
  term_1 = matrix(vs,tr$nbX,tr$nbY,byrow=TRUE) - tr$gamma[xs,ys]
  term_log = 2 - exp(term_1/tr$tau[xs,ys])
  ret = tr$alpha[xs,ys] + tr$tau[xs,ys] * log(term_log)
  return(ret)
}
#
Vcal.ETU <- function(tr, us, xs=1:tr$nbX, ys=1:tr$nbY)
{
  term_1 = us - tr$alpha[xs,ys]
  term_log = 2 - exp(term_1/tr$tau[xs,ys])
  ret = tr$gamma[xs,ys] + tr$tau[xs,ys] * log(term_log)
  return(ret)
}
#
WU.ETU <- function(tr, Us, xs=1:tr$nbX, ys=1:tr$nbY)
{
  term_1 = 2*exp( (tr$gamma[xs,ys] - Us)/tr$tau[xs,ys] )
  term_2 = exp( (tr$gamma[xs,ys] - tr$alpha[xs,ys])/tr$tau[xs,ys] )
  term_log = term_1 - term_2
  ret = -tr$tau[xs,ys] * log(term_log)
  return(ret)
}
#
WV.ETU <- function(tr, Vs, xs=1:tr$nbX, ys=1:tr$nbY)
{
  term_1 = 2*exp( (tr$alpha[xs,ys] - Vs)/tr$tau[xs,ys] )
  term_2 = exp( (tr$alpha[xs,ys] - tr$gamma[xs,ys])/tr$tau[xs,ys] )
  term_log = term_1 - term_2
  #
  ret = tr$tau[xs,ys] * log(term_log)
}
################################################################################
########################        End of ETU transfers         ###################
################################################################################


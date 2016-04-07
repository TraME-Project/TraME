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

inversePWA_old <- function(a, B, C, k = 1)
{
    nbX = length(a)
    nbY = dim(B)[2]
    #
    vals = rep(0,nbX)
    for(x in 1:nbX){    
        sortB = sort(B[x,],index.return=T)
        sigma = sortB$ix
        
        b = sortB$x
        small_C = C[x,sigma]
        
        ylow = 1
        yup = nbY
        while(yup > ylow){
            ymid = ylow + (yup - ylow) %/% 2
            lhs = k * b[ymid] + sum(small_C * pmin(b[ymid],b))
            if(lhs == a[x]){
                yup = ylow = ymid
            }else if(lhs > a[x]){
                yup = ymid
            }else{
                ylow= ymid + 1
            }
        }
        if((ylow==1) & ( k*b[ylow]+sum(small_C * pmin(b[ylow],b)) >= a[x])){
            vals[x] = a[x] / (k+sum(small_C))
        }else{
            ysincluded = which((1:nbY) <= ylow)
            vals[x]= (a[x] - sum((small_C*b)[ysincluded])) / (k + sum(small_C) - sum(small_C[ysincluded]))
        }
    }
    #
    return(vals)
}

inversePWA <- function(a, B, C, k=1.0)
{
    #
    vals <- .Call("invPWA_R", a,B,C,k, PACKAGE = "TraME")$vals
    #
    return(c(vals))
}

tests_TraME <- function(nbDraws = 1e3)
{
    ptm = proc.time()
    #
    hash_arum <- tests_arum(notifications=FALSE,nbDraws=10*nbDraws)
    hash_equilibrium <- tests_equilibrium(notifications=FALSE,nbDraws=nbDraws)
    hash_estimation <- tests_estimation(notifications=FALSE)
    #
    time = proc.time() - ptm
    message(paste0('All tests completed. Overall time elapsed = ', round(time["elapsed"],5), 's.'))
    #
    hash_vals <- c(hash_arum,hash_equilibrium,hash_estimation)
    message(compare_hashvals(hash_vals))
    return(hash_vals)
}

compare_hashvals <- function(hash_vals)
{
  true_hash <- c("659ba9a641a1415cfe0cd36ded9d097c", "70995dea67d9deeed03d717c0c29d405", "1291c1bcbcfda7db348a15b247939bee")
  #
  if(identical(hash_vals,true_hash)){ conclusion = 'Test results are correct!\n'
  }else{
    conclusion = '*** CAUTION *** There is a problem with the results of: '
    if(!identical(hash_vals[1],true_hash[1])){
      conclusion= paste0(conclusion,'arums tests; ')
    }
    if(!identical(hash_vals[2],true_hash[2])){
      conclusion = paste0(conclusion,'equilibrium tests; ')
    }
    if(!identical(hash_vals[3],true_hash[3])){
      conclusion = paste0(conclusion,'estimation tests; ')
    }
    conclusion = paste0(conclusion,'please check.\n')
  }
  return(conclusion)
}

verify_signature <- function()
{
    output_hide <- capture.output(hash_vals <- suppressMessages(tests_TraME()))
    #
    message(compare_hashvals(hash_vals))

}

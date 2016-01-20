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

inversePWA <- function(a, B, C)
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
            lhs = b[ymid] + sum(small_C * pmin(b[ymid],b))
            if(lhs == a[x]){
                yup = ylow = ymid
            }else if(lhs > a[x]){
                yup = ymid
            }else{
                ylow=ymid+1
            }
        }
        if((ylow==1) & ( b[ylow]+sum(small_C * pmin(b[ylow],b)) >= a[x])){
            vals[x] = a[x] / (1+sum(small_C))
        }else{
            ysincluded = which((1:nbY) <= ylow)
            vals[x]= (a[x] - sum((small_C*b)[ysincluded])) / (1 + sum(small_C) - sum(small_C[ysincluded]))
        }
    }
    #
    return(vals)
}

tests_TraME <- function(){
    ptm = proc.time()
    #
    tests_arum(notifications=FALSE)
    tests_equilibrium(notifications=FALSE)
    tests_estimation(notifications=FALSE)
    #
    time = proc.time() - ptm
    message(paste0('All tests completed. Overall time elapsed = ', round(time["elapsed"],5), 's.'))
}

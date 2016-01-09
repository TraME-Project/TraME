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

init_TraME <- function(nbSlaves = 0, isMaster=TRUE, withGurobi=TRUE)
{
    assign("TraME_withGurobi", withGurobi, envir = .GlobalEnv)
    assign("TraME_nbSlaves", nbSlaves, envir = .GlobalEnv)
    assign("TraME_isMaster", isMaster, envir = .GlobalEnv)
    #
    if(withGurobi){
        require('gurobi')
    }else{
        warning("Initialization without Gurobi. LP-based algorithms will not run.")
    }
    #
    if(isMaster){    # this is the master; should initialize the workers, if any
        if(nbSlaves==0){
            sfInit(parallel=FALSE)
        }else{
            sfInit(parallel=TRUE, cpus=nbSlaves, type="SOCK")
            sfExportAll()
            sfSapply(x=1:nbSlaves,fun=function(myindex){init_TraME(nbSlaves = 0,isMaster = F, withGurobi = withGurobi)})
        }
    } 
}

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
        c = C[x,sigma] # Keith: change this!!!
        
        ylow = 1
        yup = nbY
        while(yup > ylow){
            ymid = ylow + (yup - ylow) %/% 2
            lhs = b[ymid] + sum(c * pmin(b[ymid],b))
            if(lhs == a[x]){
                yup = ylow = ymid
            }else if(lhs > a[x]){
                yup = ymid
            }else{
                ylow=ymid+1
            }
        }
        if((ylow==1) & ( b[ylow]+sum(c * pmin(b[ylow],b)) >= a[x])){
            vals[x] = a[x] / (1+sum(c))
        }else{
            ysincluded = which((1:nbY) <= ylow)
            vals[x]= (a[x] - sum((c*b)[ysincluded])) / (1 + sum(c) - sum(c[ysincluded]))
        }
    }
    #
    return(vals)
}

tests_TraME <- function(withGurobi=TRUE){
    ptm = proc.time()
    tests_arum(notifications=FALSE)
    tests_equilibrium(notifications=FALSE)
    tests_estimation(notifications=FALSE)
    time = proc.time() - ptm
    message(paste0('All tests completed. Overall time elapsed=', time["elapsed"], 's.'))
}

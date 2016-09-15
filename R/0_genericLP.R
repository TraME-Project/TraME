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

genericLP <- function(obj, A, modelsense, rhs, sense, Q=NULL, lb=NULL, ub=NULL, start=NULL, Method = -1){
    #
    if(.trame_lp_options$gurobi_exists==TRUE){
        gurobiModel = list(A=A,obj=obj,modelsense=modelsense,rhs=rhs,sense=sense,
                           Q=Q,lb=lb,ub=ub,start=start)
        result = gurobi(gurobiModel, params=list(OutputFlag=0 , Method = Method))
        #
        if(result$status=="OPTIMAL"){
            solution = result$x
            objval = result$objval
            pi = result$pi
            rc = result$rc
        }else{
            stop("optimization problem with Gurobi")
        }
    }else if(.trame_lp_options$glpk_exists==TRUE){
        if(is.null(Q)==FALSE){
            warning("GLPK cannot solve quadratic programming problems.\n")
        }
        #
        bounds = list()
        if(is.null(lb)==FALSE){
            bounds$lower = list(ind = 1:length(lb), val = lb)
        }
        if(is.null(ub)==FALSE){
            bounds$upper = list(ind = 1:length(ub), val = ub)
        }
        #
        dir_glpk = sense
        if(length(dir_glpk)==1){
            dir_glpk = rep(dir_glpk,nrow(A))
        }
        if(any(dir_glpk=="=")){
            dir_glpk[dir_glpk=="="] = "=="
        }
        #
        max_glpk = ifelse(modelsense=="max",TRUE,FALSE)
        #
        result = Rglpk_solve_LP(obj=obj,mat=A,dir=dir_glpk,rhs=rhs,bounds=bounds,max=max_glpk)
        #
        if(result$status==0){
            solution = result$solution
            objval = result$optimum
            pi = result$pi
            rc = result$rc
        }else{
            stop("optimization problem with GLPK")
        }
    }else{
        stop("No LP solver available.")
    }
    #
    ret = list(solution=solution, objval=objval,
               pi=pi, rc=rc)
    #
    return(ret)
}

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

.onAttach <- function(...)
{
    .trame_lp_options <<- list()
    #
    gurobi_exists <- "gurobi" %in% rownames(installed.packages())
    glpk_exists <- "Rglpk" %in% rownames(installed.packages())
    #
    .trame_lp_options$gurobi_exists <<- gurobi_exists
    .trame_lp_options$glpk_exists <<- glpk_exists
    #
    if(gurobi_exists==TRUE||glpk_exists==TRUE){
        obj <- c(2, 4, 3)
        mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
        dir <- c("<=", "<=", "<=")
        rhs <- c(60, 40, 80)
        #
        if(gurobi_exists==TRUE){ # test a simple problem
            #
            gurobiModel = list(A=mat,obj=obj,modelsense="max",rhs=rhs,sense=dir)
            result_gurobi = try(gurobi::gurobi(gurobiModel, params=list(OutputFlag=0)),silent=TRUE)
            #
            if(class(result_gurobi)=="try-error"){
                .trame_lp_options$gurobi_works <<- FALSE
                packageStartupMessage("Gurobi seems to be present on your machine, but failed to solve a simple example.\n")
            }else if(round(result_gurobi$objval,0)==77){
                .trame_lp_options$gurobi_works <<- TRUE
                packageStartupMessage("TraME will use Gurobi to solve LP problems.")
            }
        }else if(glpk_exists==TRUE){
            #
            result_glpk = try(Rglpk::Rglpk_solve_LP(obj, mat, dir, rhs, max = TRUE),silent=TRUE)
            #
            if(class(result_glpk)=="try-error"){
                .trame_lp_options$glpk_works <<- FALSE
                packageStartupMessage("GLPK seems to be present on your machine, but failed to solve a simple example.\n")
            }else if(round(result_glpk$optimum,0)==77){
                .trame_lp_options$glpk_works <<- TRUE
                packageStartupMessage("TraME will use GLPK to solve LP problems.")
            }
        }
    }else{
        packageStartupMessage("Could not find any LP solver.\n")
    }
}
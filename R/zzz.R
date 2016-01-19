.onAttach <- function(...)
{
    .trame_gurobi_exists <<- "gurobi" %in% rownames(installed.packages())
    .trame_glpk_exists <<- "Rglpk" %in% rownames(installed.packages())
    #
    if(.trame_gurobi_exists==TRUE){ # test a simple problem
        obj <- c(2, 4, 3)
        mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
        dir <- c("<=", "<=", "<=")
        rhs <- c(60, 40, 80)
        max <- TRUE
        #
        gurobiModel = list(A=mat,obj=obj,modelsense="max",rhs=rhs,sense=dir)
        result_gurobi = try(gurobi::gurobi(gurobiModel, params=list(OutputFlag=0)),silent=TRUE)
        #
        if(class(result_gurobi)=="try-error"){
            .trame_gurobi_works <<- FALSE
            message("Gurobi seems to be present on your machine, but failed to solve a simple example.\n")
        }else if(round(result_gurobi$objval,0)==77){
            .trame_gurobi_works <<- TRUE
            message("TraME will use Gurobi to solve LP problems.")
        }
    }else if(.trame_glpk_exists==TRUE){
        obj <- c(2, 4, 3)
        mat <- matrix(c(3, 2, 1, 4, 1, 3, 2, 2, 2), nrow = 3)
        dir <- c("<=", "<=", "<=")
        rhs <- c(60, 40, 80)
        max <- TRUE
        #
        result_glpk = try(result_glpk = Rglpk::Rglpk_solve_LP(obj, mat, dir, rhs, max = max),silent=TRUE)
        #
        if(class(result_gurobi)=="try-error"){
            .trame_glpk_works <<- FALSE
            message("GLPK seems to be present on your machine, but failed to solve a simple example.\n")
        }else if(round(result_gurobi$objval,0)==77){
            .trame_glpk_works <<- TRUE
            message("TraME will use GLPK to solve LP problems.")
        }
    }else{
        message("Could not find any LP solver.\n")
    }
}
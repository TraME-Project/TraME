genericLP <- function(obj,A,modelsense,rhs,sense,lb=NULL,ub=NULL){
    #
    if(.trame_gurobi_exists==TRUE){
        gurobiModel = list(A=A,obj=obj,modelsense=modelsense,rhs=rhs,sense=sense,lb=lb,ub=ub)
        result = gurobi(gurobiModel, params=list(OutputFlag=0))
        #
        if(result$status=="OPTIMAL"){
            solution = result$x
            objval = result$objval
            pi = result$pi
            rc = result$rc
        }else{
            stop("optimization problem with Gurobi")
        }
    }else if(.trame_glpk_exists==TRUE){
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

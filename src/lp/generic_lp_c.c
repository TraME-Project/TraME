/*################################################################################
  ##
  ##   Copyright (C) 2015-2018 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
  ##
  ##   This file is part of TraME.
  ##
  ##   TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * Generic Linear/Quadratic Programming
 *
 * Keith O'Hara
 * 05/08/2016
 *
 * This version:
 * 02/26/2017
 */
 
#include <string.h>

#include "lp/generic_lp_c.h"

#if defined(TRAME_USE_GUROBI)
    #include "trame_gurobi.c"
#elif defined(TRAME_USE_GLPK)
    #include "trame_glpk.c"
#else
    
#endif

//
// Dense setup 
int generic_LP_C(int rows, int cols, double* obj, double* A, int model_opt_sense, 
                 double* rhs, char* constr_sense, double* Q, double* lb, double* ub, 
                 double* objval, double* sol_mat_X, double* sol_mat_RC, 
                 double* dual_mat_PI, double* dual_mat_SLACK)
{
    int success = 0;
#if defined(TRAME_USE_GUROBI)
    success = trame_gurobi_switch(rows,cols,obj,A,model_opt_sense,rhs,constr_sense,Q,lb,ub,objval,sol_mat_X,sol_mat_RC,dual_mat_PI,dual_mat_SLACK);
#elif defined(TRAME_USE_GLPK)
    if (Q) {
        printf("GLPK does not support quadratic programming problems\n");
    }
    success = trame_glpk(rows,cols,obj,A,model_opt_sense,rhs,constr_sense,lb,ub,objval,sol_mat_X,sol_mat_RC,dual_mat_PI,dual_mat_SLACK);
#else
    printf("No linear programming kit detected\n");
    success = 0;
#endif

    return success;
}

//
// for use with sparse matrix inputs; sparse A, dense (or empty) Q
int generic_LP_C_sparse(int rows, int cols, double* obj, int numnz, int* vbeg, int* vind, double* vval, 
                        int model_opt_sense, double* rhs, char* constr_sense, double* Q, double* lb, double* ub, 
                        double* objval, double* sol_mat_X, double* sol_mat_RC, double* dual_mat_PI, double* dual_mat_SLACK)
{
    int success = 0;
#if defined(TRAME_USE_GUROBI)
    success = trame_gurobi_sparse(rows,cols,obj,numnz,vbeg,vind,vval,model_opt_sense,rhs,constr_sense,Q,lb,ub,objval,sol_mat_X,sol_mat_RC,dual_mat_PI,dual_mat_SLACK);
#elif defined(TRAME_USE_GLPK)
    if (Q) {
        printf("GLPK does not support quadratic programming problems\n");
    }
    success = trame_glpk_sparse(rows,cols,obj,numnz,vbeg,vind,vval,model_opt_sense,rhs,constr_sense,lb,ub,objval,sol_mat_X,sol_mat_RC,dual_mat_PI,dual_mat_SLACK);
#else
    printf("No linear programming kit detected\n");
    success = 0;
#endif

    return success;
}

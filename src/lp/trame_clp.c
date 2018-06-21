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

#include <Clp_C_Interface.h>

// clang-mp-6.0 -O3 -I$COIN_UTILS_INCLUDE_PATH -I$COIN_LP_INCLUDE_PATH clp_c2.c -o clp_c2.out -L$COIN_UTILS_LIB_PATH -lCoinUtils -L$COIN_LP_LIB_PATH -lClp

int trame_clp(int n_constr, int n_vars, double* __restrict__ obj, double* __restrict__ A, int model_opt_sense, 
              double* __restrict__ rhs, char* __restrict__ constr_dir, double* lb, double* ub, 
              double* __restrict__ objval, double* __restrict__ sol_mat_X, double* __restrict__ sol_mat_RC, 
              double* __restrict__ dual_mat_PI, double* __restrict__ dual_mat_SLACK)
{
    Clp_Simplex *lp;

    lp = Clp_newModel();

    //

    if (model_opt_sense==0) { // minimize
        Clp_setOptimizationDirection(lp, 1.0); // −1: maximize, 1: minimize
    } else if (model_opt_sense==1) { // maximize
        Clp_setOptimizationDirection(lp, -1.0); // −1: maximize, 1: minimize
    } else {
        printf("unrecognized input for model_opt_sense; should be 0 (minimize) or 1 (maximize).\n");
        goto QUIT;
    }

    //

    int* row_inds = malloc(sizeof(int) * (n_var*n_constr)); // ia 
    int* col_inds = malloc(sizeof(int) * (n_var+1)); // ja

    col_inds[0] = 0;

    for (int i=0; i < n_vars; i++)
    {
        for (int j=0; j < n_constr; j++)
        {
            row_inds[j+i*n_constr] = j;
        }
        col_inds[i+1] = (i+1)*n_constr;
    }

    //

    Clp_loadProblem(lp,n_vars,n_constr,col_inds,row_inds,A,clb,NULL,obj,NULL,rub);

    //

    int ret = Clp_initialSolve(lp);

    //

    if (objval)
    {
        *objval = Clp_objectiveValue(lp);
    }

    // column primal

    sol_mat_X      = Clp_primalColumnSolution(lp);
    dual_mat_Slack = Clp_primalRowSolution(lp);

    sol_mat_RC  = Clp_dualColumnSolution(lp);
    dual_mat_PI = Clp_dualRowSolution(lp);

    //

QUIT:
    free(row_inds);
    free(col_inds);
    Clp_deleteModel(lp);

    return 0;
}

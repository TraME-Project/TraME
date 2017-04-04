/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
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
 * Generic Linear/Quadratic Programming solver using Gurobi
 *
 * Keith O'Hara
 * 05/08/2016
 *
 * This version:
 * 02/22/2017
 */

#include "trame.hpp"

extern "C" {
    #include "lp/generic_lp_c.h"
}

bool 
trame::generic_LP(int k, int n, double *obj, double* A, int model_opt_sense, double* rhs, char* constr_sense, double* Q, double* lb, double* ub, double* start, double& objval, double* sol_mat_X, double* sol_mat_RC, double* dual_mat_PI, double* dual_mat_SLACK)
{
    // k: number of constraints ('rows')
    // n: number of variables ('columns')
    //
    // Initialize
    bool success = false;

    // Call C-version of the solver
    int solved = generic_LP_C(k, n, obj, A, model_opt_sense, rhs, constr_sense, Q, lb, ub, 
                              &objval, sol_mat_X, sol_mat_RC, dual_mat_PI, dual_mat_SLACK);
    //
    //
    if (solved == 1) {
        success = true;
    }
    //
    return success;
}
    
bool 
trame::generic_LP(int k, int n, double *obj, int numnz, int* vbeg, int* vind, double* vval, int model_opt_sense, double* rhs, char* constr_sense, double* Q, double* lb, double* ub, double* start, double& objval, double* sol_mat_X, double* sol_mat_RC, double* dual_mat_PI, double* dual_mat_SLACK)
{
    // k: number of constraints ('rows')
    // n: number of variables ('columns')
    //
    // Initialize
    bool success = false;
    //
    // Call C-version of the solver
    int solved = generic_LP_C_sparse(k, n, obj, numnz, vbeg, vind, vval, model_opt_sense, rhs, constr_sense, Q, lb, ub, 
                                     &objval, sol_mat_X, sol_mat_RC, dual_mat_PI, dual_mat_SLACK);
    //
    // Put solution matrices together
    if (solved == 1) {
        success = true;
    }
    //
    return success;
}

/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
  ##      Simon Weber
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
 * Generic Linear Programming solver using GLPK
 *
 * Keith O'Hara
 * 02/24/2017
 *
 * This version:
 * 02/26/2017
 */

#include "lp/trame_glpk.h"

//
// Dense setup; to be used with Armadillo memptr-based passing 
// Armadillo uses column-major ordering, as opposed to C-standard row-major ordering
int trame_glpk(int n_constr, int n_vars, double* obj, double* A, int model_opt_sense, 
               double* rhs, char* constr_sense, double* lb, double* ub, 
               double* objval, double* sol_mat_X, double* sol_mat_RC, 
               double* dual_mat_PI, double* dual_mat_SLACK)
{
    int success = 0;
    //
    int has_lb = (lb) ? 1 : 0;
    int has_ub = (ub) ? 1 : 0;

    int i,j;
    //
    // we need this for inputting constraint matrix values
    int* row_inds = malloc(sizeof(int) * (n_constr+1));
    for (i=1; i < n_constr + 1; i++) {
        row_inds[i] = i;
    }

    glp_prob *lp;
    lp = glp_create_prob();
    
    glp_term_out(GLP_OFF);

    if (model_opt_sense==0) { // minimize
        glp_set_obj_dir(lp, GLP_MIN);
    } else if (model_opt_sense==1) { // maximize
        glp_set_obj_dir(lp, GLP_MAX);
    } else {
        printf("unrecognized input for model_opt_sense; should be 0 or 1.\n");
        goto QUIT;
    }

    glp_add_rows(lp, n_constr);
    glp_add_cols(lp, n_vars);

    for (j=0; j < n_constr; j++) {
        if (constr_sense[j]=='<') {
            glp_set_row_bnds(lp, j+1, GLP_UP, 0.0, rhs[j]);
        } else if (constr_sense[i]=='>') {
            glp_set_row_bnds(lp, j+1, GLP_LO, rhs[j], 0.0);
        } else {
            glp_set_row_bnds(lp, j+1, GLP_FX, rhs[j], rhs[j]);
        }
    }

    for (i=0; i < n_vars; i++) {
        // c'x
        glp_set_obj_coef(lp, i+1, obj[i]);

        // set bounds on x
        if (has_lb && has_ub) {
            if (lb[i] == ub[i]) {
                glp_set_col_bnds(lp, i+1, GLP_FX, lb[i], ub[i]);
            } else {
                glp_set_col_bnds(lp, i+1, GLP_DB, lb[i], ub[i]);
            }
        } else if (has_lb && !has_ub) {
            glp_set_col_bnds(lp, i+1, GLP_LO, lb[i], 0.0);
        } else if (!has_lb && has_ub) {
            glp_set_col_bnds(lp, i+1, GLP_UP, 0.0, ub[i]);
        } else if (!has_lb && !has_ub) {
            glp_set_col_bnds(lp, i+1, GLP_FR, 0.0, 0.0);
        }
        
        // input constraint matrix elements
        glp_set_mat_col(lp, i+1, n_constr, row_inds, &A[i*n_constr-1]);
    }

    glp_simplex(lp, NULL);

    int glpk_success = glp_get_status(lp);

    if (glpk_success==GLP_OPT) {
        success = 1;
    }
    // retrieve optimum
    *objval = glp_get_obj_val(lp);
    
    for (i=0; i < n_vars; i++) {
        sol_mat_X[i] = glp_get_col_prim(lp, i+1);
        sol_mat_RC[i] = glp_get_col_dual(lp, i+1);
    }
    
    for (j=0; j < n_constr; j++) {
        dual_mat_PI[j] = glp_get_row_dual(lp, j+1);
        dual_mat_SLACK[j] = glp_get_row_prim(lp, j+1);
    }
    //

QUIT:
    free(row_inds);
    glp_delete_prob(lp);

    return success;
}

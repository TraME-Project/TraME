#ifndef _generic_lp_H
#define _generic_lp_H

#include <stdlib.h>
#include <stdio.h>
#include "gurobi_c.h"

int generic_LP_C(int rows, int cols, double* obj, double* A, int model_opt_sense, 
                 double* rhs, char* sense, double* Q, double* lb, double* ub,
                 double* objval, double* sol_mat_X, double* sol_mat_RC,
                 double* dual_mat_PI, double* dual_mat_SLACK);
int generic_LP_C_switch(int rows, int cols, double* obj, double* A, int model_opt_sense,
                 double* rhs, char* sense, double* Q, double* lb, double* ub,
                 double* objval, double* sol_mat_X, double* sol_mat_RC,
                 double* dual_mat_PI, double* dual_mat_SLACK);
int generic_LP_C_sparse(int rows, int cols, double* obj, int numnz,
                 int* vbeg, int* vind, double* vval, int model_opt_sense,
                 double* rhs, char* sense, double* Q, double* lb, double* ub,
                 double* objval, double* sol_mat_X, double* sol_mat_RC,
                 double* dual_mat_PI, double* dual_mat_SLACK);

#endif

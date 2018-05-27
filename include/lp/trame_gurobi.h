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
 * Generic Linear/Quadratic Programming solver using Gurobi
 *
 * Keith O'Hara
 * 05/08/2016
 *
 * This version:
 * 02/22/2017
 */

#ifndef _trame_gurobi_H
#define _trame_gurobi_H

#include <stdlib.h>
#include <stdio.h>
#include "gurobi_c.h"

int trame_gurobi(int rows, int cols, double* obj, double* A, int model_opt_sense, 
                 double* rhs, char* constr_dir, double* Q, double* lb, double* ub,
                 double* objval, double* sol_mat_X, double* sol_mat_RC,
                 double* dual_mat_PI, double* dual_mat_SLACK);
int trame_gurobi_switch(int rows, int cols, double* obj, double* A, int model_opt_sense, 
                        double* rhs, char* constr_dir, double* Q, double* lb, double* ub,
                        double* objval, double* sol_mat_X, double* sol_mat_RC,
                        double* dual_mat_PI, double* dual_mat_SLACK);
int trame_gurobi_sparse(int rows, int cols, double* obj, int numnz,
                        int* vbeg, int* vind, double* vval, int model_opt_sense,
                        double* rhs, char* constr_dir, double* Q, double* lb, double* ub,
                        double* objval, double* sol_mat_X, double* sol_mat_RC,
                        double* dual_mat_PI, double* dual_mat_SLACK);

#endif

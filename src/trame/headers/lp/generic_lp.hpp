/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 the TraME Team:
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
 * Generic Linear/Quadratic Programming solver using Gurobi
 *
 * Keith O'Hara
 * 05/08/2016
 */

#ifndef _generic_lp_HPP
#define _generic_lp_HPP

#include "../misc/TRAME_OPTIONS.hpp"

#if (defined(WIN32) || defined(_WIN32) || defined(TRAME_USE_GUROBI_C)) && !defined(TRAME_PREDEF_GUROBI_C)
    #define TRAME_PREDEF_GUROBI_C
#endif

/*bool generic_LP(int k, int n, double *obj, double* A, int modelSense, double* rhs, char* sense, double* Q, 
                double* lb, double* ub, double* start, double& objval, arma::mat& sol_mat, arma::mat& dual_mat);

bool generic_LP(int k, int n, double *obj, int numnz, int* vbeg, int* vind, double* vval, int modelSense, double* rhs, char* sense, double* Q, 
                double* lb, double* ub, double* start, double& objval, arma::mat& sol_mat, arma::mat& dual_mat);*/

bool generic_LP(int k, int n, double *obj, double* A, int modelSense, double* rhs, char* sense, double* Q, 
                double* lb, double* ub, double* start, double& objval, double* sol_mat_X, double* sol_mat_RC, double* dual_mat_PI, double* dual_mat_SLACK);

bool generic_LP(int k, int n, double *obj, int numnz, int* vbeg, int* vind, double* vval, int modelSense, double* rhs, char* sense, double* Q, 
                double* lb, double* ub, double* start, double& objval, double* sol_mat_X, double* sol_mat_RC, double* dual_mat_PI, double* dual_mat_SLACK);

#endif

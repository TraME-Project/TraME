/*
 * Generic Linear/Quadratic Programming solver using Gurobi
 *
 * Keith O'Hara
 * 05/08/2016
 */

// Keith: put this in a separate file at some point
#include "../misc/TRAME_OPTIONS.hpp"

//#if (defined(WIN32) || defined(_WIN32) || defined(TRAME_USE_GUROBI_C)) && !defined(TRAME_PREDEF_GUROBI_C)
#if defined(WIN32) || defined(_WIN32) || defined(TRAME_USE_GUROBI_C)
    #define TRAME_PREDEF_GUROBI_C
#endif

#include <RcppArmadillo.h>

bool generic_LP(int k, int n, double *obj, double* A, int modelSense, double* rhs, char* sense, double* Q, 
                double* lb, double* ub, double* start, double& objval, arma::mat& sol_mat, arma::mat& dual_mat);

bool generic_LP(int k, int n, double *obj, int numnz, int* vbeg, int* vind, double* vval, int modelSense, double* rhs, char* sense, double* Q, 
                double* lb, double* ub, double* start, double& objval, arma::mat& sol_mat, arma::mat& dual_mat);

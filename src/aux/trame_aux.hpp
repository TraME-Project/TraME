/*
 * Find indices the correspond to maximum values
 *
 * Keith O'Hara
 * 05/08/2016
 */

#ifndef _trame_aux_HPP
#define _trame_aux_HPP

#include <RcppArmadillo.h>

arma::uvec which_max(const arma::mat* X, int which_dim);
arma::uvec uvec_linspace(int a, int b);

double pbeta (double x, double* fn_pars);
arma::vec pbeta (arma::vec x, double* fn_pars);
double qbeta (double x, double* fn_pars);
arma::vec qbeta (arma::vec x, double* fn_pars);
double iqbeta (double x, double* fn_pars);
arma::vec iqbeta (arma::vec x, double* fn_pars);
double dbeta (double x, double* fn_pars);
arma::vec dbeta (arma::vec x, double* fn_pars);

#endif

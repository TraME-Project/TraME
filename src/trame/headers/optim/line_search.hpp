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
 * Moré and Thuente line search
 *
 * Based on MINPACK fortran code and Dianne P. O'Leary's Matlab translation of MINPACK
 *
 * Keith O'Hara
 * 01/03/2017
 *
 * This version:
 * 01/11/2017
 */

#ifndef _line_search_HPP
#define _line_search_HPP

double line_search_mt(double step, arma::vec& x, arma::vec& grad, const arma::vec& direc, double* wolfe_cons_1_inp, double* wolfe_cons_2_inp, std::function<double (const arma::vec& vals_inp, arma::vec& grad, void* opt_data)> opt_objfn, void* opt_data);
double mt_sup_norm(double a, double b, double c);
int mt_step(double& st_best, double& f_best, double& d_best, double& st_other, double& f_other, double& d_other, double& step, double& f_step, double& d_step, bool& bracket, double step_min, double step_max);

#endif

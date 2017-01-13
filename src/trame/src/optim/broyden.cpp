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
 * Broyden's method for solving systems of nonlinear equations
 *
 * Keith O'Hara
 * 01/03/2017
 *
 * This version:
 * 01/11/2017
 */

#include "trame.hpp"

bool trame::broyden(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data)
{
    // notation: 'p' stands for '+1'.
    //
    bool success = false;
    int max_iter = 1000;
    double err_tol = 1e-08;
    //
    int n_vals = init_out_vals.n_elem;

    arma::vec x = init_out_vals;

    arma::mat B = arma::eye(n_vals,n_vals); // initial approx. to (inverse) Jacobian
    //
    // initialization
    double t_init = 1;

    arma::vec f_val = opt_objfn(x,opt_data);

    arma::vec d = - B*f_val;
    double t = 1;

    arma::vec x_p = x + t*d;
    arma::vec f_val_p = opt_objfn(x_p,opt_data);

    arma::vec s = x_p - x;
    arma::vec y = f_val_p - f_val;

    // update B
    B += (s - B*y) * y.t() / arma::dot(y,y);

    f_val = f_val_p;
    //
    // begin loop
    int iter = 0;
    double err = 2*err_tol;

    while (err > err_tol && iter < max_iter) {
        iter++;

        d = - B*f_val;

        x_p = x + t*d;
        f_val_p = opt_objfn(x_p,opt_data);

        s = x_p - x;
        y = f_val_p - f_val;
        // update B
        B += (s - B*y) * y.t() / arma::dot(y,y);
        //
        err = arma::as_scalar(arma::sum(arma::abs(f_val_p)));

        x = x_p;
        f_val = f_val_p;
    }
    //
    if (err <= err_tol && iter <= max_iter) {
        init_out_vals = x_p;
        success = true;
    }
    //
    return success;
}

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
 * BFGS method for quasi-Newton-based non-linear optimization
 *
 * Keith O'Hara
 * 12/23/2016
 *
 * This version:
 * 02/07/2017
 */

#include "trame.hpp"

bool trame::bfgs(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data)
{
    bool success = bfgs_int(init_out_vals,opt_objfn,opt_data,NULL,NULL,NULL);
    //
    return success;
}

bool trame::bfgs(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, double& value_out)
{
    bool success = bfgs_int(init_out_vals,opt_objfn,opt_data,&value_out,NULL,NULL);
    //
    return success;
}

bool trame::bfgs(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, double& value_out, double err_tol_inp, int max_iter_inp)
{
    bool success = bfgs_int(init_out_vals,opt_objfn,opt_data,&value_out,&err_tol_inp,&max_iter_inp);
    //
    return success;
}

bool trame::bfgs_int(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, double* value_out, double* err_tol_inp, int* max_iter_inp)
{
    // notation: 'p' stands for '+1'.
    //
    bool success = false;
    int max_iter = (max_iter_inp) ? *max_iter_inp : 1000;
    double err_tol = (err_tol_inp) ? *err_tol_inp : 1e-08;

    double wolfe_cons_1 = 1E-03; // line search tuning parameters
    double wolfe_cons_2 = 0.90;
    //
    int n_vals = init_out_vals.n_elem;
    arma::vec x = init_out_vals;

    arma::mat W = arma::eye(n_vals,n_vals); // initial approx. to (inverse) Hessian 
    //
    // initialization
    arma::vec grad(n_vals); // gradient
    opt_objfn(x,&grad,opt_data);

    double err = arma::accu(arma::abs(grad));
    if (err <= err_tol) {
        return true;
    }
    // if ||gradient(initial values)|| > tolerance, then continue
    double t_init = 1; // initial line search value
    arma::vec d = - W*grad; // direction

    arma::vec x_p = x, grad_p = grad;
    line_search_mt(t_init, x_p, grad_p, d, &wolfe_cons_1, &wolfe_cons_2, opt_objfn, opt_data);

    err = arma::accu(arma::abs(grad_p)); // check updated values
    if (err <= err_tol) {
        init_out_vals = x_p;
        return true;
    }
    // if ||gradient(x_p)|| > tolerance, then continue
    arma::vec s = x_p - x;
    arma::vec y = grad_p - grad;

    // update W
    double W_denom_term = arma::dot(y,s);

    arma::mat W_term_1 = s*y.t()*W + W*y*s.t();
    arma::mat W_term_2 = (1.0 + arma::dot(y,W*y) / W_denom_term) * s*s.t();

    W -= (W_term_1 - W_term_2) / W_denom_term; // update inverse Hessian approximation

    grad = grad_p;
    //
    // begin loop
    int iter = 0;

    while (err > err_tol && iter < max_iter) {
        iter++;
        //
        d = - W*grad;
        line_search_mt(t_init, x_p, grad_p, d, &wolfe_cons_1, &wolfe_cons_2, opt_objfn, opt_data);
        
        err = arma::accu(arma::abs(grad_p));
        if (err <= err_tol) {
            break;
        }
        // if ||gradient(x_p)|| > tolerance, then continue
        s = x_p - x;
        y = grad_p - grad;
        // update W
        W_denom_term = arma::dot(y,s);

        if (std::abs(W_denom_term) < 1E-08) {
            W_denom_term = 1E-08;
        }

        W_term_1 = s*y.t()*W + W*y*s.t();
        W_term_2 = (1.0 + arma::dot(y,W*y) / W_denom_term) * s*s.t();

        W -= (W_term_1 - W_term_2) / W_denom_term;
        //
        x = x_p;
        grad = grad_p;
    }
    //
    if (err <= err_tol && iter <= max_iter) {
        init_out_vals = x_p;
        success = true;

        if (value_out) {
            *value_out = opt_objfn(init_out_vals,NULL,opt_data);
        }
    } else {
        printf("bfgs failure: max_iter reached before convergence could be achieved.\n");
        printf("bfgs failure: best guess:\n");
        arma::cout << x_p.t() << arma::endl;
        std::cout << "error: " << err << std::endl;
    }
    //
    return success;
}

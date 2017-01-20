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
 * 01/19/2017
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
    arma::vec f_val = opt_objfn(x,opt_data);

    double err = arma::as_scalar(arma::sum(arma::abs(f_val)));
    if (err <= err_tol) {
        return true;
    }
    //
    arma::vec d = - B*f_val;

    arma::vec x_p = x + d;
    arma::vec f_val_p = opt_objfn(x_p,opt_data);

    err = arma::as_scalar(arma::sum(arma::abs(f_val_p)));
    if (err <= err_tol) {
        init_out_vals = x_p;
        return true;
    }
    //
    arma::vec s = x_p - x;
    arma::vec y = f_val_p - f_val;

    B += (s - B*y) * y.t() / arma::dot(y,y); // update B

    f_val = f_val_p;
    //
    // begin loop
    int iter = 0;

    while (err > err_tol && iter < max_iter) {
        iter++;
        //
        d = - B*f_val;

        x_p = x + d;
        f_val_p = opt_objfn(x_p,opt_data);

        err = arma::as_scalar(arma::sum(arma::abs(f_val_p)));

        if (err <= err_tol) {
            break;
        }
        //
        s = x_p - x;
        y = f_val_p - f_val;
        
        B += (s - B*y) * y.t() / arma::dot(y,y); // update B
        //
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

// derivative-free method of Li and Fukushima (2000)

bool trame::broyden_df(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data)
{
    // notation: 'p' stands for '+1'.
    //
    bool success = false;
    int max_iter = 1000;
    double err_tol = 1e-08;
    
    double rho = 0.9, sigma_1 = 0.001, sigma_2 = 0.001; // tuning parameters
    //
    int n_vals = init_out_vals.n_elem;

    arma::vec x = init_out_vals;

    arma::mat B = arma::eye(n_vals,n_vals); // initial approx. to (inverse) Jacobian
    //
    // initialization
    arma::vec f_val = opt_objfn(x,opt_data);
    double err = arma::as_scalar(arma::sum(arma::abs(f_val)));

    if (err <= err_tol) {
        return true;
    }

    double Fx = arma::norm(f_val,2);
    //
    arma::vec d = -f_val; // step 1

    arma::vec f_val_p = opt_objfn(x + d,opt_data);
    err = arma::as_scalar(arma::sum(arma::abs(f_val_p)));

    if (err <= err_tol) {
        init_out_vals = x + d;
        return true;
    }
    //
    double lambda;
    double Fx_p = arma::norm(f_val_p,2);

    if (Fx_p <= rho*Fx - sigma_2*arma::dot(d,d)) { // step 2
        lambda = 1.0;
    } else {
        lambda = df_proc_1(x,d,sigma_1,0,opt_objfn,opt_data); // step 3
    }
    //
    arma::vec x_p = x + lambda*d; // step 4

    arma::vec s = x_p - x;
    arma::vec y = f_val_p - f_val;

    B += (y - B*s) * s.t() / arma::dot(s,s); // step 5
    //
    x = x_p;
    f_val = f_val_p;
    Fx = Fx_p;
    //
    // begin loop
    int iter = 0;

    while (err > err_tol && iter < max_iter) {
        iter++;

        d = arma::solve(B,-f_val);
        f_val_p = opt_objfn(x + d,opt_data);

        err = arma::as_scalar(arma::sum(arma::abs(f_val)));

        if (err <= err_tol) {
            break;
        }
        //
        Fx_p = arma::norm(f_val_p,2);

        if (Fx_p <= rho*Fx - sigma_2*arma::dot(d,d)) {
            lambda = 1.0;
        } else {
            lambda = df_proc_1(x,d,sigma_1,iter,opt_objfn,opt_data);
        }
        //
        x_p = x + lambda*d;

        arma::vec s = x_p - x;
        arma::vec y = f_val_p - f_val;

        B += (y - B*s) * s.t() / arma::dot(s,s);
        //
        x = x_p;
        f_val = f_val_p;
        Fx = Fx_p;
    }
    //
    if (err <= err_tol && iter <= max_iter) {
        init_out_vals = x_p;
        success = true;
    }
    //
    return success;
}

double trame::df_eta(int k)
{
    return 1.0 / (k*k);
}

double trame::df_proc_1(const arma::vec& x_vals, const arma::vec& direc, double sigma_1, int k, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data)
{
    double beta = 0.9, lambda = 1.0;
    double eta_k = df_eta(k);
    //
    // check: || F(x_k + lambda*d_k) || <= ||F(x_k)||*(1+eta_k) - sigma_1*||lambda*d_k||^2
    double Fx = arma::norm(opt_objfn(x_vals,opt_data),2);
    double Fx_p = arma::norm(opt_objfn(x_vals + lambda*direc,opt_data),2);
    double direc_norm2 = arma::dot(direc,direc);

    double term_2 = sigma_1*std::pow(lambda,2)*direc_norm2;
    double term_3 = eta_k*Fx;
    //
    if (Fx_p <= Fx - term_2 + term_3) {
        return lambda;
    }
    //
    int iter = 0;
    
    while (1) {
        iter++;
        lambda *= beta; // lambda_i = beta^i;
        //
        Fx_p = arma::norm(opt_objfn(x_vals + lambda*direc,opt_data),2);
        term_2 = sigma_1*std::pow(lambda,2)*direc_norm2;

        if (Fx_p <= Fx - term_2 + term_3) {
            break;
        }
    }
    //
    return lambda;
}

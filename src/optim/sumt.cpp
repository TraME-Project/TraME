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
 * Sequential unconstrained minimization technique (SUMT)
 *
 * Keith O'Hara
 * 01/15/2016
 *
 * This version:
 * 01/17/2017
 */

#include "trame.hpp"

bool trame::sumt(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data,
                 std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* constr_data)> constr_fn, void* constr_data,
                 double* value_out)
{
    // notation: 'p' stands for '+1'.
    //
    bool success = false;
    int max_iter = 1000;
    double err_tol = 1e-08;
    //
    double eta = 2.0; // growth of penalty parameter
    arma::vec x = init_out_vals;
    //
    // lambda function that combines the objective function with the constraints
    std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* sumt_data)> sumt_objfn = [opt_objfn, opt_data, constr_fn, constr_data] (const arma::vec& vals_inp, arma::vec* grad, void* sumt_data) -> double {
        sumt_struct *d = reinterpret_cast<sumt_struct*>(sumt_data);
        double c_pen = d->c_pen;
        //
        int n_vals = vals_inp.n_elem;
        arma::vec grad_obj(n_vals), grad_constr(n_vals);
        //
        //double ret = opt_objfn(vals_inp,grad_obj,opt_data) + c_pen*(std::pow(std::max(0.0,constr_fn(vals_inp,grad_constr,constr_data)),2) / 2.0);
        double ret;
        double constr_val = constr_fn(vals_inp,&grad_constr,constr_data);

        if (constr_val < 0.0) {
            ret = opt_objfn(vals_inp,&grad_obj,opt_data);

            if (grad) {
                *grad = grad_obj;
            }
        } else {
            ret = opt_objfn(vals_inp,&grad_obj,opt_data) + c_pen*(constr_val*constr_val / 2.0);

            if (grad) {
                *grad = grad_obj + c_pen*grad_constr;
            }
        }
        //
        return ret;
    };
    //
    // initialization
    double c_pen = 1.0;

    sumt_struct sumt_data;
    sumt_data.c_pen = c_pen;

    arma::vec x_p = x;
    //
    // begin loop
    int iter = 0;
    double err = 2*err_tol;

    while (err > err_tol && iter < max_iter) {
        iter++;
        //
        bfgs(x_p,sumt_objfn,&sumt_data);
        err = arma::norm(x_p - x,2);
        //
        sumt_data.c_pen = eta*sumt_data.c_pen; // increase penalization parameter value
        x = x_p;
    }
    //
    if (err <= err_tol && iter <= max_iter) {
        init_out_vals = x_p;
        success = true;

        if (value_out) {
            *value_out = opt_objfn(x_p,NULL,opt_data);
        }
    } else {
        printf("sumt failure: max_iter reached before convergence could be achieved.\n");
        printf("sumt failure: best guess:\n");
        arma::cout << x_p.t() << arma::endl;
        std::cout << "error: " << err << std::endl;
    }
    //
    return success;
}

bool trame::sumt(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data,
                 std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* constr_data)> constr_fn, void* constr_data)
{
    bool success = sumt(init_out_vals,opt_objfn,opt_data,constr_fn,constr_data,NULL);

    return success;
}

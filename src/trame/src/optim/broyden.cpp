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

arma::vec unit_vec(int j, int n)
{
    arma::vec ret = arma::zeros(n,1);
    ret(j) = 1;

    return ret;
}

bool broyden(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data)
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

bool broyden_mt(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data)
{
    // notation: 'p' stands for '+1'.
    //
    bool success = false;
    int max_iter = 1000;
    double err_tol = 1e-08;
    double wolfe_cons_1 = 1E-03;
    double wolfe_cons_2 = 0.50;
    //
    int n_vals = init_out_vals.n_elem;

    arma::vec x = init_out_vals;

    arma::mat B = arma::eye(n_vals,n_vals); // initial approx. to (inverse) Jacobian
    //
    /*std::function<double (const arma::vec& vals_inp, arma::vec& grad, void* opt_data)> ls_objfn = [opt_objfn, &B] (const arma::vec& vals_inp, arma::vec& grad, void* opt_data) -> double {
        int n_vals = vals_inp.n_elem;
        double eps_diff = 1e-04;

        arma::vec Fx = opt_objfn(vals_inp,opt_data);
        double ret = arma::dot(Fx,Fx)/2.0;
        std::cout << "lambda ret: " << ret << std::endl;

        arma::vec Fx_p, Fx_m;
        for (int jj=0; jj < n_vals; jj++) {
            arma::vec Fx_p = opt_objfn(vals_inp + eps_diff*unit_vec(jj,n_vals),opt_data);
            arma::vec Fx_m = opt_objfn(vals_inp - eps_diff*unit_vec(jj,n_vals),opt_data);

            grad(jj) = (arma::dot(Fx_p,Fx_p)/2.0 - arma::dot(Fx_m,Fx_m)/2.0) / (2*eps_diff);
        }
        arma::cout << "lambda grad: " << grad.t() << arma::endl;
        //
        return ret;
    };*/
    std::function<double (const arma::vec& vals_inp, arma::vec& grad, void* opt_data)> ls_objfn = [opt_objfn, &B] (const arma::vec& vals_inp, arma::vec& grad, void* opt_data) -> double {
        int n_vals = vals_inp.n_elem;
        double eps_diff = 1e-04;

        arma::vec Fx = opt_objfn(vals_inp,opt_data);
        double ret = arma::dot(Fx,Fx)/2.0;
        std::cout << "lambda ret: " << ret << std::endl;

        arma::vec Fx_p, Fx_m;
        for (int jj=0; jj < n_vals; jj++) {
            arma::vec Fx_p = opt_objfn(vals_inp + eps_diff*unit_vec(jj,n_vals),opt_data);
            arma::vec Fx_m = opt_objfn(vals_inp - eps_diff*unit_vec(jj,n_vals),opt_data);
            
            grad(jj) = (arma::dot(Fx_p,Fx_p)/2.0 - arma::dot(Fx_m,Fx_m)/2.0) / (2*eps_diff);
        }
        arma::cout << "lambda grad: " << grad.t() << arma::endl;
        //
        return ret;
    };
    //
    // initialization
    double t_init = 1;

    arma::vec f_val = opt_objfn(x,opt_data);

    arma::vec d = - B*f_val;

    arma::vec x_p = x, grad_mt(n_vals);
    double t = line_search_mt(t_init, x_p, grad_mt, d, &wolfe_cons_1, &wolfe_cons_2, ls_objfn, opt_data);
    std::cout << "t: " << t << std::endl;
    arma::cout << x_p << arma::endl;

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
        t = line_search_mt(t_init, x_p, grad_mt, d, &wolfe_cons_1, &wolfe_cons_2, ls_objfn, opt_data);

        f_val_p = opt_objfn(x_p,opt_data);

        s = x_p - x;
        y = f_val_p - f_val;
        // update B
        B += (s - B*y) * y.t() / arma::dot(y,y);
        //
        err = arma::as_scalar(arma::sum(arma::abs(f_val_p)));

        x = x_p;
        f_val = f_val_p;

        std::cout << "f_val_p:" << f_val_p << std::endl;
    }
    //
    if (err <= err_tol && iter <= max_iter) {
        init_out_vals = x_p;
        success = true;
    }
    //
    return success;
}

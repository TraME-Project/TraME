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
 * Generic input to optimization routines
 *
 * Keith O'Hara
 * 01/11/2017
 *
 * This version:
 * 01/17/2017
 */

#include "trame.hpp"

bool trame::generic_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec& grad, void* opt_data)> opt_objfn, void* opt_data)
{
    bool success = bfgs(init_out_vals,opt_objfn,opt_data);
    //
    return success;
}

// box constraints

bool trame::generic_optim(arma::vec& init_out_vals, const arma::vec& lower_bounds, const arma::vec& upper_bounds, std::function<double (const arma::vec& vals_inp, arma::vec& grad, void* opt_data)> opt_objfn, void* opt_data)
{
    std::function<double (const arma::vec& vals_inp, arma::vec& grad, void* box_data)> box_objfn = [opt_objfn, lower_bounds, upper_bounds] (const arma::vec& vals_inp, arma::vec& grad, void* opt_data) -> double {
        //
        arma::vec vals_inv_trans = logit_inv_trans(vals_inp,lower_bounds,upper_bounds);
        //
        double ret = opt_objfn(vals_inv_trans,grad,opt_data);
        arma::mat jacob_correct = jacob_matrix_logit(vals_inp,lower_bounds,upper_bounds);

        grad = jacob_correct.t() * grad; // correct gradient for transformation
        //
        return ret;
    };
    
    arma::vec initial_vals = logit_trans(init_out_vals,lower_bounds,upper_bounds);
    //
    bool success = bfgs(initial_vals,box_objfn,opt_data);
    init_out_vals = logit_inv_trans(initial_vals,lower_bounds,upper_bounds);
    //
    return success;
}

arma::mat trame::jacob_matrix_logit(const arma::vec& trans_vals, const arma::vec& lower_bounds, const arma::vec& upper_bounds)
{
    int n_vals = trans_vals.n_elem;
    arma::mat ret_mat = arma::zeros(n_vals,n_vals);

    for (int i=0; i < n_vals; i++) {
        ret_mat(i,i) = std::exp(trans_vals(i))*(upper_bounds(i) - lower_bounds(i)) / std::pow(std::exp(trans_vals(i)) + 1,2);
    }
    //
    return ret_mat;
}

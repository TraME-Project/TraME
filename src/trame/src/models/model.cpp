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
 * general model class
 * template specialization
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * This version:
 * 11/27/2016
 */

#include "trame.hpp"

namespace trame
{
template<>
bool model<logit>::mme(const arma::mat& mu_hat, arma::mat& theta_hat)
{
    bool success = false;
    //
    double xtol_rel = (xtol_rel_inp) ? *xtol_rel_inp : 1E-04;
    int max_eval = (max_eval_inp) ? *max_eval_inp : 1E05;

    arma::vec theta_0;
    init_param(theta_0);

    arma::mat dtheta_Psi;
    dparam(NULL,dtheta_Psi);

    build_market_TU(theta_0);

    arma::mat kron_term = dtheta_Psi;
    arma::mat C_hat = arma::vectorise(arma::vectorise(mu_hat) * kron_term);

    //trame::mfe<trame::mmf> mkt_obj = build_market(theta_0);

    arma::mat C_hat = Phi_k(mu_hat);
    //
    double total_mass = arma::accu(n);

    arma::vec p = n / total_mass;
    arma::vec q = m / total_mass;

    arma::mat IX = arma::ones(nbX,1);
    arma::mat tIY = arma::ones(1,nbY);

    arma::mat f = p * tIY;
    arma::mat g = IX * q.t();

    arma::mat Pi_hat = mu_hat / total_mass;
    arma::mat v = arma::zeros(1,nbY);

    arma::mat phi_xy = arma::reshape(phi_xyk_aux,nbX*nbY,nbParams);
    //
    // add optimization data
    trame_nlopt_opt_data opt_data;

    opt_data.mme_woregal.nbX = nbX;
    opt_data.mme_woregal.nbY = nbY;

    opt_data.mme_woregal.dX = dX;
    opt_data.mme_woregal.dY = dY;
    
    opt_data.mme_woregal.max_iter_ipfp = max_iter_ipfp;
    opt_data.mme_woregal.tol_ipfp = tol_ipfp;

    opt_data.mme_woregal.sigma = sigma;

    opt_data.mme_woregal.p = p;
    opt_data.mme_woregal.q = q;

    opt_data.mme_woregal.IX = IX;
    opt_data.mme_woregal.tIY = tIY;

    opt_data.mme_woregal.f = f;
    opt_data.mme_woregal.g = g;

    opt_data.mme_woregal.v = v;
    opt_data.mme_woregal.Pi_hat = Pi_hat;

    opt_data.mme_woregal.phi_xy = phi_xy; // should be (nbX*nbY) x (nbParams)
    //
    arma::mat A0 = arma::zeros(dX*dY,1);
    std::vector<double> sol_vec = arma::conv_to< std::vector<double> >::from(A0);
    double obj_val = 0;

    bool success = generic_nlopt(dX*dY,sol_vec,obj_val,NULL,NULL,trame::affinity::mme_woregul_opt_objfn,opt_data);
    //
    theta_hat = arma::conv_to< arma::mat >::from(sol_vec);
    val_ret = obj_val;
    //
    return success;
}
}

/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
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
 * Constant Elasticity of Substitution (CES) Marriage Matching Functions (MMFs) class
 * Corresponds to the ETU transfers class
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 03/27/2017
 */

#include "trame.hpp"

void 
trame::mmfs::ces::build(const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = alpha_inp.n_rows;
    nbY = alpha_inp.n_cols;
    nbParams = 3*nbX*nbY;

    alpha = alpha_inp;
    gamma = gamma_inp;
    tau   = tau_inp;
    kappa = - 1.0 / tau_inp;

    aux_alpha = - elem_div(alpha_inp, tau_inp);
    aux_gamma = - elem_div(gamma_inp, tau_inp);

    aux_alpha_exp = arma::exp(- elem_div(alpha_inp, tau_inp));
    aux_gamma_exp = arma::exp(- elem_div(gamma_inp, tau_inp));
}

void 
trame::mmfs::ces::trans()
{
    int nbX_temp = nbX;

    nbX = nbY;
    nbY = nbX_temp;
    //
    arma::mat alpha_temp = alpha;

    alpha = gamma.t();
    gamma = alpha_temp.t();
    tau = tau.t();
    kappa = kappa.t();
    //
    arma::mat aux_alpha_temp = aux_alpha;
    arma::mat aux_alpha_exp_temp = aux_alpha_exp;

    arma::mat aux_alpha = aux_gamma.t();
    arma::mat aux_gamma = aux_alpha_temp.t();

    arma::mat aux_alpha_exp = aux_gamma_exp.t();
    arma::mat aux_gamma_exp = aux_alpha_exp_temp.t();
}

//
// MFE-related functions

// matching function

arma::mat 
trame::mmfs::ces::M(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    arma::mat ret = this->M(a_xs,b_ys,NULL,NULL);
    //
    return ret;
}

arma::mat 
trame::mmfs::ces::M(const arma::mat& a_xs, const arma::mat& b_ys, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = arma::exp(aux_alpha(x_ind,y_ind) + elem_prod(kappa(x_ind,y_ind), arma::log(a_xs)));
    arma::mat term_2 = arma::exp(aux_gamma(x_ind,y_ind) + arma::trans(elem_prod(arma::trans(kappa(x_ind,y_ind)), arma::log(b_ys))));

    arma::mat ret = arma::exp( - tau(x_ind,y_ind) % arma::log((term_1 + term_2)/2) ); 
    //
    return ret;
}

arma::mat 
trame::mmfs::ces::M(const double& a_xs, const arma::mat& b_ys, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, (int) nbX-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, (int) nbY-1);
    //
    arma::mat term_1 = arma::exp(aux_alpha(x_ind,y_ind) + kappa(x_ind,y_ind) * std::log(a_xs));
    arma::mat term_2 = arma::exp(aux_gamma(x_ind,y_ind) + arma::trans(elem_prod(arma::trans(kappa(x_ind,y_ind)), arma::log(b_ys))));

    arma::mat ret = arma::exp( - tau(x_ind,y_ind) % arma::log((term_1 + term_2)/2) ); 
    //
    return ret;
}

arma::mat 
trame::mmfs::ces::M(const arma::mat& a_xs, const double& b_ys, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, (int) nbX-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, (int) nbY-1);
    //
    arma::mat term_1 = arma::exp(aux_alpha(x_ind,y_ind) + elem_prod(kappa(x_ind,y_ind), arma::log(a_xs)));
    arma::mat term_2 = arma::exp(aux_gamma(x_ind,y_ind) + kappa(x_ind,y_ind) * std::log(b_ys));

    arma::mat ret = arma::exp( - tau(x_ind,y_ind) % arma::log((term_1 + term_2)/2) );
    //
    return ret;
}

//

arma::mat 
trame::mmfs::ces::dmu_x0(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    arma::mat term_1 = arma::exp(aux_alpha + elem_prod(kappa, arma::log(a_xs)));
    arma::mat term_2 = arma::exp(aux_gamma + arma::trans(elem_prod(arma::trans(kappa), arma::log(b_ys))));

    arma::mat mu_s = arma::exp( - tau % arma::log((term_1 + term_2)/2) ); 

    arma::mat ret = aux_alpha_exp % arma::exp( elem_prod(1 - kappa, arma::log(mu_s / a_xs)) ) / 2.0;
    //
    return ret;
}

arma::mat 
trame::mmfs::ces::dmu_0y(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    arma::mat term_1 = arma::exp(aux_alpha + elem_prod(kappa, arma::log(a_xs)));
    arma::mat term_2 = arma::exp(aux_gamma + arma::trans(elem_prod(arma::trans(kappa), arma::log(b_ys))));

    arma::mat mu_s = arma::exp( - tau % arma::log((term_1 + term_2)/2) ); 

    arma::mat ret = aux_gamma_exp % arma::trans(arma::exp( elem_prod(arma::trans(1 - kappa), arma::log(arma::trans(mu_s) / a_xs)) )) / 2.0;
    //
    return ret;
}

arma::mat 
trame::mmfs::ces::dparams_M(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    return this->dparams_M(a_xs,b_ys,NULL);
}

arma::mat 
trame::mmfs::ces::dparams_M(const arma::mat& a_xs, const arma::mat& b_ys, const arma::mat* delta_params_M)
const
{
    arma::mat a_xs_kappa = arma::exp(elem_prod(kappa, arma::log(a_xs))); // mux0s^kappa
    arma::mat b_ys_kappa = arma::exp(arma::trans(elem_prod(arma::trans(kappa), arma::log(b_ys)))); // t(mu0ys^tkappa))

    arma::mat term_1 = aux_alpha_exp % a_xs_kappa;
    arma::mat term_2 = aux_gamma_exp % b_ys_kappa;

    arma::mat mu_s = arma::exp( - tau % arma::log((term_1 + term_2)/2) );\
    arma::mat mu_s_kappa = arma::exp(elem_prod(1 - kappa, arma::log(mu_s)));

    arma::mat numer = (mu_s / kappa) % ( term_1 % arma::log(a_xs) + aux_gamma_exp % arma::trans(arma::log(b_ys)) % b_ys_kappa - arma::log(mu_s) );
    arma::mat denom = term_1 + term_2;

    arma::mat der_1 = mu_s_kappa % a_xs_kappa / 2.0;
    arma::mat der_2 = mu_s_kappa % b_ys_kappa / 2.0;
    arma::mat der_3 = numer / denom;

    arma::mat ret;

    if (delta_params_M) {
        arma::mat delta_params_1 = arma::reshape((*delta_params_M).rows(0,nbX*nbY-1),nbX,nbY);
        arma::mat delta_params_2 = arma::reshape((*delta_params_M).rows(nbX*nbY,2*nbX*nbY-1),nbX,nbY);
        arma::mat delta_params_3 = arma::reshape((*delta_params_M).rows(2*nbX*nbY,3*nbX*nbY-1),nbX,nbY);

        ret = arma::vectorise(delta_params_1 % der_1 + delta_params_2 % der_2 + delta_params_3 % der_3);
    } else {
        ret = arma::join_rows( arma::join_rows( arma::diagmat(arma::vectorise(der_1)), arma::diagmat(arma::vectorise(der_2)) ), arma::diagmat(arma::vectorise(der_3)) );
    }
    //
    return ret;
}

arma::mat 
trame::mmfs::ces::Mx0(const arma::mat& a_x)
const
{
    return a_x;
}

arma::mat 
trame::mmfs::ces::M0y(const arma::mat& b_y)
const
{
    return b_y;
}

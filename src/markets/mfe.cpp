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
 * Matching Function Equilibrium (MFE) market
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 03/14/2017
 */

#include "trame.hpp"

namespace trame
{

// builds

template<>
void mfe<etu>::build_ETU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, double* sigma_inp, bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    need_norm = need_norm_inp;

    outsideOption = (need_norm_inp) ? false : true;
    sigma = (sigma_inp) ? *sigma_inp : 1.0;

    trans_obj.build(alpha_inp,gamma_inp,tau_inp,need_norm_inp);
    //
}

template<>
void mfe<ltu>::build_LTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, double* sigma_inp, bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    need_norm = need_norm_inp;

    outsideOption = (need_norm_inp) ? false : true;
    sigma = (sigma_inp) ? *sigma_inp : 1.0;
    
    trans_obj.build(lambda_inp,phi_inp/sigma,need_norm_inp);
    //
}

template<>
void mfe<ntu>::build_NTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, double* sigma_inp, bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    need_norm = need_norm_inp;

    outsideOption = (need_norm_inp) ? false : true;
    sigma = (sigma_inp) ? *sigma_inp : 1.0;

    trans_obj.build(alpha_inp/sigma,gamma_inp/sigma,need_norm_inp);
    //
}

template<>
void mfe<trame::tu>::build_TU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, double* sigma_inp, bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    need_norm = need_norm_inp;

    outsideOption = (need_norm_inp) ? false : true;
    sigma = (sigma_inp) ? *sigma_inp : 1.0;

    trans_obj.build(phi_inp/sigma,need_norm_inp);
    //
}

// ipfp-related functions

template<>
arma::vec mfe<ntu>::marg_x_inv(const arma::mat& B_ys, arma::uvec* xs)
const
{
    arma::uvec temp_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    //
    arma::vec a_NTU = n.elem(temp_ind);
    arma::mat B_NTU = arma::trans( elem_prod(arma::trans(trans_obj.aux_gamma_exp.rows(temp_ind)/trans_obj.aux_alpha_exp.rows(temp_ind)), B_ys) );
    arma::mat C_NTU = trans_obj.aux_alpha_exp.rows(temp_ind);

    arma::vec the_a_xs = inv_pwa(a_NTU, B_NTU, C_NTU, 1.0);
    //
    return the_a_xs;
}

template<>
arma::vec mfe<tu>::marg_x_inv(const arma::mat& B_ys, arma::uvec* xs)
const
{
    arma::uvec temp_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    //
    arma::mat sqrt_A_xs;
    arma::mat sqrt_B_ys = arma::sqrt(B_ys);

    if (!need_norm) {
        arma::mat b = (trans_obj.aux_phi_exp.rows(temp_ind) * sqrt_B_ys) / 2;
        sqrt_A_xs = arma::sqrt(n.rows(temp_ind) + b%b) - b;
    } else{
        sqrt_A_xs = n.elem(temp_ind) / arma::vectorise(trans_obj.aux_phi_exp.rows(temp_ind) * sqrt_B_ys);
    }
        
    arma::vec the_a_xs = arma::vectorise(sqrt_A_xs % sqrt_A_xs);
    //
    return the_a_xs;
}

template<>
arma::vec mfe<ntu>::marg_y_inv(const arma::mat& A_xs, arma::uvec* ys)
const
{
    arma::uvec temp_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::vec a_NTU = m.elem(temp_ind);
    arma::mat B_NTU = arma::trans( elem_prod(trans_obj.aux_alpha_exp.cols(temp_ind)/trans_obj.aux_gamma_exp.cols(temp_ind), A_xs) );
    arma::mat C_NTU = arma::trans(trans_obj.aux_gamma_exp.cols(temp_ind));

    arma::vec the_b_ys = inv_pwa(a_NTU, B_NTU, C_NTU, 1.0);
    //
    return the_b_ys;
}

template<>
arma::vec mfe<tu>::marg_y_inv(const arma::mat& A_xs, arma::uvec* ys)
const
{
    arma::uvec temp_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat sqrt_B_ys;
    arma::mat sqrt_A_xs = arma::sqrt(A_xs);

    if (!need_norm) {
        arma::mat b = arma::trans(sqrt_A_xs.t() * trans_obj.aux_phi_exp.cols(temp_ind)) / 2; // not sure about this
        sqrt_B_ys = arma::sqrt(m.rows(temp_ind) + b%b) - b;
    } else {
        sqrt_B_ys = m.elem(temp_ind) / arma::vectorise(arma::trans(sqrt_A_xs.t() * trans_obj.aux_phi_exp.cols(temp_ind))); // not sure about this
    }
        
    arma::vec the_b_ys = arma::vectorise(sqrt_B_ys % sqrt_B_ys);
    //
    return the_b_ys;
}

}

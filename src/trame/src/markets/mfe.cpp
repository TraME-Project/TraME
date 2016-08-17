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
 * Matching Function Equilibrium (MFE) market
 *
 * Keith O'Hara
 * 08/16/2016
 */

#include "trame.hpp"

void trame::mfe::build_ETU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, arma::mat tau_inp, double* sigma_inp, bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = true;
    }

    if (sigma_inp) {
        sigma = *sigma_inp;
    } else {
        sigma = 1;
    }

    trans_obj.build_ETU(alpha_inp,gamma_inp,tau_inp);

    arums_G.build(nbX,nbY,sigma,outsideOption);
    arums_H.build(nbY,nbX,sigma,outsideOption);

    mmf_obj.build_ETU(n,m,arma::exp(-alpha_inp/tau_inp),arma::exp(-gamma_inp/tau_inp),-1/tau_inp,need_norm);
}

void trame::mfe::build_LTU(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, double* sigma_inp, bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = true;
    }

    if (sigma_inp) {
        sigma = *sigma_inp;
    } else {
        sigma = 1;
    }

    trans_obj.build_LTU(lambda_inp,phi_inp);

    arums_G.build(nbX,nbY,sigma,outsideOption);
    arums_H.build(nbY,nbX,sigma,outsideOption);

    mmf_obj.build_LTU(n,m,lambda_inp,arma::exp(phi_inp/sigma),need_norm);
}

void trame::mfe::build_NTU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, double* sigma_inp, bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = true;
    }

    if (sigma_inp) {
        sigma = *sigma_inp;
    } else {
        sigma = 1;
    }

    trans_obj.build_NTU(alpha_inp,gamma_inp);

    arums_G.build(nbX,nbY,sigma,outsideOption);
    arums_H.build(nbY,nbX,sigma,outsideOption);

    mmf_obj.build_NTU(n,m,arma::exp(alpha_inp/sigma),arma::exp(gamma_inp/sigma),need_norm);
}

void trame::mfe::build_TU(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, double* sigma_inp, bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = true;
    }

    if (sigma_inp) {
        sigma = *sigma_inp;
    } else {
        sigma = 1;
    }

    trans_obj.build_TU(phi_inp);

    arums_G.build(nbX,nbY,sigma,outsideOption);
    arums_H.build(nbY,nbX,sigma,outsideOption);

    mmf_obj.build_TU(n,m,arma::exp(phi_inp/(2*sigma)),need_norm);
}

void trame::mfe::trans()
{
    arma::vec n_temp = n;
    n = m;
    m = n_temp;
    // Keith: fill in normalization later

    mmf_obj.trans();
    trans_obj.trans();

    logit arums_G_temp = arums_G;
    arums_G = arums_H;
    arums_H = arums_G_temp;
    //
}

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
 * affinity model class
 *
 * Keith O'Hara
 * 09/20/2016
 */

#include "trame.hpp"

void trame::affinity::build(const arma::mat& X_inp, const arma::mat& Y_inp)
{
    this->build(X_inp,Y_inp,NULL,NULL,NULL);
}

void trame::affinity::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    this->build(X_inp,Y_inp,&n_inp,&m_inp,NULL);
}

void trame::affinity::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec* n_inp, const arma::vec* m_inp, const double* sigma_inp)
{
    need_norm = false;

    nbX = X_inp.n_rows;
    nbY = Y_inp.n_rows;

    dX = X_inp.n_cols;
    dY = Y_inp.n_cols;

    nbParams = dX*dY;

    if (sigma_inp) {
        sigma = *sigma_inp;
    } else {
        sigma = 1.0;
    }
    //
    if (n_inp) {
        n = *n_inp;
    } else {
        n.ones(nbX,1);
    }
    if (m_inp) {
        m = *m_inp;
    } else {
        m.ones(nbY,1);
    }

    if (arma::accu(n) != arma::accu(m)) {
        printf("Unequal mass of individuals in an affinity model.\n");
    }
    //
}

arma::mat trame::affinity::Phi_xyk()
{
    return phi_xyk_aux;
}

arma::mat trame::affinity::Phi_xy(const arma::mat& lambda)
{
    arma::mat phi_xyk_temp = arma::resize(phi_xyk_aux,nbX*nbY,nbParams); 
    arma::mat ret = phi_xyk_temp * lambda;
    //
    return ret;
}

arma::mat trame::affinity::Phi_k(const arma::mat& mu_hat)
{
    arma::mat phi_xyk_temp = arma::resize(phi_xyk_aux,nbX*nbY,nbParams); 
    //arma::mat ret = arma::vectorise(arma::trans(arma::vectorise(mu_hat))*phi_xyk_temp);
    arma::mat ret = phi_xyk_temp.t() * arma::vectorise(mu_hat);
    //
    return ret;
}

void trame::affinity::dparam(arma::mat* dparams_inp, arma::mat& dparamsPsi_out, arma::mat* dparamsG_out, arma::mat* dparamsH_out)
{
    arma::mat dparams_mat;
    
    if (dparams_inp) {
        dparams_mat = *dparams_inp;
    } else {
        dparams_mat = arma::eye(nbParams,nbParams);
    }
    //
    dparamsPsi_out = Phi_xy(dparams_mat);
    //
    if (dparamsG_out) {
        *dparamsG_out = arma::zeros(0,dparams_mat.n_cols);
    }
    if (dparamsH_out) {
        *dparamsH_out = arma::zeros(0,dparams_mat.n_cols);
    }
}

trame::mfe<trame::mmf> trame::affinity::build_market(const arma::mat& theta)
{
    arma::mat Phi_mkt = Phi_xy(arma::vectorise(theta));
    Phi_mkt.resize(nbX,nbY);

    trame::mfe<trame::mmf> mkt_ret;
    mkt_ret.build_TU(n,m,Phi_mkt,&sigma,false);
    //
    return mkt_ret;
}

void trame::affinity::init_param(arma::mat& params)
{
    param.zeros(nbParams,1);
}

bool trame::affinity::mme_woregul(const arma::mat& mu_hat, double* xtol_ret, int* max_iter, double* tol_ipfp, double* max_iter_ipfp)
{
    arma::vec theta_0;
    init_param(theta_0);

    trame::mfe<trame::mmf> mkt_obj(theta_0);

    arma::mat C_hat = Phi_k(mu_hat);
    //
    double total_mass = arma::accu(n);

    arma::vec p = n / total_mass;
    arma::vec q = m / total_mass;

    arma::mat IX = arma::ones(nbX,1);
    arma::mat tIY = arma::ones(1,nbY);

    arma::mat f = p * tIY;
    arma::mat g = IX * q.t();

    arma::mat mu_hat / total_mass
    arma::mat v = arma::zeros(nbY,1);
    //
    
}

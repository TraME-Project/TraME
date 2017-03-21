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
 * affinity model class
 *
 * Keith O'Hara
 * 09/20/2016
 *
 * This version:
 * 11/27/2016
 */

#include "trame.hpp"

void trame::affinity::build(const arma::mat& X_inp, const arma::mat& Y_inp)
{
    this->build_int(X_inp,Y_inp,NULL,NULL,NULL);
}

void trame::affinity::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    this->build_int(X_inp,Y_inp,&n_inp,&m_inp,NULL);
}

void trame::affinity::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp, double sigma_inp)
{
    this->build_int(X_inp,Y_inp,&n_inp,&m_inp,&sigma_inp);
}

void trame::affinity::build_int(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec* n_inp, const arma::vec* m_inp, const double* sigma_inp)
{
    need_norm = false;

    nbX = X_inp.n_rows;
    nbY = Y_inp.n_rows;

    dX = X_inp.n_cols;
    dY = Y_inp.n_cols;

    nbParams = dX*dY;

    sigma = (sigma_inp) ? *sigma_inp : 1.0;
    //
    n = (n_inp) ? *n_inp : arma::ones(nbX,1);
    m = (m_inp) ? *m_inp : arma::ones(nbY,1);

    phi_xyk_aux = arma::kron(Y_inp,X_inp);
    //
    if (arma::accu(n) != arma::accu(m)) {
        printf("Unequal mass of individuals in an affinity model.\n");
    }
    //
}

trame::mfe<trame::transfers::tu> trame::affinity::build_market(const arma::mat& theta)
{
    arma::mat Phi_mkt = Phi_xy(arma::vectorise(theta));
    Phi_mkt.reshape(nbX,nbY);

    trame::mfe<trame::transfers::tu> mkt_ret;
    mkt_ret.build_TU(n,m,Phi_mkt,&sigma,false);
    //
    return mkt_ret;
}

void trame::affinity::dparam(arma::mat* dparams_inp, arma::mat& dparamsPsi_out, arma::mat* dparamsG_out, arma::mat* dparamsH_out)
{
    arma::mat dparams_mat = (dparams_inp) ? *dparams_inp : arma::eye(nbParams,nbParams);

    dparamsPsi_out = Phi_xy(dparams_mat);
    //
    if (dparamsG_out) {
        *dparamsG_out = arma::zeros(0,dparams_mat.n_cols);
    }
    if (dparamsH_out) {
        *dparamsH_out = arma::zeros(0,dparams_mat.n_cols);
    }
}

bool trame::affinity::mme_woregul(const arma::mat& mu_hat, arma::mat& theta_hat, double& val_ret, double* xtol_rel_inp, int* max_eval_inp, double* tol_ipfp_inp, double* max_iter_ipfp_inp)
{
    bool success = false;
    //
    double tol_ipfp = (tol_ipfp_inp) ? *tol_ipfp_inp : 1E-14;
    int max_iter_ipfp = (max_iter_ipfp_inp) ? *max_iter_ipfp_inp : 1E05;
    //
    arma::vec theta_0;
    init_param(theta_0);

    //trame::mfe<trame::mmf> mkt_obj = build_market(theta_0);

    arma::mat C_hat = Phi_k(mu_hat);
    //
    double total_mass = arma::accu(n);

    if (std::abs(arma::accu(mu_hat) - total_mass) > 1E-06) { // we use this instead of a raw '!=' to account for rounding error
        printf("Total number of couples does not coincide with margins.\n");
        return success;
    }
    //
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
    trame_mme_woregal_opt_data opt_data;

    opt_data.max_iter_ipfp = max_iter_ipfp;
    opt_data.tol_ipfp = tol_ipfp;

    opt_data.nbX = nbX;
    opt_data.nbY = nbY;

    opt_data.dX = dX;
    opt_data.dY = dY;

    opt_data.sigma = sigma;

    opt_data.p = p;
    opt_data.q = q;

    opt_data.IX = IX;
    opt_data.tIY = tIY;

    opt_data.f = f;
    opt_data.g = g;

    opt_data.v = v;
    opt_data.Pi_hat = Pi_hat;

    opt_data.phi_xy = phi_xy; // should be (nbX*nbY) x (nbParams)
    //
    arma::vec opt_vec = arma::zeros(dX*dY,1);

    success = generic_optim(opt_vec,trame::affinity::mme_woregul_opt_objfn,&opt_data);
    //
    if (success) {
        theta_hat = opt_vec;
        val_ret = trame::affinity::mme_woregul_opt_objfn(opt_vec,NULL,&opt_data);
    }
    //
    return success;
}

bool trame::affinity::mme_regul(const arma::mat& mu_hat, const double& lambda, arma::mat& theta_hat, double& val_ret, double* xtol_rel_inp, int* max_eval_inp, double* tol_ipfp_inp, double* max_iter_ipfp_inp)
{
    bool success = false;
    //
    double xtol_rel = (xtol_rel_inp) ? *xtol_rel_inp : 1E-04;
    int max_eval = (max_eval_inp) ? *max_eval_inp : 1E05;

    double tol_ipfp = (tol_ipfp_inp) ? *tol_ipfp_inp : 1E-14;
    int max_iter_ipfp = (max_iter_ipfp_inp) ? *max_iter_ipfp_inp : 1E05;
    //
    arma::vec theta_0;
    init_param(theta_0);

    //trame::mfe<trame::mmf> mkt_obj = build_market(theta_0);

    arma::mat C_hat = Phi_k(mu_hat);
    //
    double total_mass = arma::accu(n);

    if (std::abs(arma::accu(mu_hat) - total_mass) > 1E-06) { // we use this instead of a raw '!=' to account for rounding error
        printf("Total number of couples does not coincide with margins.\n");
        return success;
    }
    //
    arma::vec p = n / total_mass;
    arma::vec q = m / total_mass;

    arma::mat IX = arma::ones(nbX,1);
    arma::mat tIY = arma::ones(1,nbY);

    arma::mat f = p * tIY;
    arma::mat g = IX * q.t();

    arma::mat Pi_hat = mu_hat / total_mass;
    arma::mat v = arma::zeros(1,nbY);
    arma::mat A = arma::zeros(dX*dY,1);

    arma::mat Phi = arma::reshape(Phi_xy(arma::vectorise(A)),nbX,nbY);
    //
    int iter_ipfp = 0, iter_count = 0;
    double err_ipfp = 2*tol_ipfp, err_val = 1.0;
    double t_k = 0.3; // step size for the prox grad algorithm (or grad descent when lambda=0)
    double alpha = 1.0; // for optimality check
    
    double the_val = 1.0, the_val_old = 1E04;
    arma::vec d, d_opt;
    arma::mat v_next = v, u, Pi, the_grad, A_mat, U, V, D, svd_mat, D_opt, opt_mat;

    while (err_val > xtol_rel && iter_count < max_eval) {
        iter_count++;

        Phi = arma::reshape(Phi_xy(arma::vectorise(A)),nbX,nbY);
        err_ipfp= 2*tol_ipfp;
        iter_ipfp = 0;

        while (err_ipfp > tol_ipfp && iter_ipfp < max_iter_ipfp) {
            iter_ipfp++;

            u = sigma * arma::log(arma::sum(g % arma::exp((Phi - IX * v)/sigma),1));
            v_next = sigma * arma::log(arma::sum(f % arma::exp((Phi - u * tIY)/sigma),0));

            err_ipfp = elem_max(arma::abs( arma::sum(g % arma::exp((Phi - IX * v_next - u*tIY)/sigma),1) - 1.0 ));
            v = v_next;
        }
        //
        Pi = f % g % arma::exp( (Phi - IX*v - u*tIY)/sigma );
        the_grad = Phi_k(Pi - Pi_hat);

        A -= t_k*the_grad;
        //
        if (lambda > 0) {
            // compute the proximal operator
            A_mat = arma::reshape(A,dX,dY);
            arma::svd(U,d,V,A_mat);
            D = arma::diagmat(elem_max(d - lambda*t_k,0.0));
            A = arma::vectorise(U * D * V.t());
        } // if lambda = 0 then we are just taking one step of gradient descent
        //
        if (iter_count % 10 == 0) {
            //alpha = 1.0;
            svd_mat = arma::reshape(A - alpha*the_grad,dX,dY);
            arma::svd(U,d_opt,V,svd_mat);
            D_opt = arma::diagmat(elem_max(d_opt - alpha*lambda, 0.0));

            opt_mat = arma::accu(arma::pow(A - arma::vectorise(U * D_opt *V.t()), 2));
        }
        //
        if (lambda > 0) {
            the_val = arma::accu(the_grad % arma::vectorise(A)) - sigma * arma::accu(Pi % arma::log(Pi)) + lambda * arma::accu(D);
        } else {
            the_val = arma::accu(the_grad % arma::vectorise(A)) - sigma * arma::accu(Pi % arma::log(Pi));
        }

        err_val = std::abs(the_val - the_val_old);
        //
        the_val_old = the_val;
    }

    if (err_val <= xtol_rel && iter_count < max_eval) {
        success = true;
    }
    //
    theta_hat = arma::vectorise(A);
    val_ret = the_val;
    //
    return success;
}

// internal

arma::mat trame::affinity::Phi_xyk()
{
    return phi_xyk_aux;
}

arma::mat trame::affinity::Phi_xy(const arma::mat& lambda)
{
    arma::mat phi_xyk_temp = arma::reshape(phi_xyk_aux,nbX*nbY,nbParams); 
    arma::mat ret = phi_xyk_temp * lambda;
    //
    return ret;
}

arma::mat trame::affinity::Phi_k(const arma::mat& mu_hat)
{
    arma::mat phi_xyk_temp = arma::reshape(phi_xyk_aux,nbX*nbY,nbParams); 
    arma::mat ret = arma::vectorise(arma::trans(arma::vectorise(mu_hat))*phi_xyk_temp);
    //arma::mat ret = phi_xyk_temp.t() * arma::vectorise(mu_hat);
    //
    return ret;
}

void trame::affinity::init_param(arma::mat& params)
{
    params.zeros(nbParams,1);
}

/*
 * optimization-related functions
 */

double trame::affinity::mme_woregul_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void* opt_data)
{
    trame_mme_woregal_opt_data *d = reinterpret_cast<trame_mme_woregal_opt_data*>(opt_data);
    //
    int nbX = d->nbX;
    int nbY = d->nbY;

    //int dX = d->dX;
    //int dY = d->dY;
    
    int max_iter_ipfp = d->max_iter_ipfp;
    double tol_ipfp = d->tol_ipfp;

    double sigma = d->sigma;

    arma::vec p = d->p;
    arma::vec q = d->q;

    arma::mat IX = d->IX;
    arma::mat tIY = d->tIY;

    arma::mat f = d->f;
    arma::mat g = d->g;

    arma::mat v = d->v;
    arma::mat Pi_hat = d->Pi_hat;

    arma::mat phi_xy = d->phi_xy; // should be (nbX*nbY) x (nbParams), replaces a member function call
    //
    arma::mat Phi = arma::reshape(phi_xy * vals_inp,nbX,nbY);
    //
    int iter_ipfp = 0;
    double err_ipfp= 2*tol_ipfp;
    arma::mat v_next = v, u;

    while (err_ipfp > tol_ipfp && iter_ipfp < max_iter_ipfp) {
        iter_ipfp++;

        u = sigma * arma::log(arma::sum(g % arma::exp((Phi - IX * v)/sigma),1));
        v_next = sigma * arma::log(arma::sum( f % arma::exp((Phi - u * tIY)/sigma),0));

        err_ipfp = elem_max(arma::abs( arma::sum(g % arma::exp((Phi - IX * v_next - u*tIY)/sigma),1) - 1.0 ));
        v = v_next;
    }
    //
    arma::mat Pi = f % g % arma::exp( (Phi - IX*v - u*tIY)/sigma );
    arma::mat the_grad = phi_xy.t() * arma::vectorise(Pi - Pi_hat);

    if (grad) {
        *grad = the_grad;
    }
    //
    // update v for the next opt call
    d->v = v;
    opt_data = reinterpret_cast<void*>(d);
    //
    double ret = arma::accu(the_grad % arma::vectorise(vals_inp)) - sigma*arma::accu(Pi%arma::log(Pi));
    //
    return ret;
}

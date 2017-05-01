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
 * nodal_newton for logit MFE case
 *
 * Keith O'Hara
 * 01/17/2016
 *
 * This version:
 * 02/15/2017
 */

// internal nodal_newton

template<typename Tt>
bool 
nodal_newton_int(const mfe<Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, double* val_out, const double* tol_inp, const int* max_iter_inp)
{
    bool success = false;
    //
    int nbX = market.nbX;
    int nbY = market.nbY;

    double sigma = market.sigma; // we don't check for arums_G.sigma == arums_H.sigma

    trame_mfe_opt_data<Tt> opt_data;
    opt_data.market = market;

    arma::vec sol_vec = -sigma*arma::join_cols(arma::log(market.n/2.0),arma::log(market.m/2.0)); // initial guess
    
    success = nodal_newton_optim(sol_vec,nodal_newton_opt_objfn<Tt>,&opt_data,nodal_newton_jacobian<Tt>,&opt_data);
    //
    // construct equilibrium objects
    arma::vec us = sol_vec.rows(0,nbX-1);
    arma::vec vs = sol_vec.rows(nbX,nbX+nbY-1);

    arma::vec mu_x0_s = arma::exp(-us/sigma);
    arma::vec mu_0y_s = arma::exp(-vs/sigma);

    arma::mat mu = market.mmf_obj.M(mu_x0_s,mu_0y_s);
    //
    // return equilibrium objects
    if (mu_out) {
        *mu_out = mu;
    }

    if (mu_x0_out) {
        *mu_x0_out = mu_x0_s;
    }
    if (mu_0y_out) {
        *mu_0y_out = mu_0y_s;
    }

    if (U_out) {
        *U_out = sigma * arma::log(elem_div(mu,mu_x0_s));
    }
    if (V_out) {
        *V_out = sigma * arma::trans(arma::log(elem_div(mu.t(),mu_0y_s)));
    }

    /*
    if (val_out) {
        *val_out = val;
    }
    */
    //
    return success;
}

// wrappers 

template<typename Tt>
bool 
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out)
{
    bool res = nodal_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    
    return res;
}

template<typename Tt>
bool 
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out, const double& tol_inp)
{
    bool res = nodal_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,NULL);
    
    return res;
}

template<typename Tt>
bool 
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out, const int& max_iter_inp)
{
    bool res = nodal_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&max_iter_inp);
    
    return res;
}

template<typename Tt>
bool 
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp)
{
    bool res = nodal_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,&max_iter_inp);
    
    return res;
}

template<typename Tt>
bool 
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    bool res = nodal_newton_int(market,&mu_out,NULL,NULL,&U_out,&V_out,NULL,NULL,NULL);
    
    return res;
}

template<typename Tt>
bool 
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, double& val_out, const double* tol_inp, const int* max_iter_inp)
{
    bool res = nodal_newton_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,&val_out,tol_inp,max_iter_inp);
    
    return res;
}

// optimization function

template<typename Tt>
arma::vec 
nodal_newton_opt_objfn(const arma::vec& vals_inp, void *opt_data)
{
    trame_mfe_opt_data<Tt> *d = reinterpret_cast<trame_mfe_opt_data<Tt>*>(opt_data);
    //
    int nbX = d->market.nbX;
    int nbY = d->market.nbY;
    double sigma = d->market.sigma;

    arma::vec us = vals_inp.rows(0,nbX-1);
    arma::vec vs = vals_inp.rows(nbX,nbX+nbY-1);

    arma::vec mu_x0_s = arma::exp(-us/sigma);
    arma::vec mu_0y_s = arma::exp(-vs/sigma);

    arma::mat mu = d->market.mmf_obj.M(mu_x0_s,mu_0y_s);
    //
    arma::vec ret = arma::join_cols(mu_x0_s + arma::sum(mu,1) - d->market.n, mu_0y_s + arma::trans(arma::sum(mu,0)) - d->market.m);
    //
    return ret;
}

template<typename Tt>
arma::mat 
nodal_newton_jacobian(const arma::vec& vals_inp, void *jacob_data)
{
    trame_mfe_opt_data<Tt> *d = reinterpret_cast<trame_mfe_opt_data<Tt>*>(jacob_data);
    //
    int nbX = d->market.nbX;
    int nbY = d->market.nbY;
    double sigma = d->market.sigma;

    arma::vec us = vals_inp.rows(0,nbX-1);
    arma::vec vs = vals_inp.rows(nbX,nbX+nbY-1);

    arma::vec mu_x0_s = arma::exp(-us/sigma);
    arma::vec mu_0y_s = arma::exp(-vs/sigma);

    arma::mat mu = d->market.mmf_obj.M(mu_x0_s,mu_0y_s);
    //
    arma::mat du_s = d->market.mmfs_obj.du_Psi(us,vs);
    arma::mat dv_s = 1.0 - du_s;

    arma::mat Delta_11 = arma::diagmat(mu_x0_s + arma::sum(mu%du_s,1));
    arma::mat Delta_12 = mu%dv_s;

    arma::mat Delta_21 = arma::trans(mu%du_s);
    arma::mat Delta_22 = arma::diagmat(mu_0y_s + arma::trans(arma::sum(mu%dv_s,0)));

    arma::mat Delta = arma::join_cols( arma::join_rows(Delta_11,Delta_12), arma::join_rows(Delta_21,Delta_22) );

    arma::mat ret = - Delta / sigma;
    //
    return ret;
}

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
 * nodal_newton for MFE case
 *
 * Keith O'Hara
 * 01/17/2016
 *
 * This version:
 * 07/26/2017
 */

// internal nodal_newton

template<typename Tt>
bool
nodal_newton_int(const mfe<Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, double* val_out, const double* err_tol_inp, const int* max_iter_inp)
{
    bool success = false;
    //
    const double err_tol = (err_tol_inp) ? *err_tol_inp : 1E-06;
    const int max_iter = (max_iter_inp) ? *max_iter_inp : 2000;

    const int nbX = market.nbX;
    const int nbY = market.nbY;

    const double sigma = market.sigma; // we don't check for arums_G.sigma == arums_H.sigma

    trame_mfe_opt_data<Tt> opt_data;
    opt_data.market = market;

    arma::vec sol_vec = -sigma*arma::join_cols(arma::log(market.n/2.0),arma::log(market.m/2.0)); // initial guess
    
    success = nodal_newton_optim(sol_vec,nodal_newton_opt_objfn<Tt>,&opt_data,nodal_newton_jacobian<Tt>,&opt_data,nullptr,&err_tol,&max_iter);

    //
    // construct equilibrium objects

    const arma::vec us = sol_vec.rows(0,nbX-1);
    const arma::vec vs = sol_vec.rows(nbX,nbX+nbY-1);

    const arma::vec mu_x0_s = arma::exp(-us/sigma);
    const arma::vec mu_0y_s = arma::exp(-vs/sigma);

    const arma::mat mu = market.mmfs_obj.M(mu_x0_s,mu_0y_s);

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
    return nodal_newton_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tt>
bool
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out, const double err_tol_inp)
{
    return nodal_newton_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,&err_tol_inp,nullptr);
}

template<typename Tt>
bool
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out, const int max_iter_inp)
{
    return nodal_newton_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,&max_iter_inp);
}

template<typename Tt>
bool
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out, const double err_tol_inp, const int max_iter_inp)
{
    return nodal_newton_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,&err_tol_inp,&max_iter_inp);
}

template<typename Tt>
bool
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    return nodal_newton_int(market,&mu_out,nullptr,nullptr,&U_out,&V_out,nullptr,nullptr,nullptr);
}

template<typename Tt>
bool
nodal_newton(const mfe<Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, double& val_out, const double* err_tol_inp, const int* max_iter_inp)
{
    return nodal_newton_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,&val_out,err_tol_inp,max_iter_inp);
}

//
// optimization functions

inline
bool
nodal_newton_optim(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data, arma::vec* value_out, const double* err_tol_inp, const int* max_iter_inp)
{
    optim::opt_settings opt_params;

    if (err_tol_inp) {
        opt_params.err_tol = *err_tol_inp;
    }

    if (max_iter_inp) {
        opt_params.iter_max = *max_iter_inp;
    }

    return optim::broyden_int(init_out_vals,opt_objfn,opt_data,&opt_params);
}

inline
bool 
nodal_newton_optim(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data,
                   std::function<arma::mat (const arma::vec& vals_inp, void* jacob_data)> jacob_objfn, void* jacob_data, arma::vec* value_out, const double* err_tol_inp, const int* max_iter_inp)
{
    optim::opt_settings opt_params;

    if (err_tol_inp) {
        opt_params.err_tol = *err_tol_inp;
    }

    if (max_iter_inp) {
        opt_params.iter_max = *max_iter_inp;
    }

    return optim::broyden_df_int(init_out_vals,opt_objfn,opt_data,jacob_objfn,jacob_data,&opt_params);
}

template<typename Tt>
arma::vec
nodal_newton_opt_objfn(const arma::vec& vals_inp, void *opt_data)
{
    trame_mfe_opt_data<Tt> *d = reinterpret_cast<trame_mfe_opt_data<Tt>*>(opt_data);
    //
    const int nbX = d->market.nbX;
    const int nbY = d->market.nbY;
    const double sigma = d->market.sigma;

    const arma::vec us = vals_inp.rows(0,nbX-1);
    const arma::vec vs = vals_inp.rows(nbX,nbX+nbY-1);

    const arma::vec mu_x0_s = arma::exp(-us/sigma);
    const arma::vec mu_0y_s = arma::exp(-vs/sigma);

    const arma::mat mu = d->market.mmfs_obj.M(mu_x0_s,mu_0y_s);
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
    const int nbX = d->market.nbX;
    const int nbY = d->market.nbY;
    const double sigma = d->market.sigma;

    const arma::vec us = vals_inp.rows(0,nbX-1);
    const arma::vec vs = vals_inp.rows(nbX,nbX+nbY-1);

    const arma::vec mu_x0_s = arma::exp(-us/sigma);
    const arma::vec mu_0y_s = arma::exp(-vs/sigma);

    const arma::mat mu = d->market.mmfs_obj.M(mu_x0_s,mu_0y_s);
    //
    // arma::mat du_s = d->market.mmfs_obj.du_Psi(us,vs);
    const arma::mat du_s = d->market.mmfs_obj.dmu_x0(mu_x0_s,mu_0y_s);
    const arma::mat dv_s = d->market.mmfs_obj.dmu_0y(mu_x0_s,mu_0y_s);

    const arma::mat Delta_11 = - arma::diagmat(mu_x0_s % (1.0 + arma::sum(du_s,1)));
    const arma::mat Delta_22 = - arma::diagmat(mu_0y_s % (1.0 + arma::trans(arma::sum(dv_s,0))));

    const arma::mat Delta_12 = - arma::trans(elem_prod(mu_0y_s,dv_s.t()));
    const arma::mat Delta_21 = - arma::trans(elem_prod(mu_x0_s,du_s));

    arma::mat ret = arma::join_cols( arma::join_rows(Delta_11,Delta_12), arma::join_rows(Delta_21,Delta_22) );
    //
    return ret;
}

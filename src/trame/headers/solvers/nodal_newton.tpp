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
 * nodal_newton
 *
 * Keith O'Hara
 * 01/17/2016
 */

// internal nodal_newton

template<typename Tm>
bool nodal_newton_int(const mfe<Tm>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, double* val_out, const double* tol_inp, const int* max_iter_inp)
{
    bool success = false;
    //
    int nbX = market.nbX;
    int nbY = market.nbY;
    double sigma = market.sigma;

    trame_mfe_opt_data<Tm> opt_data;
    opt_data.market = market;

    arma::vec sol_vec = -sigma*arma::join_cols(arma::log(market.n/2.0),arma::log(market.m/2.0));
    
    success = nodal_newton_optim(sol_vec,nodal_newton_opt_objfn<Tm>,&opt_data);
    std::cout << "nodal_newton_optim done" << std::endl;
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

    /*if (val_out) {
        *val_out = val;
    }*/
    //
    return success;
}

// wrappers 

template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out)
{
    bool res = nodal_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    
    return res;
}

template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out, const double& tol_inp)
{
    bool res = nodal_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,NULL);
    
    return res;
}

template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out, const int& max_iter_inp)
{
    bool res = nodal_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&max_iter_inp);
    
    return res;
}

template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp)
{
    bool res = nodal_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,&max_iter_inp);
    
    return res;
}

template<typename Tm>
bool nodal_newton(const mfe<Tm>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, double& val_out, const double* tol_inp, const int* max_iter_inp)
{
    bool res = nodal_newton_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,&val_out,tol_inp,max_iter_inp);
    
    return res;
}

// optimization function

template<typename Tm>
arma::vec nodal_newton_opt_objfn(const arma::vec& vals_inp, void *opt_data)
{
    trame_mfe_opt_data<Tm> *d = reinterpret_cast<trame_mfe_opt_data<Tm>*>(opt_data);
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

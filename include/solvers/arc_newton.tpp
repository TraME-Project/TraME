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
 * Arc Newton for TU_general
 *
 * Keith O'Hara
 * 01/17/2016
 *
 * This version:
 * 02/15/2017
 */

// internal arc_newton

template<typename Ta, typename Tm>
bool arc_newton_int(const dse<Ta,Tm>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, double* val_out, const double* tol_inp, const int* max_iter_inp)
{
    bool success = false;
    //
    int nbX = market.nbX;
    int nbY = market.nbY;

    trame_market_opt_data<Ta,Tm> opt_data;
    opt_data.market = market;

    arma::vec sol_vec = arma::vectorise(w_upper_bound(market));
    
    //success = arc_newton_optim(sol_vec,arc_newton_opt_objfn<Ta>,&opt_data);
    success = arc_newton_optim(sol_vec,arc_newton_opt_objfn<Ta,Tm>,&opt_data,arc_newton_jacobian<Ta,Tm>,&opt_data);
    //
    // construct equilibrium objects
    arma::mat sol_mat = arma::reshape(sol_vec,nbX,nbY);
    arma::mat U = market.trans_obj.UW(sol_mat);

    Ta* arums_G = const_cast<Ta*>(&market.arums_G); // Keith: this recast is unsafe, change later
    
    arma::mat mu_G;
    arums_G->G(market.n,U,mu_G);
    //
    // return equilibrium objects
    if (mu_out) {
        *mu_out = mu_G;
    }

    if (mu_x0_out) {
        *mu_x0_out = market.n - arma::sum(mu_G,1);
    }
    if (mu_0y_out) {
        *mu_0y_out = market.m - arma::trans(arma::sum(mu_G,0));
    }

    if (U_out) {
        *U_out = U;
    }
    if (V_out) {
        *V_out = market.trans_obj.VW(sol_vec);
    }

    /*if (val_out) {
        *val_out = val;
    }*/
    //
    return success;
}

// wrappers 

template<typename Ta, typename Tm>
bool arc_newton(const dse<Ta,Tm>& market, arma::mat& mu_out)
{
    bool res = arc_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    
    return res;
}

template<typename Ta, typename Tm>
bool arc_newton(const dse<Ta,Tm>& market, arma::mat& mu_out, const double& tol_inp)
{
    bool res = arc_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,NULL);
    
    return res;
}

template<typename Ta, typename Tm>
bool arc_newton(const dse<Ta,Tm>& market, arma::mat& mu_out, const int& max_iter_inp)
{
    bool res = arc_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&max_iter_inp);
    
    return res;
}

template<typename Ta, typename Tm>
bool arc_newton(const dse<Ta,Tm>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp)
{
    bool res = arc_newton_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,&max_iter_inp);
    
    return res;
}

template<typename Ta, typename Tm>
bool arc_newton(const dse<Ta,Tm>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    bool res = arc_newton_int(market,&mu_out,NULL,NULL,&U_out,&V_out,NULL,NULL,NULL);
    
    return res;
}

template<typename Ta, typename Tm>
bool arc_newton(const dse<Ta,Tm>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, double& val_out, const double* tol_inp, const int* max_iter_inp)
{
    bool res = arc_newton_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,&val_out,tol_inp,max_iter_inp);
    
    return res;
}

// optimization function

template<typename Ta, typename Tm>
arma::vec arc_newton_opt_objfn(const arma::vec& vals_inp, void *opt_data)
{
    trame_market_opt_data<Ta,Tm> *d = reinterpret_cast<trame_market_opt_data<Ta,Tm>*>(opt_data);
    //
    int nbX = d->market.nbX;
    int nbY = d->market.nbY;

    arma::mat inp_mat = arma::reshape(vals_inp,nbX,nbY);

    arma::mat U = d->market.trans_obj.UW(inp_mat);
    arma::mat V = d->market.trans_obj.VW(inp_mat);

    arma::mat mu_G, mu_H;
    d->market.arums_G.G(d->market.n,U,mu_G);
    d->market.arums_H.G(d->market.m,V.t(),mu_H);
    //
    arma::vec ret = arma::vectorise(mu_G - mu_H.t());
    //
    return ret;
}

template<typename Ta, typename Tm>
arma::mat arc_newton_jacobian(const arma::vec& vals_inp, void *jacob_data)
{
    trame_market_opt_data<Ta,Tm> *d = reinterpret_cast<trame_market_opt_data<Ta,Tm>*>(jacob_data);
    //
    int nbX = d->market.nbX;
    int nbY = d->market.nbY;

    arma::mat inp_mat = arma::reshape(vals_inp,nbX,nbY);

    arma::mat U = d->market.trans_obj.UW(inp_mat);
    arma::mat V = d->market.trans_obj.VW(inp_mat);

    arma::mat dwUW = d->market.trans_obj.dw_UW(inp_mat);
    arma::mat dwVW = d->market.trans_obj.dw_VW(inp_mat);

    arma::mat D2G_mat, D2H_mat;
    d->market.arums_G.D2G(D2G_mat,d->market.n,U,true);
    d->market.arums_H.D2G(D2H_mat,d->market.m,V.t(),false);

    arma::mat term_1 = elem_prod(arma::vectorise(dwUW),D2G_mat);
    arma::mat term_2 = elem_prod(arma::vectorise(dwVW),D2H_mat);

    arma::mat ret = term_1 - term_2;
    //
    return ret;
}

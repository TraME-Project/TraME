/*################################################################################
  ##
  ##   Copyright (C) 2015-2018 the TraME Team:
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
 * Arc Newton for TU_general
 *
 * Keith O'Hara
 * 01/17/2016
 *
 * This version:
 * 02/04/2018
 */

//
// internal arc_newton

template<typename Tg, typename Th, typename Tt>
bool
arc_newton_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, 
               double* val_out, const double err_tol, const uint_t max_iter)
{
    bool success = false;

    //

    trame_market_opt_data<Tg,Th,Tt> opt_data;
    opt_data.market = market;

    // optim

    optim::algo_settings_t settings;

    settings.err_tol = err_tol;
    settings.iter_max = max_iter;

    arma::vec sol_vec = arma::vectorise(w_upper_bound(market)); // initial guess

    success = arc_newton_optim(sol_vec,arc_newton_opt_objfn<Tg,Th,Tt>,&opt_data,arc_newton_jacobian<Tg,Th,Tt>,&opt_data,&settings);

    //
    // construct equilibrium objects

    arma::mat sol_mat = arma::reshape(sol_vec,market.nbX,market.nbY);
    arma::mat U = market.transfers_obj.UW(sol_mat);

    arma::mat mu_G;
    market.arums_G.G(market.n,U,mu_G);

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
        *V_out = market.transfers_obj.VW(sol_mat);
    }

    if (val_out) {
        *val_out = arma::accu(settings.zero_values);
    }

    //

    return success;
}

//
// wrappers

template<typename Tg, typename Th, typename Tt>
bool
arc_newton(const dse<Tg,Th,Tt>& market, arma::mat& mu_out)
{
    return arc_newton_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
arc_newton(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double err_tol_inp, const uint_t max_iter_inp)
{
    return arc_newton_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,err_tol_inp,max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
arc_newton(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    return arc_newton_int(market,&mu_out,nullptr,nullptr,&U_out,&V_out,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
arc_newton(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, 
           double& val_out, const double err_tol_inp, const uint_t max_iter_inp)
{
    return arc_newton_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,&val_out,err_tol_inp,max_iter_inp);
}

//
// optimization functions

inline
bool
arc_newton_optim(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data, optim::algo_settings_t* settings_inp)
{
    return optim::broyden_df_int(init_out_vals,opt_objfn,opt_data,settings_inp);
}

inline
bool
arc_newton_optim(arma::vec& init_out_vals, std::function<arma::vec (const arma::vec& vals_inp, void* opt_data)> opt_objfn, void* opt_data,
                             std::function<arma::mat (const arma::vec& vals_inp, void* jacob_data)> jacob_objfn, void* jacob_data, optim::algo_settings_t* settings_inp)
{
    return optim::broyden_df_int(init_out_vals,opt_objfn,opt_data,jacob_objfn,jacob_data,settings_inp);
}

template<typename Tg, typename Th, typename Tt>
arma::vec
arc_newton_opt_objfn(const arma::vec& vals_inp, void *opt_data)
{
    trame_market_opt_data<Tg,Th,Tt> *d = reinterpret_cast<trame_market_opt_data<Tg,Th,Tt>*>(opt_data);
    
    //

    const uint_t nbX = d->market.nbX;
    const uint_t nbY = d->market.nbY;

    // arma::mat inp_mat = arma::reshape(vals_inp,nbX,nbY);
    arma::mat inp_mat(const_cast<double*>(vals_inp.memptr()),nbX,nbY,false,true); // this is potentially very unsafe, but more efficient?

    arma::mat U = d->market.transfers_obj.UW(inp_mat);
    arma::mat V = d->market.transfers_obj.VW(inp_mat);

    arma::mat mu_G, mu_H;
    d->market.arums_G.G(d->market.n,U,mu_G);
    d->market.arums_H.G(d->market.m,V.t(),mu_H);
    
    //

    return arma::vectorise(mu_G - mu_H.t());
}

template<typename Tg, typename Th, typename Tt>
arma::mat
arc_newton_jacobian(const arma::vec& vals_inp, void *jacob_data)
{
    trame_market_opt_data<Tg,Th,Tt> *d = reinterpret_cast<trame_market_opt_data<Tg,Th,Tt>*>(jacob_data);

    //
    
    const uint_t nbX = d->market.nbX;
    const uint_t nbY = d->market.nbY;

    // arma::mat inp_mat = arma::reshape(vals_inp,nbX,nbY);
    arma::mat inp_mat(const_cast<double*>(vals_inp.memptr()),nbX,nbY,false,true); // this is potentially unsafe, but more efficient?

    arma::mat U = d->market.transfers_obj.UW(inp_mat);
    arma::mat V = d->market.transfers_obj.VW(inp_mat);

    arma::mat dwUW = d->market.transfers_obj.dw_UW(inp_mat);
    arma::mat dwVW = d->market.transfers_obj.dw_VW(inp_mat);

    arma::mat D2G_mat, D2H_mat;
    d->market.arums_G.D2G(D2G_mat,d->market.n,U,true);
    d->market.arums_H.D2G(D2H_mat,d->market.m,V.t(),false);

    // arma::mat term_1 = elem_prod(arma::vectorise(dwUW),D2G_mat);
    // arma::mat term_2 = elem_prod(arma::vectorise(dwVW),D2H_mat);

    //

    return elem_prod(arma::vectorise(dwUW),D2G_mat) - elem_prod(arma::vectorise(dwVW),D2H_mat);
}

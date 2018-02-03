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
 * max_welfare for general ARUMs with TU
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 07/26/2017
 */

// internal max_welfare

template<typename Tg, typename Th, typename Tt>
bool
max_welfare_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, double* val_out, const double* err_tol_inp, const int* max_iter_inp)
{
    printf("max_welfare only works for TU transfers.\n");
    return false;
}

template<typename Tg, typename Th>
bool
max_welfare_int(const dse<Tg,Th,transfers::tu>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, double* val_out, const double* err_tol_inp, const int* max_iter_inp)
{
    bool success = false;
    
    const double err_tol = (err_tol_inp) ? *err_tol_inp : 1E-06;
    const int max_iter = (max_iter_inp) ? *max_iter_inp : 2000;

    //

    const int nbX = market.nbX;
    const int nbY = market.nbY;

    trame_market_opt_data<Tg,Tg,transfers::tu> opt_data;
    opt_data.market = market;

    // optim

    const int optim_method = 2;

    optim::algo_settings settings;

    settings.err_tol = err_tol;
    settings.iter_max = max_iter;

    arma::vec sol_vec = arma::vectorise(market.trans_obj.phi / 2.0); // initial value

    success = max_welfare_optim(sol_vec,max_welfare_opt_objfn<Tg,Th,transfers::tu>,&opt_data,&settings,optim_method);

    //
    // construct equilibrium objects

    // arma::mat U = arma::reshape(sol_vec,nbX,nbY);
    arma::mat U(sol_vec.memptr(),nbX,nbY,false,true);

    arma::mat mu_G, mu_H;
    double val_G = market.arums_G.G(market.n,U,mu_G);
    double val_H = market.arums_H.G(market.m,arma::trans(market.trans_obj.phi - U),mu_H);

    double val = val_G + val_H;
    
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
        *V_out = market.trans_obj.phi - U;
    }

    if (val_out) {
        *val_out = val;
    }
    //
    return success;
}

// wrappers

template<typename Tg, typename Th, typename Tt>
bool
max_welfare(const dse<Tg,Th,Tt>& market, arma::mat& mu_out)
{
    return max_welfare_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
max_welfare(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double err_tol_inp)
{
    return max_welfare_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,&err_tol_inp,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
max_welfare(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const int max_iter_inp)
{
    return max_welfare_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,&max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
max_welfare(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double err_tol_inp, const int max_iter_inp)
{
    return max_welfare_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,&err_tol_inp,&max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
max_welfare(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    return max_welfare_int(market,&mu_out,nullptr,nullptr,&U_out,&V_out,nullptr,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
max_welfare(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, double& val_out, const double* err_tol_inp, const int* max_iter_inp)
{
    return max_welfare_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,&val_out,err_tol_inp,max_iter_inp);
}

//
// optimization functions

inline
bool
max_welfare_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, optim::algo_settings* settings_inp, const int optim_method)
{
    if (optim_method == 1) {
        return optim::lbfgs_int(init_out_vals,opt_objfn,opt_data,settings_inp);
    } else if (optim_method == 2) {
        return optim::bfgs_int(init_out_vals,opt_objfn,opt_data,settings_inp);
    } else {
        printf("error: unrecognized optim_method choice.\n");
        return false;
    }
}

template<typename Tg, typename Th, typename Tt>
double
max_welfare_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void *opt_data)
{
    trame_market_opt_data<Tg,Th,Tt> *d = reinterpret_cast<trame_market_opt_data<Tg,Th,Tt>*>(opt_data);
    //
    const int nbX = d->market.nbX;
    const int nbY = d->market.nbY;

    arma::mat U = arma::reshape(vals_inp,nbX,nbY);

    arma::mat mu_G, mu_H;
    const double val_G = d->market.arums_G.G(d->market.n,U,mu_G);
    const double val_H = d->market.arums_H.G(d->market.m,arma::trans(d->market.trans_obj.phi - U),mu_H);
    //
    if (grad) {
        *grad = arma::vectorise(mu_G - mu_H.t());
    }
    //
    return val_G + val_H;
}

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
 * max_welfare
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/21/2017
 */

// internal max_welfare

template<typename Ta>
bool max_welfare_int(const dse<Ta>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, double* val_out, const double* tol_inp, const int* max_iter_inp)
{
    bool success = false;
    //
    // warnings
    if (!market.trans_obj.TU) {
        printf("max_welfare only works for TU transfers.\n");
        return false;
    }
    //
    double tol = (tol_inp) ? *tol_inp : 1E-06;
    int max_iter = (max_iter_inp) ? *max_iter_inp : 2000;

    int nbX = market.nbX;
    int nbY = market.nbY;

    trame_market_opt_data<Ta> opt_data;
    opt_data.market = market;

    double obj_val = 0;
    arma::vec sol_vec = arma::vectorise(market.trans_obj.phi / 2.0);
    
    success = max_welfare_optim(sol_vec,max_welfare_opt_objfn_2<Ta>,&opt_data,&obj_val,&tol,&max_iter);
    //
    // construct equilibrium objects
    arma::mat U = arma::reshape(sol_vec,nbX,nbY);

    Ta* arums_G = const_cast<Ta*>(&market.arums_G); // Keith: this recast is unsafe, change later
    Ta* arums_H = const_cast<Ta*>(&market.arums_H);
    
    arma::mat mu_G, mu_H;
    double val_G = arums_G->G(market.n,U,mu_G);
    double val_H = arums_H->G(market.m,arma::trans(market.trans_obj.phi - U),mu_H);

    double val = val_G + val_H;
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
        *V_out = market.trans_obj.phi - U;
    }

    if (val_out) {
        *val_out = val;
    }
    //
    return success;
}

// wrappers 

template<typename Ta>
bool max_welfare(const dse<Ta>& market, arma::mat& mu_out)
{
    bool res = max_welfare_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    
    return res;
}

template<typename Ta>
bool max_welfare(const dse<Ta>& market, arma::mat& mu_out, const double& tol_inp)
{
    bool res = max_welfare_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,NULL);
    
    return res;
}

template<typename Ta>
bool max_welfare(const dse<Ta>& market, arma::mat& mu_out, const int& max_iter_inp)
{
    bool res = max_welfare_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&max_iter_inp);
    
    return res;
}

template<typename Ta>
bool max_welfare(const dse<Ta>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp)
{
    bool res = max_welfare_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,&max_iter_inp);
    
    return res;
}

template<typename Ta>
bool max_welfare(const dse<Ta>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    bool res = max_welfare_int(market,&mu_out,NULL,NULL,&U_out,&V_out,NULL,NULL,NULL);
    
    return res;
}

template<typename Ta>
bool max_welfare(const dse<Ta>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, double& val_out, const double* tol_inp, const int* max_iter_inp)
{
    bool res = max_welfare_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,&val_out,tol_inp,max_iter_inp);
    
    return res;
}

// optimization function

template<typename Ta>
double max_welfare_opt_objfn(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data)
{
    trame_market_opt_data<Ta> *d = reinterpret_cast<trame_market_opt_data<Ta>*>(opt_data);
    //
    int nbX = d->market.nbX;
    int nbY = d->market.nbY;
    
    arma::mat U = arma::conv_to< arma::mat >::from(x_inp);
    U.reshape(nbX,nbY);

    arma::mat mu_G, mu_H;
    double val_G = d->market.arums_G.G(d->market.n,U,mu_G);
    double val_H = d->market.arums_H.G(d->market.m,arma::trans(d->market.trans_obj.phi - U),mu_H);
    //
    if (!grad.empty()) {
        arma::vec grad_vec = arma::vectorise(mu_G - mu_H.t());
        grad = arma::conv_to< std::vector<double> >::from(grad_vec);
    }
    //
    double ret = val_G + val_H;
    //
    return ret;
}

template<typename Ta>
double max_welfare_opt_objfn_2(const arma::vec& vals_inp, arma::vec* grad, void *opt_data)
{
    trame_market_opt_data<Ta> *d = reinterpret_cast<trame_market_opt_data<Ta>*>(opt_data);
    //
    int nbX = d->market.nbX;
    int nbY = d->market.nbY;
    
    arma::mat U = arma::reshape(vals_inp,nbX,nbY);

    arma::mat mu_G, mu_H;
    double val_G = d->market.arums_G.G(d->market.n,U,mu_G);
    double val_H = d->market.arums_H.G(d->market.m,arma::trans(d->market.trans_obj.phi - U),mu_H);
    //
    if (grad) {
        *grad = arma::vectorise(mu_G - mu_H.t());
    }
    //
    double ret = val_G + val_H;
    //
    return ret;
}

/*
template<typename Ta>
bool max_welfare_int(const dse<Ta>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, double* val_out, const double* tol_inp, const int* max_iter_inp)
{
    bool success = false;
    //
    // warnings
    if (!market.trans_obj.TU) {
        printf("max_welfare only works for TU transfers.\n");
        return false;
    }
    //
    int nbX = market.nbX;
    int nbY = market.nbY;

    trame_market_opt_data<Ta> opt_data;
    opt_data.market = market;

    arma::vec U_init = arma::vectorise(market.trans_obj.phi / 2.0);
    int n_pars_opt = U_init.n_elem;

    std::vector<double> sol_vec = arma::conv_to< std::vector<double> >::from(U_init);
    double obj_val = 0;
    
    success = max_welfare_nlopt(n_pars_opt,sol_vec,obj_val,NULL,NULL,max_welfare_opt_objfn<Ta>,opt_data);
    //
    // construct equilibrium objects
    arma::mat U = arma::conv_to< arma::mat >::from(sol_vec);
    U.reshape(nbX,nbY);

    Ta* arums_G = const_cast<Ta*>(&market.arums_G); // Keith: this recast is unsafe, change later
    Ta* arums_H = const_cast<Ta*>(&market.arums_H);
    
    arma::mat mu_G, mu_H;
    double val_G = arums_G->G(market.n,U,mu_G);
    double val_H = arums_H->G(market.m,arma::trans(market.trans_obj.phi - U),mu_H);

    double val = val_G + val_H;
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
        *V_out = market.trans_obj.phi - U;
    }

    if (val_out) {
        *val_out = val;
    }
    //
    return success;
}
*/

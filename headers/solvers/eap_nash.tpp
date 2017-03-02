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
 * EAP-Nash
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/15/2017
 */

// internal eap_nash

template<typename Ta>
bool eap_nash_int(const dse<Ta>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* u_out, arma::mat* v_out, const bool* xFirst_inp, const double* tol_inp, const int* max_iter_inp)
{
    bool success = false;
    //
    transfers trans_obj = market.trans_obj;

    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;

    bool xFirst = (xFirst_inp) ? *xFirst_inp : true;
    double tol = (tol_inp) ? *tol_inp : 1E-12;
    int max_iter = (max_iter_inp) ? *max_iter_inp : 10000;
    //
    arma::mat v_curr, v_next, v_err;

    if (xFirst) {
        v_curr = v_from_us(trans_obj,arma::zeros(nbX,1),NULL,NULL);
    } else {
        v_curr = arma::zeros(nbY,1);
    }

    int iter = 0;
    double err = 2*tol;
    
    while (err > tol && iter < max_iter) {
        iter++;

        v_next = update_v(trans_obj,v_curr,n,m,xFirst);
        v_err = arma::abs(v_next - v_curr);
        err = elem_max(v_err);

        v_curr = v_next;
    }
    //
    arma::mat subdiff;
    arma::mat u = u_from_vs(trans_obj,v_curr,NULL,&subdiff);

    arma::vec uv_vec = arma::join_cols(arma::vectorise(u),arma::vectorise(v_curr));

    if (u_out) {
        *u_out = u;
    }
    if (v_out) {
        *v_out = v_curr;
    }
    //
    arma::vec obj_lp = arma::vectorise(subdiff);

    arma::mat A_1_lp = arma::kron(arma::ones(1,nbY), arma::eye(nbX,nbX));
    arma::mat A_2_lp = arma::kron(arma::eye(nbY,nbY), arma::ones(1,nbX));
    arma::mat A_lp = arma::join_cols(A_1_lp, A_2_lp);

    arma::vec rhs_lp = arma::join_cols(n,m);

    int k_lp = A_lp.n_rows;
    int n_lp = A_lp.n_cols;

    int jj;
    char* sense_lp = new char[k_lp];
    for (jj=0; jj<k_lp; jj++) {
        if (uv_vec(jj) - 0 < tol) {
            sense_lp[jj] = '<';
        } else {
            sense_lp[jj] = '=';
        }
    }

    int modelSense = 1; // maximize

    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    bool LP_optimal = false;
    double val_lp = 0.0;
    //
    try {
        LP_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), A_lp.memptr(), modelSense, rhs_lp.memptr(), sense_lp, NULL, NULL, NULL, NULL, val_lp, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

        if (LP_optimal) {
            arma::mat mu = arma::reshape(sol_mat.col(0),nbX,nbY);
            //
            if (mu_out) {
                *mu_out = mu;
            }
            
            if (mu_x0_out) {
                *mu_x0_out = n - arma::sum(mu,1);
            }
            if (mu_0y_out) {
                *mu_0y_out = m - arma::trans(arma::sum(mu,0));
            }
            //
            success = true;
        } else {
            std::cout << "Non-optimal value found during optimization" << std::endl;
        }
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //
    return success;
}

// wrappers

template<typename Ta>
bool eap_nash(const dse<Ta>& market, arma::mat& mu_out)
{
    bool res = eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    
    return res;
}

template<typename Ta>
bool eap_nash(const dse<Ta>& market, arma::mat& mu_out, const bool& xFirst_inp)
{
    bool res = eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,&xFirst_inp,NULL,NULL);
    
    return res;
}

template<typename Ta>
bool eap_nash(const dse<Ta>& market, arma::mat& mu_out, const double& tol_inp)
{
    bool res = eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,NULL);
    
    return res;
}

template<typename Ta>
bool eap_nash(const dse<Ta>& market, arma::mat& mu_out, const int& max_iter_inp)
{
    bool res = eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&max_iter_inp);
    
    return res;
}

template<typename Ta>
bool eap_nash(const dse<Ta>& market, arma::mat& mu_out, const bool& xFirst_inp, const double& tol_inp)
{
    bool res = eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,&xFirst_inp,&tol_inp,NULL);
    
    return res;
}

template<typename Ta>
bool eap_nash(const dse<Ta>& market, arma::mat& mu_out, const bool& xFirst_inp, const int& max_iter_inp)
{
    bool res = eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&max_iter_inp);
    
    return res;
}

template<typename Ta>
bool eap_nash(const dse<Ta>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp)
{
    bool res = eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,&max_iter_inp);
    
    return res;
}

template<typename Ta>
bool eap_nash(const dse<Ta>& market, arma::mat& mu_out, const bool& xFirst_inp, const double& tol_inp, const int& max_iter_inp)
{
    bool res = eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,&xFirst_inp,&tol_inp,&max_iter_inp);
    
    return res;
}

template<typename Ta>
bool eap_nash(const dse<Ta>& market, arma::mat& mu_out, arma::mat& u_out, arma::mat& v_out)
{
    bool res = eap_nash_int(market,&mu_out,NULL,NULL,&u_out,&v_out,NULL,NULL,NULL);
    
    return res;
}

template<typename Ta>
bool eap_nash(const dse<Ta>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::vec& u_out, arma::vec& v_out, const bool* xFirst_inp, const double* tol_inp, const int* max_iter_inp)
{
    bool res = eap_nash_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&u_out,&v_out,xFirst_inp,tol_inp,max_iter_inp);
    
    return res;
}

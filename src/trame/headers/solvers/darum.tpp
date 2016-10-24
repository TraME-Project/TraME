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
 * darum
 *
 * Keith O'Hara
 * 08/16/2016
 */

// internal darum

template<typename Ta>
bool darum_int(const dse<Ta>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, const double* tol_inp, const int* max_iter_inp)
{
    bool success = false;
    //
    // warnings
    if (!market.trans_obj.NTU) {
        printf("darum only works for NTU transfers.\n");
        return false;
    }
    //
    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;

    arma::mat alpha = market.trans_obj.alpha;
    arma::mat gamma = market.trans_obj.gamma;

    Ta* arums_G = const_cast<Ta*>(&market.arums_G); // Keith: this recast is unsafe, change later
    Ta* arums_H = const_cast<Ta*>(&market.arums_H);

    arma::mat mu_NR = arma::max(n * arma::ones(1,nbY), arma::ones(nbX,1) * m.t());
    //
    double tol = (tol_inp) ? *tol_inp : 1E-12;
    int max_iter = (max_iter_inp) ? *max_iter_inp : 10000;
    //
    int iter = 0;
    double err = 2*tol;
    
    arma::mat U_P, U_D, mu_P, mu_D, mu_diff;

    while (err > tol && iter < max_iter) {
        iter++;
        //
        arums_G->Gbar(alpha,mu_NR,n,U_P,mu_P);
        arums_H->Gbar(gamma.t(),mu_P.t(),m,U_D,mu_D);

        mu_diff = mu_P - mu_D.t();
        mu_NR -= mu_diff;
        //
        err = elem_max(arma::abs(mu_diff));
    }
    
    if (err <= tol && iter < max_iter) {
        success = true;
    }
    //
    // return equilibrium objects
    if (mu_out) {
        *mu_out = mu_D.t();
    }

    if (mu_x0_out) {
        *mu_x0_out = n - arma::sum(mu_D.t(),1);
    }
    if (mu_0y_out) {
        *mu_0y_out = m - arma::trans(arma::sum(mu_D.t(),0));
    }

    if (U_out) {
        *U_out = U_P;
    }
    if (V_out) {
        *V_out = U_D.t();
    }
    //
    return success;
}

// wrappers 

template<typename Ta>
bool darum(const dse<Ta>& market, arma::mat& mu_out)
{
    bool res = darum_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL);
    
    return res;
}

template<typename Ta>
bool darum(const dse<Ta>& market, arma::mat& mu_out, const double& tol_inp)
{
    bool res = darum_int(market,&mu_out,NULL,NULL,NULL,NULL,&tol_inp,NULL);
    
    return res;
}

template<typename Ta>
bool darum(const dse<Ta>& market, arma::mat& mu_out, const int& max_iter_inp)
{
    bool res = darum_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&max_iter_inp);
    
    return res;
}

template<typename Ta>
bool darum(const dse<Ta>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp)
{
    bool res = darum_int(market,&mu_out,NULL,NULL,NULL,NULL,&tol_inp,&max_iter_inp);
    
    return res;
}

template<typename Ta>
bool darum(const dse<Ta>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, const double* tol_inp, const int* max_iter_inp)
{
    bool res = darum_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,tol_inp,max_iter_inp);
    
    return res;
}

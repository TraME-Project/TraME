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
 * jacobi solver
 *
 * Keith O'Hara
 * 08/25/2016
 */

template<typename Ta>
bool jacobi(dse<Ta> market, bool xFirst, arma::mat* w_low_inp, arma::mat* w_up_inp, double* tol_inp, arma::mat& mu, arma::vec& mux0, arma::vec& mu0y, arma::mat& U, arma::mat& V)
{
    if (market.need_norm) {
        printf("Jacobi does not yet allow for the case without unmatched agents.\n");
        return false;
    }
    //
    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;

    Ta arums_G = market.arums_G;
    Ta arums_H = market.arums_H;

    transfers trans_obj = market.trans_obj;
    //
    printf("jacobi 1\n");
    arma::mat w;
    if (!w_up_inp) {
        w = w_upper_bound(market);
    } else {
        w = *w_up_inp;

        arma::mat temp_UW = trans_obj.UW(w,NULL,NULL);
        arma::mat temp_VW = trans_obj.VW(w,NULL,NULL);
        arums_G.U = temp_UW;
        arums_H.U = temp_VW.t();

        arums_G.G(n);
        arums_H.G(m);

        arma::mat Z = arums_G.mu_sol - arma::trans(arums_H.mu_sol);

        if (elem_min(Z) < 0) {
            printf("jacobi: w_up provided not an actual upper bound.\n");
            return false;
        }
    }
    //
    printf("jacobi 2\n");
    arma::mat w_low;
    if (!w_low_inp) {
        dse<Ta> market_trans = market;
        market_trans.trans();

        w_low = -arma::trans(w_upper_bound(market_trans));
    } else {
        w_low = *w_low_inp;

        arma::mat temp_UW = trans_obj.UW(w_low,NULL,NULL);
        arma::mat temp_VW = trans_obj.VW(w_low,NULL,NULL);
        arums_G.U = temp_UW;
        arums_H.U = temp_VW.t();

        arums_G.G(n);
        arums_H.G(m);

        arma::mat Z = arums_G.mu_sol - arma::trans(arums_H.mu_sol);

        if (elem_max(Z) > 0) {
            printf("jacobi: w_low provided not an actual upper bound.\n");
            return false;
        }
    }
    //
    U = trans_obj.UW(w,NULL,NULL);
    V = trans_obj.VW(w,NULL,NULL);

    arums_G.U = U;
    arums_H.U = V.t();

    arums_G.G(n);
    arums_H.G(m);

    arma::mat Z = arums_G.mu_sol - arma::trans(arums_H.mu_sol);
    //
    bool success = false;

    int iter = 0;
    int max_iter = 10000;

    double tol = 1E-04;
    double err = 2*tol;

    


    //
    return true;
}

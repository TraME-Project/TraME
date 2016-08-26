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

template<typename Ta>
bool darum(dse<Ta> market, bool xFirst, double* tol_inp, arma::mat& mu, arma::vec& mux0, arma::vec& mu0y, arma::mat& U, arma::mat& V)
{
    if (!market.trans_obj.NTU) {
        printf("darum only works for NTU transfers.\n");
        return false;
    }
    //
    arma::mat alpha = market.trans_obj.alpha;
    arma::mat gamma = market.trans_obj.gamma;

    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;

    arma::mat mu_NR = arma::max(n * arma::ones(1,nbY), arma::ones(nbX,1) * m.t());
    //
    double tol;
    if (tol_inp) {
        tol = *tol_inp;
    } else {
        tol = 1E-12;
    }

    bool success = false;
    double err = 2*tol;
    int iter = 0;
    int max_iter = 10000;
    //
    arma::mat U_P, U_D, mu_P, mu_D, mu_diff;

    while (err > tol && iter < max_iter) {
        iter++;
        //
        market.arums_G.Gbar(alpha,mu_NR,n,U_P,mu_P);
        market.arums_H.Gbar(gamma.t(),mu_P.t(),m,U_D,mu_D);

        mu_diff = mu_P - mu_D.t();
        mu_NR -= mu_diff;
        //
        err = arma::as_scalar(arma::max(arma::max(arma::abs(mu_diff))));
    }
    //
    if (err <= tol && iter < max_iter) {
        success = true;
    }
    mu = mu_D.t();

    mux0 = n - arma::sum(mu_D.t(),1);
    mu0y = m - arma::trans(arma::sum(mu_D.t(),0));

    U = U_P;
    V = U_D.t();
    //
    return success;
}

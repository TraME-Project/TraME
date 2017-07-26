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
 * auxiliary functions for equilibrium solvers
 *
 * Keith O'Hara
 * 08/23/2016
 */

/*
 * w upper bound
 * used in jacobi and arc_newton solvers
 *
 * Keith O'Hara
 * 03/22/2017
 *
 * This version:
 * 07/26/2017
 */

template<typename Tg, typename Th, typename Tt>
arma::mat
w_upper_bound(const dse<Tg,Th,Tt>& market)
{
    const int nbX = market.nbX;
    const int nbY = market.nbY;

    int transfers_type = market.trans_obj.transfers_type;
    //
    int iter = 0;
    const int max_iter = 1000;
    
    double k = 1.0, Z_min_val = -10.0;
    arma::uvec x_ind(1);
    arma::vec mu_cond_x(nbY,1);

    arma::mat U_star_x, w(nbX,nbY), U(nbX,nbY), V(nbX,nbY);

    while (Z_min_val < 0 && iter < max_iter) {
        iter ++;

        for (int x=0; x < nbX; x++) {
            x_ind(0) = x;

            if (transfers_type == 1) {
                double mu_fill = 1.0 / (std::pow(2.0,-k) + (double) nbY);
                mu_cond_x.fill(mu_fill);

                market.arums_G.Gstarx(mu_cond_x,U_star_x,x);
                U.row(x) = U_star_x.t();

                w.row(x) = market.trans_obj.WU(U_star_x.t(),&x_ind,NULL);
                V.row(x) = market.trans_obj.VW(w.row(x),&x_ind,NULL);
            } else if (transfers_type == 2) {
                w.row(x).fill(std::pow(2,k));

                U.row(x) = market.trans_obj.UW(w.row(x),&x_ind,NULL);
                V.row(x) = market.trans_obj.VW(w.row(x),&x_ind,NULL);
            } else {
                printf("w_upper_bound error: unrecognized transfers_type");
                return w;
            }
        }
        //
        arma::mat mu_G, mu_H;
        market.arums_G.G(market.n,U,mu_G);
        market.arums_H.G(market.m,V.t(),mu_H);

        arma::mat Z = mu_G - mu_H.t();
        Z_min_val = elem_min(Z);
        //
        k *= 2;
    }
    //
    return w;
}

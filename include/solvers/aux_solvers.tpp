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
 * auxiliary functions for equilibrium solvers
 *
 * Keith O'Hara
 * 08/23/2016
 *
 * This version:
 * 02/04/2018
 */

/*
 * w upper bound
 * used in jacobi and arc_newton solvers
 */

template<typename Tg, typename Th, typename Tt>
arma::mat
w_upper_bound(const dse<Tg,Th,Tt>& market)
{
    const uint_t nbX = market.nbX;
    const uint_t nbY = market.nbY;

    const uint_t transfers_type = market.trans_obj.transfers_type;

    //

    uint_t iter = 0;
    const uint_t max_iter = 1000;
    
    double k = 1.0, Z_min_val = -10.0;
    arma::uvec x_ind(1);
    arma::vec mu_cond_x(nbY,1);

    arma::mat U_star_x, w(nbX,nbY), U(nbX,nbY), V(nbX,nbY);

    while (Z_min_val < 0 && iter < max_iter) {
        iter++;

        for (uint_t x=0; x < nbX; x++)
        {
            x_ind(0) = x;

            if (transfers_type == 1) {
                double mu_fill = 1.0 / (std::pow(2.0,-k) + static_cast<double>(nbY));
                mu_cond_x.fill(mu_fill);

                market.arums_G.Gstarx(mu_cond_x,U_star_x,x);
                U.row(x) = U_star_x.t();

                w.row(x) = market.trans_obj.WU(U_star_x.t(),&x_ind,nullptr);
                V.row(x) = market.trans_obj.VW(w.row(x),&x_ind,nullptr);
            } else if (transfers_type == 2) {
                w.row(x).fill(std::pow(2,k));

                U.row(x) = market.trans_obj.UW(w.row(x),&x_ind,nullptr);
                V.row(x) = market.trans_obj.VW(w.row(x),&x_ind,nullptr);
            } else {
                printf("w_upper_bound error: unrecognized transfers_type");
                return w;
            }
        }
        
        //

        arma::mat mu_G, mu_H;
        market.arums_G.G(market.n,U,mu_G);
        market.arums_H.G(market.m,V.t(),mu_H);

        Z_min_val = elem_min(mu_G - mu_H.t());

        //

        k *= 2;
    }
    //
    return w;
}

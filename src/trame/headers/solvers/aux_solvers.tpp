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
 * auxiliary functions for equilibrium solvers
 *
 * Keith O'Hara
 * 08/23/2016
 */

/*
 * w upper bound
 * used in 'jacobi' solver
 *
 * Keith O'Hara
 * 08/25/2016
 */

template<typename Ta>
arma::mat w_upper_bound(dse<Ta> market)
{
    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;

    Ta arums_G = market.arums_G;
    Ta arums_H = market.arums_H;

    transfers trans_obj = market.trans_obj;
    int transfers_type = trans_obj.transfers_type;
    //
    //bool success = false;
    int iter = 0;
    int max_iter = 1000;
    double Z_min_val = -10.0;
    
    int x;
    double k = 1.0;
    double mu_fill;
    arma::uvec x_ind(1);
    arma::vec mu_cond_x(nbY,1);
    arma::mat U_star_x, Z;

    arma::mat w(nbX,nbY);
    arma::mat U(nbX,nbY);
    arma::mat V(nbX,nbY);

    while (Z_min_val < 0 && iter < max_iter) {
        iter ++;

        for (x=0; x < nbX; x++) {
            x_ind(0) = x;

            if (transfers_type==1) {
                mu_fill = 1.0 / (std::pow(2.0,-k) + (double) nbY);
                mu_cond_x.fill(mu_fill);

                arums_G.Gstarx(U_star_x,mu_cond_x,x);
                U.row(x) = arma::trans(U_star_x);

                w.row(x) = trans_obj.WU(U_star_x.t(),&x_ind,NULL);
                V.row(x) = trans_obj.VW(w.row(x),&x_ind,NULL);
            } else if (transfers_type==2) {
                w.row(x).fill(std::pow(2,k));

                U.row(x) = trans_obj.UW(w.row(x),&x_ind,NULL);
                V.row(x) = trans_obj.VW(w.row(x),&x_ind,NULL);
            } else {
                printf("w_upper_bound error: unrecognized transfers_type");
                return w;
            }
        }
        //
        arums_G.U = U;
        arums_H.U = V.t();

        arums_G.G(n);
        arums_H.G(m);

        Z = arums_G.mu_sol - arma::trans(arums_H.mu_sol);
        Z_min_val = elem_min(Z);
        //
        k *= 2;
    }
    //
    return w;
}

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
 * inverse PWA function
 *
 * Keith O'Hara
 * 08/24/2016
 *
 * Modified: 
 *   09/15/16, Alfred Galichon
 *   07/26/17, Keith O'Hara
 */

#include "ancillary/ancillary.hpp"

arma::vec 
trame::inv_pwa(const arma::vec& a, const arma::mat& B, const arma::mat& C, const double k)
{
    const int nb_X = a.n_elem;
    const int nb_Y = B.n_cols;
    //
    arma::vec vals = arma::zeros(nb_X,1);
    
#ifdef TRAME_USE_OMP
    #pragma omp parallel for
#endif
    for (int x=0; x < nb_X; x++) {
        const arma::vec Bx = B.row(x).t();  // transpose to ensure that Bx is a column-vector
        const arma::uvec Bx_sort_ind = arma::sort_index(Bx);
        const arma::vec b = Bx.elem(Bx_sort_ind);
        //
        arma::vec small_C = C.row(x).t();
        small_C = small_C.elem(Bx_sort_ind);

        const double small_C_sum = arma::accu(small_C);
        //
        int y_low = 1;
        int y_up  = nb_Y;

        arma::vec b_low(nb_Y), b_mid(nb_Y);
        
        while (y_up > y_low + 1) {
            int y_mid = y_low + std::floor((y_up - y_low)/2.0);
            
            b_low.fill(b(y_low-1));
            b_mid.fill(b(y_mid-1));
            
            const double lhs = k * b(y_mid-1) + arma::accu(small_C % arma::min(b_mid,b));
            //
            if (lhs == a(x)) {
                y_low = y_mid;
                y_up = y_mid;
            } else if (lhs > a(x)) {
                y_up = y_mid;
            } else {
                y_low = y_mid + 1;
            }
        }
        //
        if ((y_low==1) && ( k * b(y_low-1) + arma::accu(small_C % arma::min(b_low,b)) >= a(x) )) {
            vals(x) = a(x) / (k + small_C_sum);
        } else {
            arma::vec b_up(nb_Y);
            b_up.fill(b(y_up-1));

            if ((y_up==nb_Y) && ( k * b(y_up-1) + arma::accu(small_C % arma::min(b_up,b)) <= a(x) )) {
                vals(x) = (a(x) - arma::accu(small_C % b)) / k ;            
            } else {
                int y_incl = y_low - 1;

                const arma::vec temp_vec = small_C % b;
                const double term_1 = a(x) - arma::accu(temp_vec.rows(0,y_incl));
                const double term_2 = k + small_C_sum - arma::accu(small_C.rows(0,y_incl));
                
                vals(x) = term_1 / term_2;
            }
        }
    }
    //
    return vals;
}

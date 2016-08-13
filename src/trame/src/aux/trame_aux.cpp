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
 * Aux functions
 *
 * Keith O'Hara
 * 08/08/2016
 */

#include "aux/trame_aux.hpp"

// Find indices that correspond to maximum values 
arma::uvec which_max(const arma::mat* X, int which_dim)
{
    int i,j;
    
    int n = X->n_rows;
    int k = X->n_cols;

    arma::uvec max_vec;
     
    int max_ind = 0;
    double max_val = 0;
     
    if (which_dim==0) { // each column
        max_vec.set_size(k);
        for (j=0; j<k; j++) {
            max_val = (*X)(0,j);
            max_ind = 0;
            for (i=1; i<n; i++) {
                if ((*X)(i,j)>max_val) {
                    max_val = (*X)(i,j);
                    max_ind = i;
                }
            }
            max_vec(j) = max_ind;
        }
    } else { // each row
        max_vec.set_size(n);
        for (i=0; i<n; i++) {
            max_val = (*X)(i,0);
            max_ind = 0;
            for (j=0; j<k; j++) {
                if ((*X)(i,j) > max_val) {
                    max_val = (*X)(i,j);
                    max_ind = j;
                }
            }
            max_vec(i) = max_ind;
        }
    }
    //
    return max_vec;
}

// Generate an equally-spaced grid of integers
arma::uvec uvec_linspace (int a, int b)
{
    int i;
    int n_points = b - a + 1;
    
    arma::uvec ret(n_points);
    //
    for (i = 0; i<n_points; i++) {
        ret(i) = a + i;
    }
    //
    return ret;
}

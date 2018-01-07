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
 * test auxiliary functions
 *
 * Keith O'Hara
 * 05/17/2016
 *
 * This version:
 * 08/20/2017
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    const int dim_1 = 5;
    const int dim_2 = 4;
    const int dim_3 = 3;

    //
    printf("\n*===================   Start of trame_aux test   ===================*\n");
    printf("\n");

    // LogSumExp

    arma::vec lse_vec(3);
    lse_vec(0) = 2.0; lse_vec(1) = 3.0; lse_vec(2) = 4.0;

    double lse_val = trame::lse(lse_vec);
    std::cout << "lse_val = " << lse_val << std::endl;
    
    // which max

    arma::mat X = arma::randu(dim_1,dim_2);

    trame::which_max(X, static_cast<int>(0));
    trame::which_max(X, static_cast<int>(1));

    // unit_vec

    trame::unit_vec(static_cast<int>(1), static_cast<int>(5));

    // uvec_linspace

    trame::uvec_linspace(-1,1);

    // byrow

    arma::mat br_vec_1 = arma::randu(dim_1,1);
    arma::mat br_vec_2 = arma::randu(1,dim_2);

    arma::mat br_X = arma::randu(dim_1,dim_2);

    trame::byrow(br_vec_1,dim_1,dim_2);
    trame::byrow(br_vec_2,dim_1,dim_2);
    trame::byrow(br_X,dim_1,dim_2);

    trame::byrow(br_vec_1,dim_1+1,dim_2+1); // error

    // elem-by-elem operations

    arma::mat eX_mat_1 = arma::randu(dim_1,dim_2);
    arma::mat eX_mat_2 = arma::randu(dim_1,dim_2);

    arma::mat eX_vec_1 = arma::randu(dim_1,1);
    arma::mat eX_vec_2 = arma::randu(1,dim_2);

    trame::elem_add(eX_mat_1,eX_mat_2);
    trame::elem_add(eX_vec_1,eX_mat_2);
    trame::elem_add(eX_mat_2,eX_vec_1);
    trame::elem_add(eX_mat_1,eX_vec_2);
    trame::elem_add(eX_vec_2,eX_mat_1);
    trame::elem_add(eX_vec_1,eX_vec_2);

    trame::elem_sub(eX_mat_1,eX_mat_2);
    trame::elem_sub(eX_vec_1,eX_mat_2);
    trame::elem_sub(eX_mat_2,eX_vec_1);
    trame::elem_sub(eX_mat_1,eX_vec_2);
    trame::elem_sub(eX_vec_2,eX_mat_1);
    trame::elem_sub(eX_vec_1,eX_vec_2);

    trame::elem_prod(eX_mat_1,eX_mat_2);
    trame::elem_prod(eX_vec_1,eX_mat_2);
    trame::elem_prod(eX_mat_2,eX_vec_1);
    trame::elem_prod(eX_mat_1,eX_vec_2);
    trame::elem_prod(eX_vec_2,eX_mat_1);
    trame::elem_prod(eX_vec_1,eX_vec_2);

    trame::elem_div(eX_mat_1,eX_mat_2);
    trame::elem_div(eX_vec_1,eX_mat_2);
    trame::elem_div(eX_mat_2,eX_vec_1);
    trame::elem_div(eX_mat_1,eX_vec_2);
    trame::elem_div(eX_vec_2,eX_mat_1);
    trame::elem_div(eX_vec_1,eX_vec_2);

    trame::elem_min(eX_mat_1);
    trame::elem_min(eX_mat_1,eX_mat_2);
    trame::elem_min(eX_mat_1,0.5);
    trame::elem_min(0.5,eX_mat_2);
    trame::elem_min(0.4,0.5);

    trame::elem_max(eX_mat_1);
    trame::elem_max(eX_mat_1,eX_mat_2);
    trame::elem_max(eX_mat_1,0.5);
    trame::elem_max(0.5,eX_mat_2);
    trame::elem_max(0.4,0.5);

    // cube_sum

    arma::cube cs_X = arma::randu(dim_1,dim_2,dim_3);

    trame::cube_sum(cs_X,0);
    trame::cube_sum(cs_X,1);
    trame::cube_sum(cs_X,2); // error

    // cube_to_mat

    arma::cube ctm_X = arma::randu(dim_1,dim_2,dim_3);

    trame::cube_to_mat(ctm_X);

    //
    printf("\n*===================   End of trame_aux test   ===================*\n");
    printf("\n");
    //

    std::chrono::time_point<std::chrono::system_clock>end = std::chrono::system_clock::now();
        
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        
    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
    //
    return 0;
}

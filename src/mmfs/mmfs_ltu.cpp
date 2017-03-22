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
 * Linearly Transferable Utility (LTU) Marriage Matching Functions (MMFs) class
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 03/12/2017
 */

#include "trame.hpp"

void trame::mmfs::ltu::build(const arma::mat& lambda_LTU, const arma::mat& phi_LTU, bool need_norm_LTU)
{
    need_norm = need_norm_LTU;

    nbX = lambda_LTU.n_rows;
    nbY = lambda_LTU.n_cols;
    nbParams = 2*nbX*nbY;

    lambda = lambda_LTU;
    phi = phi_LTU;

    aux_zeta = 1 - lambda_LTU;
    aux_phi_exp = arma::exp(phi_LTU);
}

void trame::mmfs::ltu::trans()
{
    int nbX_temp = nbX;

    nbX = nbY;
    nbY = nbX_temp;
    //
    arma::mat lambda_temp;

    lambda = aux_zeta.t();
    aux_zeta = lambda_temp.t();

    phi = phi.t();
    aux_phi_exp = aux_phi_exp.t();
    //
}

//
// MFE-related functions

arma::mat trame::mmfs::ltu::M(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    arma::mat ret = this->M(a_xs,b_ys,NULL,NULL);
    //
    return ret;
}

arma::mat trame::mmfs::ltu::M(const arma::mat& a_xs, const arma::mat& b_ys, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = arma::exp(elem_prod(lambda(x_ind,y_ind), arma::log(a_xs)));
    arma::mat term_2 = arma::trans(arma::exp( elem_prod(arma::trans(aux_zeta(x_ind,y_ind)), arma::log(b_ys)) ));
    arma::mat term_3 = aux_phi_exp(x_ind,y_ind);

    arma::mat ret = term_1 % term_2 % term_3;
    //
    return ret;
}

arma::mat trame::mmfs::ltu::M(const double& a_xs, const arma::mat& b_ys, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, (int) nbX-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, (int) nbY-1);
    //
    arma::mat term_1 = arma::exp(lambda(x_ind,y_ind) * std::log(a_xs));
    arma::mat term_2 = arma::trans(arma::exp( elem_prod(arma::trans(aux_zeta(x_ind,y_ind)), arma::log(b_ys)) ));
    arma::mat term_3 = aux_phi_exp(x_ind,y_ind);

    arma::mat ret = term_1 % term_2 % term_3;
    //
    return ret;
}

arma::mat trame::mmfs::ltu::M(const arma::mat& a_xs, const double& b_ys, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, (int) nbX-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, (int) nbY-1);
    //
    arma::mat term_1 = arma::exp(elem_prod(lambda(x_ind,y_ind), arma::log(a_xs)));
    arma::mat term_2 = arma::exp(aux_zeta(x_ind,y_ind) * std::log(b_ys));
    arma::mat term_3 = aux_phi_exp(x_ind,y_ind);
        
    arma::mat ret = term_1 % term_2 % term_3;
    //
    return ret;
}

arma::mat trame::mmfs::ltu::Mx0(const arma::mat& a_x)
const
{
    return a_x;
}

arma::mat trame::mmfs::ltu::M0y(const arma::mat& b_y)
const
{
    return b_y;
}

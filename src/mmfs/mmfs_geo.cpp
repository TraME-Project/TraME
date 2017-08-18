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
 * Geometric Matching Market Functions (MMFs) class
 * Corresponds to TU transfers class
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 07/24/2017
 */

#include "ancillary/ancillary.hpp"
#include "mmfs/mmfs.hpp"

void 
trame::mmfs::geo::build(const arma::mat& phi_inp, const bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = phi_inp.n_rows;
    nbY = phi_inp.n_cols;
    dim_params = nbX*nbY;

    phi = phi_inp;
    aux_phi_exp = arma::exp(phi_inp / 2.0);
}

void 
trame::mmfs::geo::trans()
{
    std::swap(nbX,nbY);
    //
    arma::inplace_trans(phi);
    arma::inplace_trans(aux_phi_exp);
}

//
// MFE-related functions
//

//
// matching function

arma::mat 
trame::mmfs::geo::M(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    return this->M(a_xs,b_ys,nullptr,nullptr);
}

arma::mat 
trame::mmfs::geo::M(const arma::mat& a_xs, const arma::mat& b_ys, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1); 
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    const arma::mat term_1 = aux_phi_exp(x_ind,y_ind);
    const arma::mat term_2 = arma::sqrt(a_xs * b_ys.t());
    //
    return term_1 % term_2;
}

arma::mat 
trame::mmfs::geo::M(const double a_xs, const arma::mat& b_ys, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    const arma::mat term_1 = aux_phi_exp(x_ind,y_ind);
    const arma::mat term_2 = arma::sqrt(a_xs * b_ys.t());
    //
    return term_1 % term_2;
}

arma::mat 
trame::mmfs::geo::M(const arma::mat& a_xs, const double b_ys, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    const arma::mat term_1 = aux_phi_exp(x_ind,y_ind);
    const arma::mat term_2 = arma::sqrt(a_xs * b_ys);
    //
    return term_1 % term_2;
}

//

arma::mat 
trame::mmfs::geo::dmu_x0(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    const arma::mat term_1 = aux_phi_exp / 2.0;
    const arma::mat term_2 = arma::sqrt( (1/a_xs) * b_ys.t());
    //
    return term_1 % term_2;
}

arma::mat 
trame::mmfs::geo::dmu_0y(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    const arma::mat term_1 = aux_phi_exp / 2.0;
    const arma::mat term_2 = arma::sqrt( a_xs * (1/b_ys.t()) );
    //
    return term_1 % term_2;
}

arma::mat 
trame::mmfs::geo::dparams_M(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    return this->dparams_M(a_xs,b_ys,nullptr);
}

arma::mat 
trame::mmfs::geo::dparams_M(const arma::mat& a_xs, const arma::mat& b_ys, const arma::mat* delta_params_M)
const
{
    const arma::mat der_1 = arma::sqrt( a_xs * b_ys.t() );

    arma::mat ret;

    if (delta_params_M) {
        const arma::mat delta_params_1 = arma::reshape((*delta_params_M),nbX,nbY);

        ret = arma::vectorise(delta_params_1 % der_1);
    } else {
        ret = arma::diagmat(arma::vectorise(der_1));
    }
    //
    return ret;
}

arma::mat 
trame::mmfs::geo::Mx0(const arma::mat& a_x)
const
{
    return a_x;
}

arma::mat 
trame::mmfs::geo::M0y(const arma::mat& b_y)
const
{
    return b_y;
}

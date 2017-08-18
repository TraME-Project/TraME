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
 * Non-Transferable Utility (NTU) transfers class
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 07/25/2017
 */

#include "trame.hpp"

void
trame::transfers::ntu::build(const arma::mat& alpha_inp, const arma::mat& gamma_inp, const bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = alpha_inp.n_rows;
    nbY = alpha_inp.n_cols;
    dim_params = 2*nbX*nbY;

    alpha = alpha_inp;
    gamma = gamma_inp;

    aux_alpha_exp = arma::exp(alpha_inp);
    aux_gamma_exp = arma::exp(gamma_inp);
}

void
trame::transfers::ntu::trans()
{
    std::swap(nbX,nbY);
    //
    arma::inplace_trans(alpha);
    arma::inplace_trans(gamma);
    alpha.swap(gamma);

    arma::inplace_trans(aux_alpha_exp);
    arma::inplace_trans(aux_gamma_exp);
    aux_alpha_exp.swap(aux_gamma_exp);
}

void
trame::transfers::ntu::gen_mmf(mmfs::min& mmf_obj)
const
{
    mmf_obj.build(alpha,gamma,need_norm);
}

trame::mmfs::min
trame::transfers::ntu::gen_mmf()
const
{
    mmfs::min mmf_obj;
    this->gen_mmf(mmf_obj);
    //
    return mmf_obj;
}

//
// DSE-related functions
//

arma::mat 
trame::transfers::ntu::Psi(const arma::mat& U, const arma::mat& V)
const
{
    return this->Psi(U,V,nullptr,nullptr);
}

//
// Implicit Parameterization

arma::mat 
trame::transfers::ntu::Psi(const arma::mat& U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return elem_max(U - alpha(x_ind,y_ind), V - gamma(x_ind,y_ind));
}

arma::mat 
trame::transfers::ntu::Psi(const double U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return elem_max(U - alpha(x_ind,y_ind), V - gamma(x_ind,y_ind));
}

arma::mat 
trame::transfers::ntu::Psi(const arma::mat& U, const double V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return elem_max(U - alpha(x_ind,y_ind), V - gamma(x_ind,y_ind));
}

double 
trame::transfers::ntu::Psi(const double U, const double V, const int x_ind, const int y_ind)
const
{
    return std::max(U - alpha(x_ind,y_ind), V - gamma(x_ind,y_ind));
}

// Derivative of Psi wrt u

arma::mat 
trame::transfers::ntu::du_Psi(const arma::mat& U, const arma::mat& V)
const
{
    return this->du_Psi(U,V,nullptr,nullptr);
}

arma::mat 
trame::transfers::ntu::du_Psi(const arma::mat& U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = arma::zeros(x_ind.n_elem,y_ind.n_elem);
    ret.elem( arma::find(U - alpha(x_ind,y_ind) >= V - gamma(x_ind,y_ind)) ).ones();
    //
    return ret;
}

arma::mat 
trame::transfers::ntu::du_Psi(const double U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = arma::zeros(x_ind.n_elem,y_ind.n_elem);
    ret.elem( arma::find(U - alpha(x_ind,y_ind) >= V - gamma(x_ind,y_ind)) ).ones();
    //
    return ret;
}

arma::mat 
trame::transfers::ntu::du_Psi(const arma::mat& U, const double V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = arma::zeros(x_ind.n_elem,y_ind.n_elem);
    ret.elem( arma::find(U - alpha(x_ind,y_ind) >= V - gamma(x_ind,y_ind)) ).ones();
    //
    return ret;
}

// dparams

arma::mat 
trame::transfers::ntu::dparams_Psi(const arma::mat& U, const arma::mat& V)
const
{
    return this->dparams_Psi(U,V,nullptr);
}

arma::mat 
trame::transfers::ntu::dparams_Psi(const arma::mat& U, const arma::mat& V, const arma::mat* dparams)
const
{
    arma::mat ret(nbX,nbY);
    //
    arma::vec du_psi = arma::vectorise(du_Psi(U,V,nullptr,nullptr));

    if (!dparams) {
        ret = - arma::join_rows(arma::diagmat(du_psi),arma::diagmat(1 - du_psi));
    } else {
        const arma::mat dparams_1 = dparams->rows(0,nbX*nbY-1);
        const arma::mat dparams_2 = dparams->rows(nbX*nbY,2*nbX*nbY-1);

        ret = - arma::vectorise(du_psi % dparams_1 + (1 - du_psi) % dparams_2);
    }
    //
    return ret;
}

//
// Explicit Parameterization

// Ucal and Vcal

arma::mat 
trame::transfers::ntu::Ucal(const arma::mat& vs)
const
{
    return this->Ucal(vs,nullptr,nullptr);
}

arma::mat 
trame::transfers::ntu::Ucal(const arma::mat& vs, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = arma::zeros(x_ind.n_elem,y_ind.n_elem); // shouldn't be zero, should return an error
    //
    printf("no Ucal defined for NTU case\n");
    //
    return ret;
}

double 
trame::transfers::ntu::Ucal(const double vs, const int xs, const int ys)
const
{
    double ret = 0; // shouldn't be zero, should return an error
    //
    printf("no Ucal defined for NTU case\n");
    //
    return ret;
}

arma::mat 
trame::transfers::ntu::Vcal(const arma::mat& us)
const
{
    return this->Vcal(us,nullptr,nullptr);
}

arma::mat 
trame::transfers::ntu::Vcal(const arma::mat& us, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret(x_ind.n_elem,y_ind.n_elem); // shouldn't be zero, should return an error
    //
    printf("no Vcal defined for NTU case\n");
    //
    return ret;
}

double 
trame::transfers::ntu::Vcal(const double us, const int xs, const int ys)
const
{
    double ret = 0; // shouldn't be zero, should return an error
    //
    printf("no Vcal defined for NTU case\n");
    //
    return ret;
}

// UW and VW

arma::mat 
trame::transfers::ntu::UW(const arma::mat& Ws)
const
{
    return this->UW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::ntu::UW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return - Psi(0.0,-Ws,&x_ind,&y_ind);
}

double 
trame::transfers::ntu::UW(const double Ws, const int x_ind, const int y_ind)
const
{
    return - Psi((double) 0.0,-Ws,x_ind,y_ind);
}

arma::mat 
trame::transfers::ntu::VW(const arma::mat& Ws)
const
{
    return this->VW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::ntu::VW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return - Psi(Ws,0.0,&x_ind,&y_ind);
}

double 
trame::transfers::ntu::VW(const double Ws, const int x_ind, const int y_ind)
const
{
    return - Psi(Ws,(double) 0.0,x_ind,y_ind);
}

// dw

arma::mat 
trame::transfers::ntu::dw_UW(const arma::mat& Ws)
const
{
    return this->dw_UW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::ntu::dw_UW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return 1.0 - du_Psi(0.0,-Ws,&x_ind,&y_ind);
}

arma::mat 
trame::transfers::ntu::dw_VW(const arma::mat& Ws)
const
{
    return this->dw_VW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::ntu::dw_VW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return - du_Psi(Ws,0.0,&x_ind,&y_ind);
}

// WU and WV

arma::mat 
trame::transfers::ntu::WU(const arma::mat& Us)
const
{
    return this->WU(Us,nullptr,nullptr);
}

arma::mat 
trame::transfers::ntu::WU(const arma::mat& Us, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return Us - alpha(x_ind,y_ind);
}

arma::mat 
trame::transfers::ntu::WV(const arma::mat& Vs)
const
{
    return this->WV(Vs,nullptr,nullptr);
}

arma::mat 
trame::transfers::ntu::WV(const arma::mat& Vs, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return gamma(x_ind,y_ind) - Vs;
}

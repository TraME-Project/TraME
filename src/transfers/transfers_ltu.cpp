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
 * Linearly Transferable Utility (LTU) transfers class
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 07/25/2017
 */

#include "trame.hpp"

void
trame::transfers::ltu::build(const arma::mat& lambda_inp, const arma::mat& phi_inp, const bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = lambda_inp.n_rows;
    nbY = lambda_inp.n_cols;
    dim_params = 2*nbX*nbY;

    lambda = lambda_inp;
    phi = phi_inp;

    aux_zeta = 1 - lambda_inp;
    aux_phi_exp = arma::exp(phi_inp);
}

void
trame::transfers::ltu::trans()
{
    std::swap(nbX,nbY);
    //
    arma::inplace_trans(lambda);
    arma::inplace_trans(aux_zeta);
    lambda.swap(aux_zeta);

    arma::inplace_trans(phi);
    arma::inplace_trans(aux_phi_exp);
}

void
trame::transfers::ltu::gen_mmf(mmfs::cd& mmf_obj)
const
{
    mmf_obj.build(lambda,phi,need_norm);
}

trame::mmfs::cd
trame::transfers::ltu::gen_mmf()
const
{
    mmfs::cd mmf_obj;
    this->gen_mmf(mmf_obj);
    //
    return mmf_obj;
}

//
// DSE-related functions

arma::mat 
trame::transfers::ltu::Psi(const arma::mat& U, const arma::mat& V)
const
{
    return this->Psi(U,V,nullptr,nullptr);
}

// Implicit Parameterization
arma::mat 
trame::transfers::ltu::Psi(const arma::mat& U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = elem_prod(lambda(x_ind,y_ind), U) + elem_prod(aux_zeta(x_ind,y_ind), V) - phi(x_ind,y_ind);
    //
    return ret;
}

arma::mat 
trame::transfers::ltu::Psi(const double U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = lambda(x_ind,y_ind)*U + elem_prod(aux_zeta(x_ind,y_ind), V) - phi(x_ind,y_ind);
    //
    return ret;
}

arma::mat 
trame::transfers::ltu::Psi(const arma::mat& U, const double V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = elem_prod(lambda(x_ind,y_ind), U) + aux_zeta(x_ind,y_ind)*V - phi(x_ind,y_ind);
    //
    return ret;
}

double 
trame::transfers::ltu::Psi(const double U, const double V, const int x_ind, const int y_ind)
const
{
    return lambda(x_ind,y_ind) * U + aux_zeta(x_ind,y_ind)*V - phi(x_ind,y_ind);
}

arma::mat 
trame::transfers::ltu::du_Psi(const arma::mat& U, const arma::mat& V)
const
{
    return this->du_Psi(U,V,nullptr,nullptr);
}

arma::mat 
trame::transfers::ltu::du_Psi(const arma::mat& U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return lambda(x_ind,y_ind);
}

arma::mat 
trame::transfers::ltu::du_Psi(const double U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return lambda(x_ind,y_ind);
}

arma::mat 
trame::transfers::ltu::du_Psi(const arma::mat& U, const double V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return lambda(x_ind,y_ind);
}

arma::mat 
trame::transfers::ltu::dparams_Psi(const arma::mat& U, const arma::mat& V, const arma::mat& dparams)
const
{
    return this->dparams_Psi(U,V,&dparams);
}

arma::mat 
trame::transfers::ltu::dparams_Psi(const arma::mat& U, const arma::mat& V, const arma::mat* dparams)
const
{
    arma::mat ret(nbX,nbY);
    //
    const arma::vec U_minus_V = arma::vectorise(U - V);

    if (!dparams) {
        ret = arma::join_rows(arma::diagmat(U_minus_V),-arma::eye(nbX*nbY,nbX*nbY));
    } else {
        const arma::mat dparams_1 = dparams->rows(0,nbX*nbY-1);
        const arma::mat dparams_2 = dparams->rows(nbX*nbY,2*nbX*nbY-1);
        //
        ret = arma::vectorise(U_minus_V % dparams_1 - dparams_2);
    }
    //
    return ret;
}

// Explicit Parameterization
arma::mat 
trame::transfers::ltu::Ucal(const arma::mat& vs)
const
{
    return this->Ucal(vs,nullptr,nullptr);
}

arma::mat 
trame::transfers::ltu::Ucal(const arma::mat& vs, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    const arma::mat term_1 = phi(x_ind,y_ind) - elem_prod(aux_zeta(x_ind,y_ind), vs.t());
    const arma::mat term_2 = lambda(x_ind,y_ind);
        
    return term_1 / term_2;
}

double 
trame::transfers::ltu::Ucal(const double vs, const int xs, const int ys)
const
{
    const double term_1 = phi(xs,ys) - aux_zeta(xs,ys) * vs;
    const double term_2 = lambda(xs,ys);

    return term_1 / term_2;
}

arma::mat 
trame::transfers::ltu::Vcal(const arma::mat& us)
const
{
    return this->Vcal(us,nullptr,nullptr);
}

arma::mat 
trame::transfers::ltu::Vcal(const arma::mat& us, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    const arma::mat term_1 = phi(x_ind,y_ind) - elem_prod(lambda(x_ind,y_ind), us);
    const arma::mat term_2 = aux_zeta(x_ind,y_ind);

    return term_1 / term_2;
}

double 
trame::transfers::ltu::Vcal(const double us, const int xs, const int ys)
const
{
    const double term_1 = phi(xs,ys) - lambda(xs,ys) * us;
    const double term_2 = aux_zeta(xs,ys);

    return term_1 / term_2;
}

arma::mat 
trame::transfers::ltu::UW(const arma::mat& Ws)
const
{
    return this->UW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::ltu::UW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return - Psi(0.0,-Ws,&x_ind,&y_ind);
}

double 
trame::transfers::ltu::UW(const double Ws, const int x_ind, const int y_ind)
const
{
    return - Psi((double) 0.0,-Ws,x_ind,y_ind);
}

arma::mat 
trame::transfers::ltu::VW(const arma::mat& Ws)
const
{
    return this->VW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::ltu::VW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return - Psi(Ws,0.0,&x_ind,&y_ind);
}

double 
trame::transfers::ltu::VW(const double Ws, const int x_ind, const int y_ind)
const
{
    return - Psi(Ws,(double) 0.0,x_ind,y_ind);
}

arma::mat 
trame::transfers::ltu::dw_UW(const arma::mat& Ws)
const
{
    return this->dw_UW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::ltu::dw_UW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return 1.0 - du_Psi(0.0,-Ws,&x_ind,&y_ind);
}

arma::mat 
trame::transfers::ltu::dw_VW(const arma::mat& Ws)
const
{
    return this->dw_VW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::ltu::dw_VW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return - du_Psi(Ws,0.0,&x_ind,&y_ind);
}

arma::mat 
trame::transfers::ltu::WU(const arma::mat& Us)
const
{
    return this->WU(Us,nullptr,nullptr);
}

arma::mat 
trame::transfers::ltu::WU(const arma::mat& Us, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return (Us - phi(x_ind,y_ind)) / aux_zeta(x_ind,y_ind);
}

arma::mat 
trame::transfers::ltu::WV(const arma::mat& Vs)
const
{
    return this->WV(Vs,nullptr,nullptr);
}

arma::mat 
trame::transfers::ltu::WV(const arma::mat& Vs, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return (phi(x_ind,y_ind) - Vs) / lambda(x_ind,y_ind);
}

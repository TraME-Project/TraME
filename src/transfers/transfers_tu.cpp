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
 * Transferable Utility (TU) transfers class
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 07/25/2017
 */

#include "trame.hpp"

void 
trame::transfers::tu::build(const arma::mat& phi_inp, const bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = phi_inp.n_rows;
    nbY = phi_inp.n_cols;
    dim_params = nbX*nbY;

    phi = phi_inp;
    aux_phi_exp = arma::exp(phi_inp / 2.0);
}

void 
trame::transfers::tu::trans()
{
    std::swap(nbX,nbY);
    //
    arma::inplace_trans(phi);
    arma::inplace_trans(aux_phi_exp);
}

void
trame::transfers::tu::gen_mmf(mmfs::geo& mmf_obj)
const
{
    mmf_obj.build(phi,need_norm);
}

trame::mmfs::geo
trame::transfers::tu::gen_mmf()
const
{
    mmfs::geo mmf_obj;
    this->gen_mmf(mmf_obj);
    //
    return mmf_obj;
}

//
// DSE-related functions
//

//
// Implicit Parameterization

arma::mat 
trame::transfers::tu::Psi(const arma::mat& U, const arma::mat& V)
const
{
    return this->Psi(U,V,nullptr,nullptr);
}

arma::mat 
trame::transfers::tu::Psi(const arma::mat& U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return (U + V - phi(x_ind,y_ind)) / 2.0;
}

arma::mat 
trame::transfers::tu::Psi(const double U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return (U + V - phi(x_ind,y_ind)) / 2.0;
}

arma::mat 
trame::transfers::tu::Psi(const arma::mat& U, const double V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return (U + V - phi(x_ind,y_ind)) / 2.0;
}

double 
trame::transfers::tu::Psi(const double U, const double V, const int x_ind, const int y_ind)
const
{
    return (U + V - phi(x_ind,y_ind)) / 2.0;
}

// Derivative of Psi wrt u

arma::mat 
trame::transfers::tu::du_Psi(const arma::mat& U, const arma::mat& V)
const
{
    return this->du_Psi(U,V,nullptr,nullptr);
}

arma::mat 
trame::transfers::tu::du_Psi(const arma::mat& U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    ret.fill(0.5);
    //
    return ret;
}

arma::mat 
trame::transfers::tu::du_Psi(const double U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    ret.fill(0.5);
    //
    return ret;
}

arma::mat 
trame::transfers::tu::du_Psi(const arma::mat& U, const double V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    ret.fill(0.5);
    //
    return ret;
}

// dparams

arma::mat 
trame::transfers::tu::dparams_Psi(const arma::mat& U, const arma::mat& V)
const
{
    return this->dparams_Psi(U,V,nullptr);
}

arma::mat 
trame::transfers::tu::dparams_Psi(const arma::mat& U, const arma::mat& V, const arma::mat* dparams)
const
{
    arma::mat ret;
    //
    if (!dparams) {
        ret = - 0.5*arma::eye(nbX*nbY,nbX*nbY);
    } else {
        ret = - (*dparams)/2;
    }
    //
    return ret;
}

//
// Explicit Parameterization

// Ucal and Vcal

arma::mat 
trame::transfers::tu::Ucal(const arma::mat& vs)
const
{
    return this->Ucal(vs,nullptr,nullptr);
}

arma::mat 
trame::transfers::tu::Ucal(const arma::mat& vs, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return phi(x_ind,y_ind) - byrow(vs,x_ind.n_elem,y_ind.n_elem);
}

double 
trame::transfers::tu::Ucal(const double vs, const int xs, const int ys)
const
{
    return phi(xs,ys) - vs;
}

arma::mat 
trame::transfers::tu::Vcal(const arma::mat& us)
const
{
    return this->Vcal(us,nullptr,nullptr);
}

arma::mat 
trame::transfers::tu::Vcal(const arma::mat& us, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return phi(x_ind,y_ind) - us;
}

double 
trame::transfers::tu::Vcal(const double us, const int xs, const int ys)
const
{
    return phi(xs,ys) - us;
}

// UW and VW

arma::mat 
trame::transfers::tu::UW(const arma::mat& Ws)
const
{
    return this->UW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::tu::UW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return - Psi(0.0,-Ws,&x_ind,&y_ind);
}

double 
trame::transfers::tu::UW(const double Ws, const int x_ind, const int y_ind)
const
{
    return - Psi((double) 0.0,-Ws,x_ind,y_ind);
}

arma::mat 
trame::transfers::tu::VW(const arma::mat& Ws)
const
{
    return this->VW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::tu::VW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return - Psi(Ws,0.0,&x_ind,&y_ind);
}

double 
trame::transfers::tu::VW(const double Ws, const int x_ind, const int y_ind)
const
{
    return - Psi(Ws,(double) 0.0,x_ind,y_ind);
}

// dw

arma::mat 
trame::transfers::tu::dw_UW(const arma::mat& Ws)
const
{
    return this->dw_UW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::tu::dw_UW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return 1.0 - du_Psi(0.0,-Ws,&x_ind,&y_ind);
}

arma::mat 
trame::transfers::tu::dw_VW(const arma::mat& Ws)
const
{
    return this->dw_VW(Ws,nullptr,nullptr);
}

arma::mat 
trame::transfers::tu::dw_VW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return - du_Psi(Ws,0.0,&x_ind,&y_ind);
}

// WU and WV

arma::mat 
trame::transfers::tu::WU(const arma::mat& Us)
const
{
    return this->WU(Us,nullptr,nullptr);
}

arma::mat 
trame::transfers::tu::WU(const arma::mat& Us, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return 2*Us - phi(x_ind,y_ind);
}

arma::mat 
trame::transfers::tu::WV(const arma::mat& Vs)
const
{
    return this->WV(Vs,nullptr,nullptr);
}

arma::mat 
trame::transfers::tu::WV(const arma::mat& Vs, const arma::uvec* xs, const arma::uvec* ys)
const
{
    const arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    const arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    return phi(x_ind,y_ind) - 2*Vs;
}

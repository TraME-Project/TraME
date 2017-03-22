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
 * Linearly Transferable Utility (LTU) transfers class
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 03/12/2017
 */

#include "trame.hpp"

void trame::transfers::ltu::build(const arma::mat& lambda_LTU, const arma::mat& phi_LTU, bool need_norm_LTU)
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

void trame::transfers::ltu::trans()
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
// DSE-related functions

arma::mat trame::transfers::ltu::Psi(const arma::mat& U, const arma::mat& V)
const
{
    arma::mat ret = this->Psi(U,V,NULL,NULL);
    //
    return ret;
}

// Implicit Parameterization
arma::mat trame::transfers::ltu::Psi(const arma::mat& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = elem_prod(lambda(x_ind,y_ind), U) + elem_prod(aux_zeta(x_ind,y_ind), V) - phi(x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::ltu::Psi(const double& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = lambda(x_ind,y_ind)*U + elem_prod(aux_zeta(x_ind,y_ind), V) - phi(x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::ltu::Psi(const arma::mat& U, const double& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = elem_prod(lambda(x_ind,y_ind), U) + aux_zeta(x_ind,y_ind)*V - phi(x_ind,y_ind);
    //
    return ret;
}

double trame::transfers::ltu::Psi(const double& U, const double& V, int x_ind, int y_ind)
const
{
    double ret = lambda(x_ind,y_ind) * U + aux_zeta(x_ind,y_ind)*V - phi(x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::ltu::du_Psi(const arma::mat& U, const arma::mat& V)
const
{
    arma::mat ret = this->du_Psi(U,V,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::ltu::du_Psi(const arma::mat& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = lambda(x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::ltu::du_Psi(const double& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = lambda(x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::ltu::du_Psi(const arma::mat& U, const double& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = lambda(x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::ltu::dtheta_Psi(const arma::mat& U, const arma::mat& V, const arma::mat& dtheta)
{
    arma::mat ret = this->dtheta_Psi(U,V,&dtheta);
    return ret;
}

arma::mat trame::transfers::ltu::dtheta_Psi(const arma::mat& U, const arma::mat& V, const arma::mat* dtheta)
{
    arma::mat ret(nbX,nbY);
    //
    arma::vec U_minus_V = arma::vectorise(U - V);

    if (!dtheta) {
        ret = arma::join_rows(arma::diagmat(U_minus_V),-arma::eye(nbX*nbY,nbX*nbY));
    } else {
        arma::mat dtheta_1 = dtheta->rows(0,nbX*nbY-1);
        arma::mat dtheta_2 = dtheta->rows(nbX*nbY,2*nbX*nbY-1);
        //
        ret = arma::vectorise(U_minus_V % dtheta_1 - dtheta_2);
    }
    //
    return ret;
}

// Explicit Parameterization
arma::mat trame::transfers::ltu::Ucal(const arma::mat& vs)
const
{
    arma::mat ret = this->Ucal(vs,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::ltu::Ucal(const arma::mat& vs, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = phi(x_ind,y_ind) - elem_prod(aux_zeta(x_ind,y_ind), vs.t());
    arma::mat term_2 = lambda(x_ind,y_ind);
        
    arma::mat ret = term_1 / term_2;
    //
    return ret;
}

double trame::transfers::ltu::Ucal(const double& vs, int xs, int ys)
const
{
    double term_1 = phi(xs,ys) - aux_zeta(xs,ys) * vs;
    double term_2 = lambda(xs,ys);

    double ret = term_1 / term_2;
    //
    return ret;
}

arma::mat trame::transfers::ltu::Vcal(const arma::mat& us)
const
{
    arma::mat ret = this->Vcal(us,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::ltu::Vcal(const arma::mat& us, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = phi(x_ind,y_ind) - elem_prod(lambda(x_ind,y_ind), us);
    arma::mat term_2 = aux_zeta(x_ind,y_ind);

    arma::mat ret = term_1 / term_2;
    //
    return ret;
}

double trame::transfers::ltu::Vcal(const double& us, int xs, int ys)
const
{
    double term_1 = phi(xs,ys) - lambda(xs,ys) * us;
    double term_2 = aux_zeta(xs,ys);

    double ret = term_1 / term_2;
    //
    return ret;
}

arma::mat trame::transfers::ltu::UW(const arma::mat& Ws)
const
{
    arma::mat ret = this->UW(Ws,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::ltu::UW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - Psi(0.0,-Ws,&x_ind,&y_ind);
    //
    return ret;
}

double trame::transfers::ltu::UW(const double& Ws, int x_ind, int y_ind)
const
{
    double ret = - Psi((double) 0.0,-Ws,x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::ltu::VW(const arma::mat& Ws)
const
{
    arma::mat ret = this->VW(Ws,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::ltu::VW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - Psi(Ws,0.0,&x_ind,&y_ind);
    //
    return ret;
}

double trame::transfers::ltu::VW(const double& Ws, int x_ind, int y_ind)
const
{
    double ret = - Psi(Ws,(double) 0.0,x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::ltu::dw_UW(const arma::mat& Ws)
const
{
    arma::mat ret = this->dw_UW(Ws,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::ltu::dw_UW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = 1.0 - du_Psi(0.0,-Ws,&x_ind,&y_ind);
    //
    return ret;
}

arma::mat trame::transfers::ltu::dw_VW(const arma::mat& Ws)
const
{
    arma::mat ret = this->dw_VW(Ws,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::ltu::dw_VW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - du_Psi(Ws,0.0,&x_ind,&y_ind);

    return ret;
}

arma::mat trame::transfers::ltu::WU(const arma::mat& Us)
{
    arma::mat ret = this->WU(Us,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::ltu::WU(const arma::mat& Us, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = (Us - phi(x_ind,y_ind)) / aux_zeta(x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::ltu::WV(const arma::mat& Vs)
{
    arma::mat ret = this->WV(Vs,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::ltu::WV(const arma::mat& Vs, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = (phi(x_ind,y_ind) - Vs) / lambda(x_ind,y_ind);
    //
    return ret;
}

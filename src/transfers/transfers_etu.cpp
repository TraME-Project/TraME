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
 * ETU transfers class
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 03/12/2017
 */

#include "trame.hpp"

void trame::etu::build(const arma::mat& alpha_ETU, const arma::mat& gamma_ETU, const arma::mat& tau_ETU, bool need_norm_ETU)
{
    need_norm = need_norm_ETU;

    nbX = alpha_ETU.n_rows;
    nbY = alpha_ETU.n_cols;
    nbParams = 3*nbX*nbY;

    alpha = alpha_ETU;
    gamma = gamma_ETU;
    tau   = tau_ETU;
    kappa = - 1.0 / tau_ETU;

    aux_alpha = - elem_div(alpha_ETU, tau_ETU);
    aux_gamma = - elem_div(gamma_ETU, tau_ETU);

    aux_alpha_exp = arma::exp(- elem_div(alpha_ETU, tau_ETU));
    aux_gamma_exp = arma::exp(- elem_div(gamma_ETU, tau_ETU));
}

void trame::etu::trans()
{
    int nbX_temp = nbX;

    nbX = nbY;
    nbY = nbX_temp;
    //
    arma::mat alpha_temp = alpha;

    alpha = gamma.t();
    gamma = alpha_temp.t();
    tau = tau.t();
    kappa = kappa.t();
    //
    arma::mat aux_alpha_temp = aux_alpha;
    arma::mat aux_alpha_exp_temp = aux_alpha_exp;

    arma::mat aux_alpha = aux_gamma.t();
    arma::mat aux_gamma = aux_alpha_temp.t();

    arma::mat aux_alpha_exp = aux_gamma_exp.t();
    arma::mat aux_gamma_exp = aux_alpha_exp_temp.t();
    //
}

//
// MFE-related functions

arma::mat trame::etu::M(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    arma::mat ret = this->M(a_xs,b_ys,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::M(const arma::mat& a_xs, const arma::mat& b_ys, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = arma::exp(aux_alpha(x_ind,y_ind) + elem_prod(kappa(x_ind,y_ind), arma::log(a_xs)));
    arma::mat term_2 = arma::exp(aux_gamma(x_ind,y_ind) + arma::trans(elem_prod(arma::trans(kappa(x_ind,y_ind)), arma::log(b_ys))));

    arma::mat ret = arma::exp( - tau(x_ind,y_ind) % arma::log((term_1 + term_2)/2) ); 
    //
    return ret;
}

arma::mat trame::etu::M(const double& a_xs, const arma::mat& b_ys, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, (int) nbX-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, (int) nbY-1);
    //
    arma::mat term_1 = arma::exp(aux_alpha(x_ind,y_ind) + kappa(x_ind,y_ind) * std::log(a_xs));
    arma::mat term_2 = arma::exp(aux_gamma(x_ind,y_ind) + arma::trans(elem_prod(arma::trans(kappa(x_ind,y_ind)), arma::log(b_ys))));

    arma::mat ret = arma::exp( - tau(x_ind,y_ind) % arma::log((term_1 + term_2)/2) ); 
    //
    return ret;
}

arma::mat trame::etu::M(const arma::mat& a_xs, const double& b_ys, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, (int) nbX-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, (int) nbY-1);
    //
    arma::mat term_1 = arma::exp(aux_alpha(x_ind,y_ind) + elem_prod(kappa(x_ind,y_ind), arma::log(a_xs)));
    arma::mat term_2 = arma::exp(aux_gamma(x_ind,y_ind) + kappa(x_ind,y_ind) * std::log(b_ys));

    arma::mat ret = arma::exp( - tau(x_ind,y_ind) % arma::log((term_1 + term_2)/2) );
    //
    return ret;
}

arma::mat trame::etu::Mx0(const arma::mat& a_x)
const
{
    return a_x;
}

arma::mat trame::etu::M0y(const arma::mat& b_y)
const
{
    return b_y;
}

//
// DSE-related functions

arma::mat trame::etu::Psi(const arma::mat& U, const arma::mat& V)
const
{
    arma::mat ret = this->Psi(U,V,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::Psi(const arma::mat& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat temp_1 = elem_prod(arma::exp(elem_div(U, tau(x_ind,y_ind))), aux_alpha_exp(x_ind,y_ind));
    arma::mat temp_2 = elem_prod(arma::exp(elem_div(V, tau(x_ind,y_ind))), aux_gamma_exp(x_ind,y_ind));

    arma::mat ret =  elem_prod(tau(x_ind,y_ind), arma::log(0.5 * (temp_1 + temp_2)));
    //
    return ret;
}

arma::mat trame::etu::Psi(const double& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat temp_1 = elem_prod(arma::exp(U/tau(x_ind,y_ind)), aux_alpha_exp(x_ind,y_ind));
    arma::mat temp_2 = elem_prod(arma::exp(elem_div(V, tau(x_ind,y_ind))), aux_gamma_exp(x_ind,y_ind));
        
    arma::mat ret =  elem_prod(tau(x_ind,y_ind), arma::log(0.5 * (temp_1 + temp_2)));
    //
    return ret;
}

arma::mat trame::etu::Psi(const arma::mat& U, const double& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat temp_1 = elem_prod(arma::exp(elem_div(U, tau(x_ind,y_ind))), aux_alpha_exp(x_ind,y_ind));
    arma::mat temp_2 = elem_prod(arma::exp(V/tau(x_ind,y_ind)), aux_gamma_exp(x_ind,y_ind));

    arma::mat ret =  elem_prod(tau(x_ind,y_ind), arma::log(0.5 * (temp_1 + temp_2)));
    //
    return ret;
}

double trame::etu::Psi(const double& U, const double& V, int x_ind, int y_ind)
const
{
    double temp_1 = std::exp(U/tau(x_ind,y_ind)) * aux_alpha_exp(x_ind,y_ind);
    double temp_2 = std::exp(V/tau(x_ind,y_ind)) * aux_gamma_exp(x_ind,y_ind);

    double ret =  tau(x_ind,y_ind) * std::log(0.5 * (temp_1 + temp_2));
    //
    return ret;
}

arma::mat trame::etu::du_Psi(const arma::mat& U, const arma::mat& V)
const
{
    arma::mat ret = this->du_Psi(U,V,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::du_Psi(const arma::mat& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret =  1 / (1 + arma::exp((V - U + alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind)));
    //
    return ret;
}

arma::mat trame::etu::du_Psi(const double& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret =  1 / (1 + arma::exp((V - U + alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind)));
    //
    return ret;
}

arma::mat trame::etu::du_Psi(const arma::mat& U, const double& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret =  1 / (1 + arma::exp((V - U + alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind)));
    //
    return ret;
}

arma::mat trame::etu::dtheta_Psi(const arma::mat& U, const arma::mat& V, const arma::mat& dtheta)
{
    arma::mat ret = this->dtheta_Psi(U,V,&dtheta);
    return ret;
}

arma::mat trame::etu::dtheta_Psi(const arma::mat& U, const arma::mat& V, const arma::mat* dtheta)
{
    arma::mat ret(nbX,nbY);
    //
    arma::mat du_psi_mat = du_Psi(U,V,NULL,NULL);
    arma::vec du_psi = arma::vectorise(du_psi_mat);

    if (!dtheta) {
        arma::mat term_1 = (U - alpha) % du_psi_mat;
        arma::mat term_2 = (V - gamma) % (1 - du_psi_mat);

        arma::mat dsigma_psi_mat = (Psi(U,V,NULL,NULL) - term_1 - term_2)/tau;
        arma::vec dsigma_psi = arma::vectorise(dsigma_psi_mat);
        //
        ret = arma::join_rows(arma::diagmat(-du_psi),arma::join_rows(arma::diagmat(du_psi-1),arma::diagmat(dsigma_psi)));
    } else {
        arma::mat dtheta_1 = dtheta->rows(0,nbX*nbY-1);
        arma::mat dtheta_2 = dtheta->rows(nbX*nbY,2*nbX*nbY-1);
        arma::mat dtheta_3 = dtheta->rows(2*nbX*nbY,3*nbX*nbY-1);

        arma::vec dsigma_psi_dtheta = arma::zeros(dtheta_3.n_rows,1);
        double min_check = elem_min(dtheta_3);

        if (min_check!=0) {
            arma::mat term_1 = (U - alpha) % du_psi_mat;
            arma::mat term_2 = (V - gamma) % (1 - du_psi_mat);

            arma::mat dsigma_psi_mat = (Psi(U,V,NULL,NULL) - term_1 - term_2)/tau;

            dsigma_psi_dtheta = dtheta_3 % arma::vectorise(dsigma_psi_mat);
        }
        //
        ret = arma::vectorise(-du_psi % dtheta_1 - (1-du_psi) % dtheta_2 + dsigma_psi_dtheta);
    }
    //
    return ret;
}

arma::mat trame::etu::Ucal(const arma::mat& vs)
const
{
    arma::mat ret = this->Ucal(vs,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::Ucal(const arma::mat& vs, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = vs - gamma(x_ind,y_ind);
    arma::mat term_log = 2 - arma::exp(elem_div(term_1,tau(x_ind,y_ind)));

    arma::mat ret = alpha(x_ind,y_ind) + elem_prod(tau(x_ind,y_ind), arma::log(term_log));
    //
    return ret;
}

double trame::etu::Ucal(const double& vs, int xs, int ys)
const
{
    double term_1 = vs - gamma(xs,ys);
    double term_log = 2 - std::exp(term_1/tau(xs,ys));

    double ret = alpha(xs,ys) + tau(xs,ys) * std::log(term_log);
    //
    return ret;
}

arma::mat trame::etu::Vcal(const arma::mat& us)
const
{
    arma::mat ret = this->Vcal(us,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::Vcal(const arma::mat& us, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = us - alpha(x_ind,y_ind);
    arma::mat term_log = 2 - arma::exp(elem_div(term_1,tau(x_ind,y_ind)));

    arma::mat ret = gamma(x_ind,y_ind) + elem_prod(tau(x_ind,y_ind), arma::log(term_log));
    //
    return ret;
}

double trame::etu::Vcal(const double& us, int xs, int ys)
const
{
    double term_1 = us - alpha(xs,ys);
    double term_log = 2 - std::exp(term_1/tau(xs,ys));

    double ret = gamma(xs,ys) + tau(xs,ys) * std::log(term_log);
    //
    return ret;
}

arma::mat trame::etu::UW(const arma::mat& Ws)
const
{
    arma::mat ret = this->UW(Ws,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::UW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - Psi(0.0,-Ws,&x_ind,&y_ind);

    return ret;
}

double trame::etu::UW(const double& Ws, int x_ind, int y_ind)
const
{
    double ret = - Psi((double) 0.0,-Ws,x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::etu::VW(const arma::mat& Ws)
const
{
    arma::mat ret = this->VW(Ws,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::VW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - Psi(Ws,0.0,&x_ind,&y_ind);

    return ret;
}

double trame::etu::VW(const double& Ws, int x_ind, int y_ind)
const
{
    double ret = - Psi(Ws,(double) 0.0,x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::etu::dw_UW(const arma::mat& Ws)
const
{
    arma::mat ret = this->dw_UW(Ws,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::dw_UW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = 1.0 - du_Psi(0.0,-Ws,&x_ind,&y_ind);

    return ret;
}

arma::mat trame::etu::dw_VW(const arma::mat& Ws)
const
{
    arma::mat ret = this->dw_VW(Ws,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::dw_VW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - du_Psi(Ws,0.0,&x_ind,&y_ind);

    return ret;
}

arma::mat trame::etu::WU(const arma::mat& Us)
{
    arma::mat ret = this->WU(Us,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::WU(const arma::mat& Us, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = 2 * arma::exp( (gamma(x_ind,y_ind) - Us)/tau(x_ind,y_ind) );
    arma::mat term_2 = arma::exp( (gamma(x_ind,y_ind) - alpha(x_ind,y_ind))/tau(x_ind,y_ind) );
    arma::mat term_log = term_1 - term_2;

    arma::mat ret = - tau(x_ind,y_ind) % arma::log(term_log);
    //
    return ret;
}

arma::mat trame::etu::WV(const arma::mat& Vs)
{
    arma::mat ret = this->WV(Vs,NULL,NULL);
    //
    return ret;
}

arma::mat trame::etu::WV(const arma::mat& Vs, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = 2 * arma::exp( (alpha(x_ind,y_ind) - Vs)/tau(x_ind,y_ind) );
    arma::mat term_2 = arma::exp( (alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind) );
    arma::mat term_log = term_1 - term_2;

    ret = - tau(x_ind,y_ind) % arma::log(term_log);
    //
    return ret;
}

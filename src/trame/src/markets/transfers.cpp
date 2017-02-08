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
 * transfers class
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 11/01/2016
 */

#include "trame.hpp"

void trame::transfers::build_ETU(const arma::mat& alpha_ETU, const arma::mat& gamma_ETU, const arma::mat& tau_ETU)
{
    alpha = alpha_ETU;
    gamma = gamma_ETU;
    tau   = tau_ETU;

    nbX = alpha_ETU.n_rows;
    nbY = alpha_ETU.n_cols;
    nbParams = 3*nbX*nbY;

    aux_exp_alphaovertau = arma::exp(- elem_div(alpha_ETU, tau_ETU));
    aux_exp_gammaovertau = arma::exp(- elem_div(gamma_ETU, tau_ETU));

    ETU = true;
    transfers_type = 2;
}

void trame::transfers::build_LTU(const arma::mat& lambda_LTU, const arma::mat& phi_LTU)
{
    lambda = lambda_LTU;
    phi    = phi_LTU;

    nbX = lambda_LTU.n_rows;
    nbY = lambda_LTU.n_cols;
    nbParams = 2*nbX*nbY;

    aux_zeta = 1 - lambda_LTU;

    LTU = true;
    transfers_type = 1;
}

void trame::transfers::build_NTU(const arma::mat& alpha_NTU, const arma::mat& gamma_NTU)
{
    alpha = alpha_NTU;
    gamma = gamma_NTU;

    nbX = alpha_NTU.n_rows;
    nbY = alpha_NTU.n_cols;
    nbParams = 2*nbX*nbY;

    NTU = true;
    transfers_type = 2;
}

void trame::transfers::build_TU(const arma::mat& phi_TU)
{
    phi = phi_TU;

    nbX = phi_TU.n_rows;
    nbY = phi_TU.n_cols;
    nbParams = nbX*nbY;

    TU = true;
    transfers_type = 1;
}

void trame::transfers::trans()
{
    int nbX_temp;
    //
    if (ETU) {
        nbX_temp = nbX;
        nbX = nbY;
        nbY = nbX_temp;

        arma::mat alpha_temp = alpha;

        alpha = gamma.t();
        gamma = alpha_temp.t();
        tau   = tau.t();

        arma::mat aux_exp_alphaovertau_temp = aux_exp_alphaovertau;

        aux_exp_alphaovertau = aux_exp_gammaovertau.t();
        aux_exp_gammaovertau = aux_exp_alphaovertau_temp.t();
    }

    if (LTU) {
        nbX_temp = nbX;
        nbX = nbY;
        nbY = nbX_temp;

        arma::mat lambda_temp = lambda;

        lambda   = aux_zeta.t();
        phi      = phi.t();
        aux_zeta = lambda_temp.t();
    }

    if (NTU) {
        nbX_temp = nbX;
        nbX = nbY;
        nbY = nbX_temp;

        arma::mat alpha_temp = alpha;

        alpha = gamma.t();
        gamma = alpha_temp.t();
    }

    if (TU) {
        nbX_temp = nbX;
        nbX = nbY;
        nbY = nbX_temp;

        phi = phi.t();
    }
}

// Psi functuons are const because of calls to const market pointers in solvers
arma::mat trame::transfers::Psi(const arma::mat& U, const arma::mat& V)
const
{
    arma::mat ret = this->Psi(U,V,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::Psi(const arma::mat& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);

    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    //
    if (ETU) {
        arma::mat temp_1 = elem_prod(arma::exp(elem_div(U,tau(x_ind,y_ind))), aux_exp_alphaovertau(x_ind,y_ind));
        arma::mat temp_2 = elem_prod(arma::exp(elem_div(V,tau(x_ind,y_ind))), aux_exp_gammaovertau(x_ind,y_ind));

        ret =  elem_prod(tau(x_ind,y_ind), arma::log(0.5 * (temp_1 + temp_2)));
        goto finished;
    }

    if (LTU) {
        ret = elem_prod(lambda(x_ind,y_ind), U) + elem_prod(aux_zeta(x_ind,y_ind), V) - phi(x_ind,y_ind);
        goto finished;
    }

    if (NTU) {
        ret = arma::max(U - alpha(x_ind,y_ind), V - gamma(x_ind,y_ind));
        goto finished;
    }

    if (TU) {
        ret = (U + V - phi(x_ind,y_ind)) / 2;
        goto finished;
    }
    //
finished:
    return ret;
}

arma::mat trame::transfers::Psi(const double& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);

    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    //
    if (ETU) {
        arma::mat temp_1 = elem_prod(arma::exp(U/tau(x_ind,y_ind)), aux_exp_alphaovertau(x_ind,y_ind));
        arma::mat temp_2 = elem_prod(arma::exp(elem_div(V,tau(x_ind,y_ind))), aux_exp_gammaovertau(x_ind,y_ind));
        
        ret =  elem_prod(tau(x_ind,y_ind), arma::log(0.5 * (temp_1 + temp_2)));
        goto finished;
    }

    if (LTU) {
        ret = lambda(x_ind,y_ind)*U + elem_prod(aux_zeta(x_ind,y_ind), V) - phi(x_ind,y_ind);
        goto finished;
    }

    if (NTU) {
        ret = arma::max(U - alpha(x_ind,y_ind), V - gamma(x_ind,y_ind));
        goto finished;
    }

    if (TU) {
        ret = (U + V - phi(x_ind,y_ind)) / 2;
        goto finished;
    }
    //
finished:
    return ret;
}

arma::mat trame::transfers::Psi(const arma::mat& U, const double& V, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);

    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    //
    if (ETU) {
        arma::mat temp_1 = elem_prod(arma::exp(elem_div(U,tau(x_ind,y_ind))), aux_exp_alphaovertau(x_ind,y_ind));
        arma::mat temp_2 = elem_prod(arma::exp(V/tau(x_ind,y_ind)), aux_exp_gammaovertau(x_ind,y_ind));

        ret =  elem_prod(tau(x_ind,y_ind), arma::log(0.5 * (temp_1 + temp_2)));
        goto finished;
    }

    if (LTU) {
        ret = elem_prod(lambda(x_ind,y_ind), U) + aux_zeta(x_ind,y_ind)*V - phi(x_ind,y_ind);
        goto finished;
    }

    if (NTU) {
        ret = arma::max(U - alpha(x_ind,y_ind), V - gamma(x_ind,y_ind));
        goto finished;
    }

    if (TU) {
        ret = (U + V - phi(x_ind,y_ind)) / 2;
        goto finished;
    }
    //
finished:
    return ret;
}

double trame::transfers::Psi(const double& U, const double& V, int x_ind, int y_ind)
const
{
    double ret = 0.0;
    //
    if (ETU) {
        double temp_1 = std::exp(U/tau(x_ind,y_ind)) * aux_exp_alphaovertau(x_ind,y_ind);
        double temp_2 = std::exp(V/tau(x_ind,y_ind)) * aux_exp_gammaovertau(x_ind,y_ind);

        ret =  tau(x_ind,y_ind) * std::log(0.5 * (temp_1 + temp_2));
        goto finished;
    }

    if (LTU) {
        ret = lambda(x_ind,y_ind) * U + aux_zeta(x_ind,y_ind)*V - phi(x_ind,y_ind);
        goto finished;
    }

    if (NTU) {
        ret = std::max(U - alpha(x_ind,y_ind), V - gamma(x_ind,y_ind));
        goto finished;
    }

    if (TU) {
        ret = (U + V - phi(x_ind,y_ind)) / 2;
        goto finished;
    }
    //
finished:
    return ret;
}

arma::mat trame::transfers::du_Psi(const arma::mat& U, const arma::mat& V)
{
    arma::mat ret = this->du_Psi(U,V,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::du_Psi(const arma::mat& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);

    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    //
    if (ETU) {
        ret =  1 / (1 + arma::exp((V - U + alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind)));
        goto finished;
    }

    if (LTU) {
        ret = lambda(x_ind,y_ind);
        goto finished;
    }

    if (NTU) {
        ret.elem( arma::find(U - alpha(x_ind,y_ind) >= V - gamma(x_ind,y_ind)) ).ones();
        goto finished;
    }

    if (TU) {
        ret.fill(0.5);
        goto finished;
    }
    //
finished:
    return ret;
}

arma::mat trame::transfers::du_Psi(const double& U, const arma::mat& V, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);

    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    //
    if (ETU) {
        ret =  1 / (1 + arma::exp((V - U + alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind)));
        goto finished;
    }

    if (LTU) {
        ret = lambda(x_ind,y_ind);
        goto finished;
    }

    if (NTU) {
        ret.elem( arma::find(U - alpha(x_ind,y_ind) >= V - gamma(x_ind,y_ind)) ).ones();
        goto finished;
    }

    if (TU) {
        ret.fill(0.5);
        goto finished;
    }
    //
finished:
    return ret;
}

arma::mat trame::transfers::du_Psi(const arma::mat& U, const double& V, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);

    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    //
    if (ETU) {
        ret =  1 / (1 + arma::exp((V - U + alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind)));
        goto finished;
    }

    if (LTU) {
        ret = lambda(x_ind,y_ind);
        goto finished;
    }

    if (NTU) {
        ret.elem( arma::find(U - alpha(x_ind,y_ind) >= V - gamma(x_ind,y_ind)) ).ones();
        goto finished;
    }

    if (TU) {
        ret.fill(0.5);
        goto finished;
    }
    //
finished:
    return ret;
}

arma::mat trame::transfers::dtheta_Psi(const arma::mat& U, const arma::mat& V, arma::mat* dtheta)
{
    arma::mat ret(nbX,nbY);
    //
    if (ETU) {
        arma::mat du_psi_mat = du_Psi(U,V,NULL,NULL);
        arma::vec du_psi = arma::vectorise(du_psi_mat);

        if (!dtheta) {
            arma::mat term_1, term_2;
            term_1 = (U - alpha) % du_psi_mat;
            term_2 = (V - gamma) % (1 - du_psi_mat);

            arma::mat dsigmapsi_mat = (Psi(U,V,NULL,NULL) - term_1 - term_2)/tau;
            arma::vec dsigmapsi = arma::vectorise(dsigmapsi_mat);
            //
            ret = arma::join_rows(arma::diagmat(-du_psi),arma::join_rows(arma::diagmat(du_psi-1),arma::diagmat(dsigmapsi)));
            goto finished;
        } else {
            arma::mat dtheta_1 = dtheta->rows(0,nbX*nbY-1);
            arma::mat dtheta_2 = dtheta->rows(nbX*nbY,2*nbX*nbY-1);
            arma::mat dtheta_3 = dtheta->rows(2*nbX*nbY,3*nbX*nbY-1);

            arma::vec dsigmapsidtheta = arma::zeros(dtheta_3.n_rows,1);
            double min_check = arma::as_scalar(arma::min(arma::min(dtheta_3)));

            if (min_check!=0) {
                arma::mat term_1, term_2;
                term_1 = (U - alpha) % du_psi_mat;
                term_2 = (V - gamma) % (1 - du_psi_mat);

                arma::mat dsigmapsi_mat = (Psi(U,V,NULL,NULL) - term_1 - term_2)/tau;

                dsigmapsidtheta = dtheta_3 % arma::vectorise(dsigmapsi_mat);
            }
            //
            ret = arma::vectorise(-du_psi % dtheta_1 - (1-du_psi) % dtheta_2 + dsigmapsidtheta);
            goto finished;
        }
    }

    if (LTU) {
        arma::vec U_minus_V = arma::vectorise(U - V);

        if (!dtheta) {
            ret = arma::join_rows(arma::diagmat(U_minus_V),-arma::eye(nbX*nbY,nbX*nbY));
            goto finished;
        } else {
            arma::mat dtheta_1 = dtheta->rows(0,nbX*nbY-1);
            arma::mat dtheta_2 = dtheta->rows(nbX*nbY,2*nbX*nbY-1);
            //
            ret = arma::vectorise(U_minus_V % dtheta_1 - dtheta_2);
            goto finished;
        }
    }

    if (NTU) {
        arma::vec du_psi = arma::vectorise(du_Psi(U,V,NULL,NULL));

        if (!dtheta) {
            ret = - arma::join_rows(arma::diagmat(du_psi),arma::diagmat(1 - du_psi));
            goto finished;
        } else {
            arma::mat dtheta_1 = dtheta->rows(0,nbX*nbY-1);
            arma::mat dtheta_2 = dtheta->rows(nbX*nbY,2*nbX*nbY-1);

            ret = - arma::vectorise(du_psi % dtheta_1 + (1 - du_psi) % dtheta_2);
            goto finished;
        }
    }

    if (TU) {
        if (!dtheta) {
            ret = - 0.5*arma::eye(nbX*nbY,nbX*nbY);
            goto finished;
        } else {
            ret = - (*dtheta)/2;
            goto finished;
        }
    }
    //
finished:
    return ret;
}

arma::mat trame::transfers::Ucal(const arma::mat& vs)
const
{
    arma::mat ret = this->Ucal(vs,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::Ucal(const arma::mat& vs, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);

    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    //
    if (ETU) {
        arma::mat term_1 = vs - gamma(x_ind,y_ind);
        arma::mat term_log = 2 - arma::exp(elem_div(term_1,tau(x_ind,y_ind)));

        ret = alpha(x_ind,y_ind) + elem_prod(tau(x_ind,y_ind), arma::log(term_log));
        goto finished;
    }

    if (LTU) {
        arma::mat term_1 = phi(x_ind,y_ind) - elem_prod(aux_zeta(x_ind,y_ind), vs.t());
        arma::mat term_2 = lambda(x_ind,y_ind);
        
        ret = term_1 / term_2;
        goto finished;
    }
    /* no Ucal for NTU?
    if (NTU) {

    }
    */
    if (TU) {
        ret = phi(x_ind,y_ind) - vs;
        goto finished;
    }
    //
finished:
    return ret;
}

double trame::transfers::Ucal(const double& vs, int xs, int ys)
const
{
    double ret = 0;
    //
    if (ETU) {
        double term_1 = vs - gamma(xs,ys);
        double term_log = 2 - std::exp(term_1/tau(xs,ys));

        ret = alpha(xs,ys) + tau(xs,ys) * std::log(term_log);
        goto finished;
    }

    if (LTU) {
        double term_1 = phi(xs,ys) - aux_zeta(xs,ys) * vs;
        double term_2 = lambda(xs,ys);

        ret = term_1 / term_2;
        goto finished;
    }
    /* no Ucal for NTU?
    if (NTU) {

    }
    */
    if (TU) {
        ret = phi(xs,ys) - vs;
        goto finished;
    }
    //
finished:
    return ret;
}

arma::mat trame::transfers::Vcal(const arma::mat& us)
const
{
    arma::mat ret = this->Vcal(us,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::Vcal(const arma::mat& us, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);

    arma::mat ret(x_ind.n_elem,y_ind.n_elem);
    //
    if (ETU) {
        arma::mat term_1 = us - alpha(x_ind,y_ind);
        arma::mat term_log = 2 - arma::exp(elem_div(term_1,tau(x_ind,y_ind)));

        ret = gamma(x_ind,y_ind) + elem_prod(tau(x_ind,y_ind), arma::log(term_log));
        goto finished;
    }

    if (LTU) {
        arma::mat term_1 = phi(x_ind,y_ind) - elem_prod(lambda(x_ind,y_ind), us);
        arma::mat term_2 = aux_zeta(x_ind,y_ind);

        ret = term_1 / term_2;
        goto finished;
    }
    /* no Ucal for NTU?
    if (NTU) {

    }
    */
    if (TU) {
        ret = phi(x_ind,y_ind) - us;
        goto finished;
    }
    //
finished:
    return ret;
}

double trame::transfers::Vcal(const double& us, int xs, int ys)
const
{
    double ret = 0;
    //
    if (ETU) {
        double term_1 = us - alpha(xs,ys);
        double term_log = 2 - std::exp(term_1/tau(xs,ys));

        ret = gamma(xs,ys) + tau(xs,ys) * std::log(term_log);
        goto finished;
    }

    if (LTU) {
        double term_1 = phi(xs,ys) - lambda(xs,ys) * us;
        double term_2 = aux_zeta(xs,ys);

        ret = term_1 / term_2;
        goto finished;
    }
    /* no Ucal for NTU?
    if (NTU) {

    }
    */
    if (TU) {
        ret = phi(xs,ys) - us;
        goto finished;
    }
    //
finished:
    return ret;
}

arma::mat trame::transfers::UW(const arma::mat& Ws)
const
{
    arma::mat ret = this->UW(Ws,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::UW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - Psi(0.0,-Ws,&x_ind,&y_ind);

    return ret;
}

double trame::transfers::UW(const double& Ws, int x_ind, int y_ind)
const
{
    double ret = - Psi((double) 0.0,-Ws,x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::VW(const arma::mat& Ws)
const
{
    arma::mat ret = this->VW(Ws,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::VW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - Psi(Ws,0.0,&x_ind,&y_ind);

    return ret;
}

double trame::transfers::VW(const double& Ws, int x_ind, int y_ind)
const
{
    double ret = - Psi(Ws,(double) 0.0,x_ind,y_ind);
    //
    return ret;
}

arma::mat trame::transfers::du_UW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = 1 - du_Psi(0.0,-Ws,&x_ind,&y_ind);

    return ret;
}

arma::mat trame::transfers::du_VW(const arma::mat& Ws, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - du_Psi(Ws,0.0,&x_ind,&y_ind);

    return ret;
}

arma::mat trame::transfers::WU(const arma::mat& Us)
{
    arma::mat ret = this->WU(Us,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::WU(const arma::mat& Us, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);

    arma::mat ret;
    //
    if (ETU) {
        arma::mat term_1 = 2 * arma::exp( (gamma(x_ind,y_ind) - Us)/tau(x_ind,y_ind) );
        arma::mat term_2 = arma::exp( (gamma(x_ind,y_ind) - alpha(x_ind,y_ind))/tau(x_ind,y_ind) );
        arma::mat term_log = term_1 - term_2;

        ret = - tau(x_ind,y_ind) % arma::log(term_log);
        goto finished;
    }

    if (LTU) {
        ret = (Us - phi(x_ind,y_ind)) / aux_zeta(x_ind,y_ind);
        goto finished;
    }

    if (NTU) {
        ret = Us - alpha(x_ind,y_ind);
        goto finished;
    }

    if (TU) {
        ret = 2*Us - phi(x_ind,y_ind);
        goto finished;
    }
    //
finished:
    return ret;
}

arma::mat trame::transfers::WV(const arma::mat& Vs)
{
    arma::mat ret = this->WV(Vs,NULL,NULL);
    //
    return ret;
}

arma::mat trame::transfers::WV(const arma::mat& Vs, arma::uvec* xs, arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    
    arma::mat ret;
    //
    if (ETU) {
        arma::mat term_1 = 2 * arma::exp( (alpha(x_ind,y_ind) - Vs)/tau(x_ind,y_ind) );
        arma::mat term_2 = arma::exp( (alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind) );
        arma::mat term_log = term_1 - term_2;

        ret = - tau(x_ind,y_ind) % arma::log(term_log);
        goto finished;
    }

    if (LTU) {
        ret = (phi(x_ind,y_ind) - Vs) / lambda(x_ind,y_ind);
        goto finished;
    }

    if (NTU) {
        ret = gamma(x_ind,y_ind) - Vs;
        goto finished;
    }

    if (TU) {
        ret = phi(x_ind,y_ind) - 2*Vs;
        goto finished;
    }
    //
finished:
    return ret;
}

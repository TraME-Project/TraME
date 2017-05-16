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
 * Exponentially Transferable Utility (ETU) transfers class
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 03/12/2017
 */

#include "trame.hpp"

void 
trame::transfers::etu::build(const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = alpha_inp.n_rows;
    nbY = alpha_inp.n_cols;
    nbParams = 3*nbX*nbY;

    alpha = alpha_inp;
    gamma = gamma_inp;
    tau   = tau_inp;

    kappa = - 1.0 / tau_inp;

    aux_alpha = - elem_div(alpha_inp, tau_inp);
    aux_gamma = - elem_div(gamma_inp, tau_inp);

    aux_alpha_exp = arma::exp(- elem_div(alpha_inp, tau_inp));
    aux_gamma_exp = arma::exp(- elem_div(gamma_inp, tau_inp));
}

void 
trame::transfers::etu::trans()
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

void
trame::transfers::etu::gen_mmf(mmfs::ces& mmf_obj)
const
{
    mmf_obj.build(alpha,gamma,tau,need_norm);
}

trame::mmfs::ces
trame::transfers::etu::gen_mmf()
const
{
    mmfs::ces mmf_obj;
    this->gen_mmf(mmf_obj);
    //
    return mmf_obj;
}

//
// DSE-related functions
// Implicit Parameterization

arma::mat 
trame::transfers::etu::Psi(const arma::mat& U, const arma::mat& V)
const
{
    return this->Psi(U,V,NULL,NULL);
}

arma::mat 
trame::transfers::etu::Psi(const arma::mat& U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = elem_prod(arma::exp(elem_div(U, tau(x_ind,y_ind))), aux_alpha_exp(x_ind,y_ind));
    arma::mat term_2 = elem_prod(arma::exp(elem_div(V, tau(x_ind,y_ind))), aux_gamma_exp(x_ind,y_ind));

    arma::mat ret =  elem_prod(tau(x_ind,y_ind), arma::log(0.5 * (term_1 + term_2)));
    //
    return ret;
}

arma::mat 
trame::transfers::etu::Psi(const double& U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = elem_prod(arma::exp(U/tau(x_ind,y_ind)), aux_alpha_exp(x_ind,y_ind));
    arma::mat term_2 = elem_prod(arma::exp(elem_div(V, tau(x_ind,y_ind))), aux_gamma_exp(x_ind,y_ind));
        
    arma::mat ret =  elem_prod(tau(x_ind,y_ind), arma::log(0.5 * (term_1 + term_2)));
    //
    return ret;
}

arma::mat 
trame::transfers::etu::Psi(const arma::mat& U, const double& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = elem_prod(arma::exp(elem_div(U, tau(x_ind,y_ind))), aux_alpha_exp(x_ind,y_ind));
    arma::mat term_2 = elem_prod(arma::exp(V/tau(x_ind,y_ind)), aux_gamma_exp(x_ind,y_ind));

    arma::mat ret =  elem_prod(tau(x_ind,y_ind), arma::log(0.5 * (term_1 + term_2)));
    //
    return ret;
}

double 
trame::transfers::etu::Psi(const double& U, const double& V, int x_ind, int y_ind)
const
{
    double term_1 = std::exp(U/tau(x_ind,y_ind)) * aux_alpha_exp(x_ind,y_ind);
    double term_2 = std::exp(V/tau(x_ind,y_ind)) * aux_gamma_exp(x_ind,y_ind);

    double ret =  tau(x_ind,y_ind) * std::log(0.5 * (term_1 + term_2));
    //
    return ret;
}

//
// Derivative of Psi wrt u

arma::mat
trame::transfers::etu::du_Psi(const arma::mat& U, const arma::mat& V)
const
{
    return this->du_Psi(U,V,NULL,NULL);
}

arma::mat
trame::transfers::etu::du_Psi(const arma::mat& U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret =  1 / (1 + arma::exp((V - U + alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind)));
    //
    return ret;
}

arma::mat 
trame::transfers::etu::du_Psi(const double& U, const arma::mat& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret =  1 / (1 + arma::exp((V - U + alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind)));
    //
    return ret;
}

arma::mat 
trame::transfers::etu::du_Psi(const arma::mat& U, const double& V, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret =  1 / (1 + arma::exp((V - U + alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind)));
    //
    return ret;
}

arma::mat 
trame::transfers::etu::dparams_Psi(const arma::mat& U, const arma::mat& V, const arma::mat& dparams)
{
    return this->dparams_Psi(U,V,&dparams);
}

arma::mat 
trame::transfers::etu::dparams_Psi(const arma::mat& U, const arma::mat& V, const arma::mat* dparams)
{
    arma::mat ret(nbX,nbY);
    //
    arma::mat du_psi_mat = du_Psi(U,V,NULL,NULL);
    arma::vec du_psi = arma::vectorise(du_psi_mat);

    if (!dparams) {
        arma::mat term_1 = (U - alpha) % du_psi_mat;
        arma::mat term_2 = (V - gamma) % (1 - du_psi_mat);

        arma::mat dsigma_psi_mat = (Psi(U,V,NULL,NULL) - term_1 - term_2)/tau;
        arma::vec dsigma_psi = arma::vectorise(dsigma_psi_mat);
        //
        ret = arma::join_rows(arma::diagmat(-du_psi),arma::join_rows(arma::diagmat(du_psi-1),arma::diagmat(dsigma_psi)));
    } else {
        arma::mat dparams_1 = dparams->rows(0,nbX*nbY-1);
        arma::mat dparams_2 = dparams->rows(nbX*nbY,2*nbX*nbY-1);
        arma::mat dparams_3 = dparams->rows(2*nbX*nbY,3*nbX*nbY-1);

        arma::vec dsigma_psi_dparams = arma::zeros(dparams_3.n_rows,1);
        double min_check = elem_min(dparams_3);

        if (min_check!=0) {
            arma::mat term_1 = (U - alpha) % du_psi_mat;
            arma::mat term_2 = (V - gamma) % (1 - du_psi_mat);

            arma::mat dsigma_psi_mat = (Psi(U,V,NULL,NULL) - term_1 - term_2)/tau;

            dsigma_psi_dparams = dparams_3 % arma::vectorise(dsigma_psi_mat);
        }
        //
        ret = arma::vectorise(-du_psi % dparams_1 - (1-du_psi) % dparams_2 + dsigma_psi_dparams);
    }
    //
    return ret;
}

//
// Explicit Parameterization

arma::mat 
trame::transfers::etu::Ucal(const arma::mat& vs)
const
{
    return this->Ucal(vs,NULL,NULL);
}

arma::mat 
trame::transfers::etu::Ucal(const arma::mat& vs, const arma::uvec* xs, const arma::uvec* ys)
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

double 
trame::transfers::etu::Ucal(const double& vs, int xs, int ys)
const
{
    double term_1 = vs - gamma(xs,ys);
    double term_log = 2 - std::exp(term_1/tau(xs,ys));

    double ret = alpha(xs,ys) + tau(xs,ys) * std::log(term_log);
    //
    return ret;
}

arma::mat 
trame::transfers::etu::Vcal(const arma::mat& us)
const
{
    return this->Vcal(us,NULL,NULL);
}

arma::mat 
trame::transfers::etu::Vcal(const arma::mat& us, const arma::uvec* xs, const arma::uvec* ys)
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

double 
trame::transfers::etu::Vcal(const double& us, int xs, int ys)
const
{
    double term_1 = us - alpha(xs,ys);
    double term_log = 2 - std::exp(term_1/tau(xs,ys));

    double ret = gamma(xs,ys) + tau(xs,ys) * std::log(term_log);
    //
    return ret;
}

arma::mat 
trame::transfers::etu::UW(const arma::mat& Ws)
const
{
    return this->UW(Ws,NULL,NULL);
}

arma::mat 
trame::transfers::etu::UW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - Psi(0.0,-Ws,&x_ind,&y_ind);

    return ret;
}

double 
trame::transfers::etu::UW(const double& Ws, int x_ind, int y_ind)
const
{
    double ret = - Psi((double) 0.0,-Ws,x_ind,y_ind);
    //
    return ret;
}

arma::mat 
trame::transfers::etu::VW(const arma::mat& Ws)
const
{
    return this->VW(Ws,NULL,NULL);
}

arma::mat 
trame::transfers::etu::VW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - Psi(Ws,0.0,&x_ind,&y_ind);

    return ret;
}

double 
trame::transfers::etu::VW(const double& Ws, int x_ind, int y_ind)
const
{
    double ret = - Psi(Ws,(double) 0.0,x_ind,y_ind);
    //
    return ret;
}

arma::mat 
trame::transfers::etu::dw_UW(const arma::mat& Ws)
const
{
    return this->dw_UW(Ws,NULL,NULL);
}

arma::mat 
trame::transfers::etu::dw_UW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = 1.0 - du_Psi(0.0,-Ws,&x_ind,&y_ind);

    return ret;
}

arma::mat 
trame::transfers::etu::dw_VW(const arma::mat& Ws)
const
{
    return this->dw_VW(Ws,NULL,NULL);
}

arma::mat 
trame::transfers::etu::dw_VW(const arma::mat& Ws, const arma::uvec* xs, const arma::uvec* ys)
const
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat ret = - du_Psi(Ws,0.0,&x_ind,&y_ind);

    return ret;
}

arma::mat 
trame::transfers::etu::WU(const arma::mat& Us)
{
    return this->WU(Us,NULL,NULL);
}

arma::mat 
trame::transfers::etu::WU(const arma::mat& Us, const arma::uvec* xs, const arma::uvec* ys)
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

arma::mat 
trame::transfers::etu::WV(const arma::mat& Vs)
{
    return this->WV(Vs,NULL,NULL);
}

arma::mat 
trame::transfers::etu::WV(const arma::mat& Vs, const arma::uvec* xs, const arma::uvec* ys)
{
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, nbX-1);
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat term_1 = 2 * arma::exp( (alpha(x_ind,y_ind) - Vs)/tau(x_ind,y_ind) );
    arma::mat term_2 = arma::exp( (alpha(x_ind,y_ind) - gamma(x_ind,y_ind))/tau(x_ind,y_ind) );
    arma::mat term_log = term_1 - term_2;

    arma::mat ret = - tau(x_ind,y_ind) % arma::log(term_log);
    //
    return ret;
}

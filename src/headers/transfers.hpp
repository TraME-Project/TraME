/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 Alfred Galichon and the TraME Team
  ##
  ##   This file is part of the R package TraME.
  ##
  ##   The R package TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

// transfers

class transfers
{
    public:
        // build objects
        bool ETU = false;
        bool LTU = false;
        bool NTU = false;
        bool TU  = false;

        int nbX;
        int nbY;
        int nbParams;

        arma::mat alpha;
        arma::mat gamma;
        arma::mat lambda
        arma::mat tau;

        arma::mat phi;

        // member functions
        void build_ETU(arma::mat alpha_ETU, arma::mat gamma_ETU, arma::mat tau_ETU);
        void build_LTU(arma::mat lambda_LTU, arma::mat phi_LTU);
        void build_NTU(arma::mat alpha_NTU, arma::mat gamma_NTU);
        void build_TU(arma::mat phi_TU);

    private:
        arma::mat aux_exp_alphaovertau;
        arma::mat aux_exp_gammaovertau;

        arma::mat aux_zeta;
}

void transfers::build_ETU(arma::mat alpha_ETU, arma::mat gamma_ETU, arma::mat tau_ETU)
{
    alpha = alpha_ETU;
    gamma = gamma_ETU;
    tau   = tau_ETU;

    nbX = alpha_ETU.n_rows;
    nbY = alpha_ETU.n_cols;
    nbParams = 3*nbX*nbY;

    aux_exp_alphaovertau = arma::exp(- alpha / tau);
    aux_exp_gammaovertau = arma::exp(- gamma / tau);

    ETU = true;
}

void transfers::build_LTU(arma::mat lambda_LTU, arma::mat phi_LTU)
{
    lambda = lambda_LTU;
    phi    = phi_LTU;

    nbX = lambda_LTU.n_rows;
    nbY = lambda_LTU.n_cols;
    nbParams = 2*nbX*nbY;

    aux_zeta = 1 - lambda;

    LTU = true;
}

void transfers::build_NTU(arma::mat alpha_NTU, arma::mat gamma_NTU)
{
    alpha = alpha_NTU;
    gamma = gamma_NTU;

    nbX = alpha_NTU.n_rows;
    nbY = alpha_NTU.n_cols;
    nbParams = 2*nbX*nbY;

    NTU = true;
}

void transfers::build_TU(arma::mat phi_TU)
{
    phi = phi_TU;

    nbX = phi_TU.n_rows;
    nbY = phi_TU.n_cols;
    nbParams = nbX*nbY;

    TU = true;
}

void transfers::trans()
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

        aux_exp_alphaovertau_temp = aux_exp_alphaovertau;

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

arma::mat transfers::Psi(arma::mat U, arma::mat V)
{
    arma::mat ret(nbX,nbY);
    //
    if (ETU) {
        ret =  tau % arma::log(0.5 * (arma::exp(U/tau) % aux_exp_alphaovertau + arma::exp(V/tau) % aux_exp_gammaovertau));
        goto finished;
    }

    if (LTU) {
        ret = lambda % U + aux_zeta % V - phi;
        goto finished;
    }

    if (NTU) {
        ret = arma::max(U - alpha, V - gamma);
        goto finished;
    }

    if (TU) {
        ret = (U + V - phi) / 2;
        goto finished;
    }
    //
finished:
    return ret
}

arma::mat transfers::du_Psi(arma::mat U, arma::mat V)
{
    arma::mat ret(nbX,nbY);
    //
    if (ETU) {
        ret =  1 / (1 + arma::exp((V - U + alpha - gamma)/tau));
        goto finished;
    }

    if (LTU) {
        ret = lambda;
        goto finished;
    }

    if (NTU) {
        if (U - alpha >= V - gamma) {
            ret.ones();
        } else {
            ret.zeros();
        }
        goto finished;
    }

    if (TU) {
        ret.fill(0.5);
        goto finished;
    }
    //
finished:
    return ret
}



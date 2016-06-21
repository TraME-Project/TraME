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

// MMF class

class MMF
{
    public:
        // build objects
        bool ETU = false;
        bool LTU = false;
        bool NTU = false;
        bool TU  = false;

        int n;
        int m;

        double kappa;
        double lambda;

        arma::mat A;
        arma::mat B;
        arma::mat C;
        arma::mat D;
        arma::mat K;

        // member functions

    private:
        arma::mat aux_log_C;
        arma::mat aux_log_D;
}

void MMF::build(int n_ETU, int m_ETU, arma::mat C_ETU, arma::mat D_ETU, double kappa_ETU, bool need_norm_ETU)
{
    n = n_ETU;
    m = m_ETU;

    C = C_ETU;
    D = D_ETU;

    aux_log_C = arma::log(C);
    aux_log_D = arma::log(D);

    kappa = kappa_ETU;
    
    need_norm = need_norm_ETU;

    ETU = true;
}

void MMF::build(int n_LTU, int m_LTU, arma::mat K_LTU, double lambda_LTU, bool need_norm_LTU)
{
    n = n_LTU;
    m = m_LTU;
    K = K_LTU;

    lambda = lambda_LTU;
    
    need_norm = need_norm_LTU;

    LTU = true;
}

void MMF::build(int n_NTU, int m_NTU, arma::mat A_NTU, arma::mat B_NTU, bool need_norm_NTU)
{
    n = n_NTU;
    m = m_NTU;
    
    A = A_NTU;
    B = B_NTU;
    
    need_norm = need_norm_NTU;

    NTU = true;
}

void MMF::build(int n_TU, int m_TU, arma::mat K_TU, bool need_norm_TU)
{
    n = n_TU;
    m = m_TU;
    K = K_TU;

    need_norm = need_norm_TU;

    TU = true;
}

arma::mat MMF::M(arma::mat a_xs, arma::mat b_ys, arma::uvec* xs, arma::uvec* ys)
{
    if (ETU) {
        arma::mat term_1 = arma::exp(aux_log_C(xs,ys) + kappa(xs,ys) % arma::log(a_xs));
        arma::mat term_2 = arma::exp(aux_log_D(xs,ys) + arma::trans(arma::trans(kappa(xs,ys)) % arma::log(b_ys)));

        return ret = arma::exp((1/kappa(xs,ys)) % arma::log((term_1+term_2)/2)); 
    }
    //
    if (LTU) {
        arma::mat term_1 = arma::exp(lambda(xs,ys) % arma::log(a_xs));
        arma::mat term_2 = arma::exp(arma::trans(aux_zeta(xs,ys)) % arma::log(b_ys));
    }
    //
    if (NTU) {

    }
    //
    if (TU) {

    }
}

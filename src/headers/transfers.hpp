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
        arma::mat tau;

        arma::mat phi;

        // member functions
        void build_ETU(arma::mat alpha_ETU, arma::mat gamma_ETU, arma::mat tau_ETU);
        void build_LTU(arma::vec n_LTU, arma::vec m_LTU, arma::mat K_LTU, arma::mat lambda_LTU, bool need_norm_LTU);
        void build_NTU(arma::vec n_NTU, arma::vec m_NTU, arma::mat A_NTU, arma::mat B_NTU, bool need_norm_NTU);
        void build_TU(arma::vec n_TU, arma::vec m_TU, arma::mat K_TU, bool need_norm_TU);


    private:
        aux_exp_alphaovertau;
        aux_exp_gammaovertau;
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
}

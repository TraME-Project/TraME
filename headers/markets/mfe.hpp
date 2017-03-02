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
 * Matching Function Equilibrium (MFE) market
 *
 * Keith O'Hara
 * 08/16/2016
 */

template<class Tm>
class mfe
{
    public:
        // build objects
        bool ETU = false;
        bool LTU = false;
        bool NTU = false;
        bool TU  = false;

        bool need_norm;
        bool outsideOption;

        int nbX;
        int nbY;

        double sigma;

        arma::vec n;
        arma::vec m;

        logit arums_G;
        logit arums_H;

        transfers trans_obj;
        Tm mmf_obj;

        // member functions
        void build_ETU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, double* sigma_inp, bool need_norm_inp);
        void build_LTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, double* sigma_inp, bool need_norm_inp);
        void build_NTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, double* sigma_inp, bool need_norm_inp);
        void build_TU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, double* sigma_inp, bool need_norm_inp);

        void trans();

        bool solve(arma::mat& mu_sol);
        bool solve(arma::mat& mu_sol, const char* solver);
};

#include "mfe.tpp"

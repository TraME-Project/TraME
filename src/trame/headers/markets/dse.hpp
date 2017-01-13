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
 * Demand-Supply Equilibrium (DSE) market
 *
 * Keith O'Hara
 * 08/17/2016
 *
 * This version:
 * 11/02/2016
 */

template<class Ta>
class dse
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

        //mmf mmf_obj;
        transfers trans_obj;

        Ta arums_G;
        Ta arums_H;

        // member functions
        void build_TU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, bool need_norm_inp);
        void build_TU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, const Ta& arums_G_inp, const Ta& arums_H_inp, bool need_norm_inp);
        template<typename T> void build_TU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, T arums_G_inp, T arums_H_inp, int nbDraws, int seed, bool need_norm_inp);

        void build_NTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, bool need_norm_inp);
        void build_NTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const Ta& arums_G_inp, const Ta& arums_H_inp, bool need_norm_inp);
        template<typename T> void build_NTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, T arums_G_inp, T arums_H_inp, int nbDraws, int seed, bool need_norm_inp);

        void build_LTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda, const arma::mat& phi, bool need_norm_inp);
        void build_LTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, const Ta& arums_G_inp, const Ta& arums_H_inp, bool need_norm_inp);
        template<typename T> void build_LTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, T arums_G_inp, T arums_H_inp, int nbDraws, int seed, bool need_norm_inp);

        void trans();

        bool solve(arma::mat& mu_sol, const char* solver = "none");

    private:
        bool arum_none;
        bool arum_empirical;
        bool arum_general; // need to finish this later
};

#include "dse.tpp"

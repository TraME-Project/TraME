/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 the TraME Team:
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
 */

template <class Ta>
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
        void build_TU(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp);
        void build_TU(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, Ta arums_G_inp, Ta arums_H_inp, bool need_norm_inp);
        template <typename T> void build_TU(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, T arums_G_inp, T arums_H_inp, int nbDraws, int seed, bool need_norm_inp);

        void build_NTU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp);
        void build_NTU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, Ta arums_G_inp, Ta arums_H_inp, bool need_norm_inp);
        template <typename T> void build_NTU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, T arums_G_inp, T arums_H_inp, int nbDraws, int seed, bool need_norm_inp);

        void build_LTU(arma::vec n_inp, arma::vec m_inp, arma::mat lambda, arma::mat phi, bool need_norm_inp);
        void build_LTU(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, Ta arums_G_inp, Ta arums_H_inp, bool need_norm_inp);
        template <typename T> void build_LTU(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, T arums_G_inp, T arums_H_inp, int nbDraws, int seed, bool need_norm_inp);

        void trans();

    private:
        bool arum_none;
        bool arum_empirical;
        bool arum_general; // need to finish this later
};

#include "dse.tpp"

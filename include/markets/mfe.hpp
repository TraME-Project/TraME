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
 * Matching Function Equilibrium (MFE) market
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 07/26/2017
 */

#ifndef _trame_mfe_market_HPP
#define _trame_mfe_market_HPP

template<class Tt>
class mfe
{
    public:
        // build objects
        bool need_norm = false;
        bool outside_option = true;

        int nbX;
        int nbY;

        double sigma = 1.0;

        arma::vec n;
        arma::vec m;

        Tt mmfs_obj;

        // member functions
        ~mfe(){};
         mfe(){};
        explicit mfe(const arma::vec& n_inp, const arma::vec& m_inp);
        explicit mfe(const double sigma_inp, const bool need_norm_inp);
        explicit mfe(const arma::vec& n_inp, const arma::vec& m_inp, const double sigma_inp, const bool need_norm_inp);

        void build(const arma::vec& n_inp, const arma::vec& m_inp);
        void build(const double sigma_inp, const bool need_norm_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, const double sigma_inp, const bool need_norm_inp);

        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& mmfs_params_inp_1);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& mmfs_params_inp_1, const arma::mat& mmfs_params_inp_2);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& mmfs_params_inp_1, const arma::mat& mmfs_params_inp_2, const arma::mat& mmfs_params_inp_3);

        void trans(); // market transpose

        arma::vec marg_x_inv(const arma::mat& B_ys) const;
        arma::vec marg_x_inv(const arma::mat& B_ys, const arma::uvec* xs) const;
        arma::vec marg_y_inv(const arma::mat& A_xs) const;
        arma::vec marg_y_inv(const arma::mat& A_xs, const arma::uvec* ys) const;

        bool solve(arma::mat& mu_sol);
        bool solve(arma::mat& mu_sol, const char* solver);
        bool solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver);

    private:
        static double marg_x_inv_fn(const double z, void* opt_data);
        static double marg_y_inv_fn(const double z, void* opt_data);
};

template<class Tt>
struct trame_mfe_zeroin_data {
    bool coeff;

    int x_ind;
    int y_ind;

    arma::mat A_xs;
    arma::mat B_ys;

    mfe<Tt> mfe_obj;
};

#include "mfe.tpp"

#endif

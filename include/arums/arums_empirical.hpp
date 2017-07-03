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
 * empirical additive random utility model (ARUM) class
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 07/03/2017
 */

#ifndef _trame_arums_empirical_HPP
#define _trame_arums_empirical_HPP

class empirical
{
    public:
        // build objects
        int nbX;
        int nbY;

        int dim_params;
        int aux_nbDraws;
        int nb_options;

        bool x_homogeneous;
        bool outside_option;

        arma::cube atoms;

        // input objects
        arma::mat U;
        arma::mat mu;

        // equilibrium objects
        arma::mat U_sol;
        arma::mat mu_sol;

        // member functions
        ~empirical(){};
         empirical(){};
        explicit empirical(const int nbX_inp, const int nbY_inp);
        explicit empirical(const int nbX_inp, const int nbY_inp, const arma::cube& atoms_inp, const bool x_homogeneous_inp, const bool outside_option_inp);

        void build(const int nbX_inp, const int nbY_inp);
        void build(const int nbX_inp, const int nbY_inp, const arma::cube& atoms_inp, const bool x_homogeneous_inp, const bool outside_option_inp);

        double G(const arma::vec& n);
        double G(const arma::vec& n, const arma::mat& U_inp, arma::mat& mu_out) const;
        double Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, const int x) const;

        double Gstar(const arma::vec& n);
        double Gstar(const arma::vec& n, const arma::mat& mu_inp, arma::mat& U_out) const;
        double Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x) const;

        double Gbar(const arma::mat& Ubar, const arma::mat& mubar, const arma::vec& n, arma::mat& U_out, arma::mat& mu_out) const;
        double Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, const int x) const;

    private:
        /*
         * these private member objects are used by the LP solvers
         */
        void presolve_LP_Gstar();
        void presolve_LP_Gbar();

        bool TRAME_PRESOLVED_GSTAR = false; // initialization requires C++11
        bool TRAME_PRESOLVED_GBAR  = false;

        int k_Gstar;
        int n_Gstar;
        int numnz_Gstar;
        int* vind_Gstar;
        int* vbeg_Gstar;
        double* vval_Gstar;

        int k_Gbar;
        int n_Gbar;
        int numnz_Gbar;
        int* vind_Gbar;
        int* vbeg_Gbar;
        double* vval_Gbar;
};

#endif

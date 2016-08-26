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
 * empirical class
 *
 * Keith O'Hara
 * 08/08/2016
 */

class empirical
{
    public:
        // build objects
        int nbX;
        int nbY;

        int nbParams;
        int aux_nbDraws;
        int nbOptions;

        bool xHomogenous;
        bool outsideOption;

        arma::cube atoms;

        // input objects
        arma::mat U;
        arma::mat mu;

        // equilibrium objects
        arma::mat U_sol;
        arma::mat mu_sol;

        // member functions
        void build(int nbX_inp, int nbY_inp, arma::cube atoms_inp, bool xHomogenous_inp, bool outsideOption_inp);

        double G(arma::vec n);
        double Gx(arma::mat Ux, arma::mat& mu_x_inp, int x);

        double Gstar(arma::vec n);
        double Gstar(arma::mat& U_inp, arma::vec n);
        double Gstarx(arma::mat& Ux_inp, arma::mat mu_x, int x);

        double Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_inp, arma::mat& mu_inp);
        double Gbarx(arma::vec Ubarx, arma::vec mubarx, arma::mat& Ux_inp, arma::mat& mu_x_inp, int x);

    private:
        /*
         * these private member objects are mostly for use with Gurobi
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

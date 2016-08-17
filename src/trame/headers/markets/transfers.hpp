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
 * transfers class
 *
 * Keith O'Hara
 * 08/16/2016
 */

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
        arma::mat lambda;
        arma::mat tau;

        arma::mat phi;

        // member functions
        void build_ETU(arma::mat alpha_ETU, arma::mat gamma_ETU, arma::mat tau_ETU);
        void build_LTU(arma::mat lambda_LTU, arma::mat phi_LTU);
        void build_NTU(arma::mat alpha_NTU, arma::mat gamma_NTU);
        void build_TU(arma::mat phi_TU);

        void trans();

        arma::mat Psi(arma::mat U, arma::mat V);
        arma::mat Psi(double U, arma::mat V, arma::uvec xs, arma::uvec ys);
        arma::mat Psi(arma::mat U, double V, arma::uvec xs, arma::uvec ys);
        arma::mat Psi(arma::mat U, arma::mat V, arma::uvec xs, arma::uvec ys);

        arma::mat du_Psi(arma::mat U, arma::mat V);
        arma::mat du_Psi(double U, arma::mat V, arma::uvec xs, arma::uvec ys);
        arma::mat du_Psi(arma::mat U, double V, arma::uvec xs, arma::uvec ys);
        arma::mat du_Psi(arma::mat U, arma::mat V, arma::uvec xs, arma::uvec ys);

        arma::mat dtheta_Psi(arma::mat U, arma::mat V, arma::mat* dtheta);

        arma::mat Ucal(arma::mat vs, arma::uvec xs, arma::uvec ys);
        arma::mat Vcal(arma::mat us, arma::uvec xs, arma::uvec ys);

        arma::mat UW(arma::mat Ws, arma::uvec xs, arma::uvec ys);
        arma::mat VW(arma::mat Ws, arma::uvec xs, arma::uvec ys);
        arma::mat du_UW(arma::mat Ws, arma::uvec xs, arma::uvec ys);
        arma::mat du_VW(arma::mat Ws, arma::uvec xs, arma::uvec ys);

        arma::mat WU(arma::mat Us, arma::uvec xs, arma::uvec ys);
        arma::mat WV(arma::mat Vs, arma::uvec xs, arma::uvec ys);

    private:
        arma::mat aux_exp_alphaovertau;
        arma::mat aux_exp_gammaovertau;

        arma::mat aux_zeta;
};

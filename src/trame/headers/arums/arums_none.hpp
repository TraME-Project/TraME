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
 * none class
 *
 * Keith O'Hara
 * 08/08/2016
 */

class none
{
    public:
        // build objects
        int nbX;
        int nbY;
        int nbParams;
        
        // input objects
        arma::mat mu;
        arma::mat U;
        
        // equilibrium objects
        arma::mat mu_sol;
        arma::mat U_sol;
        
        // member functions
        void build(int nbX_inp, int nbY_inp);
        double G(arma::vec n);
        double Gx(arma::mat Ux, arma::mat& mu_x_inp);

        double Gstar(arma::mat& U_inp, arma::vec n);
        
        double Gbar(arma::mat Ubarx, arma::mat mubarx, arma::vec n, arma::mat& Ux_inp, arma::mat& mu_x_inp);
        double Gbarx(arma::mat Ubarx, arma::mat mubarx, arma::mat& Ux_inp, arma::mat& mu_x_inp);
        
        arma::vec dtheta_NablaGstar();
        
        void simul(empirical &ret, int nbDraws, int seed);
};

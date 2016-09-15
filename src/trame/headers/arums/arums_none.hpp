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
        ~none(){};
         none(){};
        explicit none(int nbX_inp, int nbY_inp);
        
        void build(int nbX_inp, int nbY_inp);
        
        double G(arma::vec n);
        double G(arma::vec n, const arma::mat& U_inp, arma::mat& mu_out);
        double Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out);
        double Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, int x);

        double Gstar(arma::vec n);
        double Gstar(arma::vec n, const arma::mat& mu_inp, arma::mat& U_out);
        double Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x);
        
        double Gbar(arma::mat Ubarx, arma::mat mubarx, arma::vec n, arma::mat& U_x_out, arma::mat& mu_x_out);
        double Gbarx(arma::mat Ubarx, arma::mat mubarx, arma::mat& U_x_out, arma::mat& mu_x_out);
        
        arma::vec dtheta_NablaGstar();
        
        empirical simul();
        empirical simul(int* nbDraws, int* seed);
        void simul(empirical& obj_out);
        void simul(empirical& obj_out, int* nbDraws, int* seed);
};

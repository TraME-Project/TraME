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
 * logit class
 *
 * Keith O'Hara
 * 08/08/2016
 */

class logit
{
    public:
        // build objects
        int nbX;
        int nbY;
        int nbParams;
        double sigma;
        bool outsideOption;
        
        // input objects
        arma::mat mu;
        arma::mat U;
        
        // equilibrium objects
        arma::mat mu_sol;
        arma::mat U_sol;
        
        // member functions
        void build(int nbX_inp, int nbY_inp, double sigma_inp, bool outsideOption_inp);
        double G(arma::vec n);
        double G(arma::mat &mu_ret, arma::vec n);
        
        double Gstar(arma::vec n);
        double Gstar(arma::mat &U_ret, arma::vec n);
        double Gstarx(arma::mat &U_x_ret, arma::mat mu_x);
        
        double Gbar(arma::mat U_bar, arma::mat mu_bar, arma::vec n, arma::mat& U_inp, arma::mat& mu_inp);
        double Gbarx(arma::mat U_bar_x, arma::mat mu_bar_x, arma::mat& U_x_inp, arma::mat& mu_x_inp);
        
        void D2G(arma::mat &H, arma::vec n, bool xFirst);
        void D2Gstar(arma::mat &H, arma::vec n, bool xFirst);
        void dtheta_NablaGstar(arma::mat &ret, arma::vec n, arma::mat dtheta, bool xFirst);
        
        void simul(empirical &ret, int nbDraws, int seed);
};

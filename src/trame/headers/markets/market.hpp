/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 Alfred Galichon and the TraME Team
  ##
  ##   This file is part of the R package TraME.
  ##
  ##   The R package TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

// logit class
class market
{
    public:
        // market objects
        arma::vec types; // these are binary analogs of "itu-rum" and "mfe"
        int n;
        int m;
        void neededNorm;
        
        double sigma;
        int outsideOption; // 1 = yes, 0 = no
        
        // equilibrium objects
        arma::mat mu;
        arma::mat mux;
        arma::mat U;
        arma::mat Ux;
        
        // member functions
        void build(int nbX_b, int nbY_b, int nbParams_b, double sigma_b, int outsideOption_b);
        void G(double &val, arma::vec n);
        void G(double &val, arma::mat &mu_ret, arma::vec n);
        void Gstar(double &val, arma::vec n);
        void Gstar(double &val, arma::mat &U_ret, arma::vec n);
        void Gstarx(double &valx, arma::mat mux);
        void Gstarx(double &valx, arma::mat &Ux_ret, arma::mat mux);
        void D2G(arma::mat &H, arma::vec n, int xFirst);
        void D2Gstar(arma::mat &H, arma::vec n, int xFirst);
        void dtheta_NablaGstar(arma::mat &ret, arma::vec n, arma::mat dtheta, int xFirst);
        void Gbarx(double &valx, arma::mat Ubarx, arma::mat mubarx);
        void simul(empirical &ret, int nbDraws, int seed);
        
        
    //private:
};
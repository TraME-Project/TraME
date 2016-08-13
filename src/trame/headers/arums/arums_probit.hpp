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
 * probit class
 *
 * Keith O'Hara
 * 08/08/2016
 */

class probit
{
    public:
        // build objects
        int nbX;
        int nbY;
        int nbParams;
        int aux_nbOptions;
        
        bool outsideOption;
        
        double rho;
        
        arma::cube Covar;
        
        // member functions
        void build(int nbX_inp, int nbY_inp, bool outsideOption_inp);
        void simul(empirical &ret, int nbDraws, int seed);
        void unifCorrelCovMatrices();
        arma::cube unifCorrelCovMatrices(double rho);
};

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
 * probit additive random utility model (ARUM) class
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 06/10/2017
 */

#ifndef _trame_arums_probit_HPP
#define _trame_arums_probit_HPP

class probit
{
    public:
        // build objects
        int nbX;
        int nbY;
        int dim_params;
        int aux_nb_options;
        
        bool outside_option;
        
        double rho;
        
        arma::cube Covar;
        
        // member functions
        ~probit(){};
         probit(){};
        explicit probit(int nbX_inp, int nbY_inp);
        explicit probit(int nbX_inp, int nbY_inp, bool outside_option_inp);
        explicit probit(int nbX_inp, int nbY_inp, double rho_inp, bool outside_option_inp);

        void build(int nbX_inp, int nbY_inp);
        void build(int nbX_inp, int nbY_inp, bool outside_option_inp);
        void build(int nbX_inp, int nbY_inp, double rho_inp, bool outside_option_inp);

        void unifCorrelCovMatrices();
        void unifCorrelCovMatrices(double rho_inp);

        empirical simul() const;
        empirical simul(int nbDraws, int seed) const;
        void simul(empirical& obj_out) const;
        void simul(empirical& obj_out, int nbDraws, int seed) const;
    
    private:
         void build_int(int nbX_inp, int nbY_inp, double* rho_inp, bool outside_option_inp);

         void simul_int(empirical& obj_out, int* nbDraws, int* seed) const;
};

#endif

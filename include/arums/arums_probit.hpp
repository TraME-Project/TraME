/*################################################################################
  ##
  ##   Copyright (C) 2015-2018 the TraME Team:
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
 * 02/04/2018
 */

#ifndef _trame_arums_probit_HPP
#define _trame_arums_probit_HPP

class probit
{
    public:
        // build objects
        uint_t nbX;
        uint_t nbY;
        uint_t dim_params;
        uint_t aux_nb_options;
        
        bool outside_option;
        
        double rho;
        
        arma::cube Covar;
        
        // member functions
        ~probit(){};
         probit(){};
        explicit probit(const uint_t nbX_inp, const uint_t nbY_inp);
        explicit probit(const uint_t nbX_inp, const uint_t nbY_inp, const bool outside_option_inp);
        explicit probit(const uint_t nbX_inp, const uint_t nbY_inp, const double rho_inp, const bool outside_option_inp);

        void build(const uint_t nbX_inp, const uint_t nbY_inp);
        void build(const uint_t nbX_inp, const uint_t nbY_inp, const bool outside_option_inp);
        void build(const uint_t nbX_inp, const uint_t nbY_inp, const double rho_inp, const bool outside_option_inp);

        void unifCorrelCovMatrices();
        void unifCorrelCovMatrices(const double rho_inp);

        empirical simul() const;
        empirical simul(const uint_t n_draws, const uint_t seed) const;
        void simul(empirical& obj_out) const;
        void simul(empirical& obj_out, const uint_t n_draws, const uint_t seed) const;
    
    protected:
         void build_int(const uint_t nbX_inp, const uint_t nbY_inp, const double* rho_inp, const bool outside_option_inp);

         void simul_int(empirical& obj_out, const uint_t* n_draws_inp, const uint_t* seed) const;
};

#endif

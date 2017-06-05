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
 * Constant Elasticity of Substitution (CES) Matching Market Functions (MMFs) class
 * (Corresponds to the ETU transfers class)
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 04/27/2017
 */

// some functions are const restricted because of calls to const market pointers in equilibrium solvers

#ifndef _trame_mmfs_ces_HPP
#define _trame_mmfs_ces_HPP

class ces
{
    public:
        // build objects
        bool need_norm;

        int nbX;
        int nbY;
        int nbParams;

        arma::mat alpha;
        arma::mat gamma;
        arma::mat tau;
        
        arma::mat kappa; // kappa = -1/tau_inp

        arma::mat aux_alpha; // - alpha / tau = log(C)
        arma::mat aux_gamma; // - gamma / tau = log(D)

        arma::mat aux_alpha_exp; // exp(- alpha / tau), also labelled C
        arma::mat aux_gamma_exp; // exp(- gamma / tau), also labelled D

        // member functions
        ~ces(){};
         ces(){};

        void build(const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, bool need_norm_inp);

        void trans();

        //
        // MFE-related functions
        arma::mat M(const arma::mat& a_xs, const arma::mat& b_ys) const;
        arma::mat M(const arma::mat& a_xs, const arma::mat& b_ys, const arma::uvec* xs, const arma::uvec* ys) const;
        arma::mat M(const double& a_xs, const arma::mat& b_ys, const arma::uvec* xs, const arma::uvec* ys) const;
        arma::mat M(const arma::mat& a_xs, const double& b_ys, const arma::uvec* xs, const arma::uvec* ys) const;

        arma::mat dmu_x0(const arma::mat& a_xs, const arma::mat& b_ys) const;
        arma::mat dmu_0y(const arma::mat& a_xs, const arma::mat& b_ys) const;

        arma::mat dparams_M(const arma::mat& a_xs, const arma::mat& b_ys) const;
        arma::mat dparams_M(const arma::mat& a_xs, const arma::mat& b_ys, const arma::mat* delta_params_M) const;

        arma::mat Mx0(const arma::mat& a_x) const;
        arma::mat M0y(const arma::mat& b_y) const;
};

#endif

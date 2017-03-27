/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
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
 * Constant Elasticity of Substitution (CES) Marriage Matching Functions (MMFs) class
 * Corresponds to the ETU transfers class
 *
 * Keith O'Hara
 * 08/15/2016
 *
 * This version:
 * 03/27/2017
 */

// some functions are const restricted because of calls to const market pointers in equilibrium solvers

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
        arma::mat kappa; // kappa_ETU = -1/tau_inp

        arma::mat aux_alpha; // - alpha_ETU / tau_ETU = log(C_ETU)
        arma::mat aux_gamma; // - gamma_ETU / tau_ETU = log(D_ETU)

        arma::mat aux_alpha_exp; // exp(- alpha_ETU / tau_ETU), also labelled C_ETU
        arma::mat aux_gamma_exp; // exp(- gamma_ETU / tau_ETU), also labelled D_ETU

        // member functions
        void build(const arma::mat& alpha_ETU, const arma::mat& gamma_ETU, const arma::mat& tau_ETU, bool need_norm_ETU);
        
        void trans();

        //
        // MFE-related functions
        arma::mat M(const arma::mat& a_xs, const arma::mat& b_ys) const;
        arma::mat M(const arma::mat& a_xs, const arma::mat& b_ys, arma::uvec* xs, arma::uvec* ys) const;
        arma::mat M(const double& a_xs, const arma::mat& b_ys, arma::uvec* xs, arma::uvec* ys) const;
        arma::mat M(const arma::mat& a_xs, const double& b_ys, arma::uvec* xs, arma::uvec* ys) const;

        arma::mat Mx0(const arma::mat& a_x) const;
        arma::mat M0y(const arma::mat& b_y) const;
};

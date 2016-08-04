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

// guard against double inclusion
#pragma once

//
// for use in zeroin()
struct trame_opt_data {
    arma::mat* exp_U_bar_X;
    arma::mat* mu_bar_X;
};

class trame_zeroin_data
{
    public:
        //
        // for logit class
        arma::mat exp_U_bar_X;
        arma::mat mu_bar_X;
        //
        // for mmf class
        arma::mat A_xs;
        arma::mat B_ys;

        int x_ind;
        int y_ind;

        bool coeff;
};

//
// nlopt structures
typedef struct {
    int x;
    int nbY;
    arma::vec Ubar_x;
    arma::mat zeta;
    arma::mat aux_DinvPsigma;
    arma::mat aux_Psigma;
    arma::mat aux_Influence_lhs;
    arma::mat aux_Influence_rhs;
    arma::vec (*pot_eps_vec)(arma::vec x, double* dist_pars);
    arma::vec (*quantile_eps_vec)(arma::vec x, double* dist_pars);
    double* dist_pars;
} trame_nlopt_opt_data;

typedef struct {
    int nbY;
} trame_nlopt_constr_data;

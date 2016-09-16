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
 * mmf class
 *
 * Keith O'Hara
 * 08/15/2016
 */

class mmf
{
    public:
        // build objects
        bool ETU = false;
        bool LTU = false;
        bool NTU = false;
        bool TU  = false;

        bool need_norm = false;

        arma::vec n;
        arma::vec m;

        arma::mat kappa;
        arma::mat lambda;

        arma::mat A;
        arma::mat B;
        arma::mat C;
        arma::mat D;
        arma::mat K;

        // member functions
        void build_ETU(arma::vec n_ETU, arma::vec m_ETU, arma::mat C_ETU, arma::mat D_ETU, arma::mat kappa_ETU, bool need_norm_ETU);
        void build_LTU(arma::vec n_LTU, arma::vec m_LTU, arma::mat lambda_LTU, arma::mat K_LTU, bool need_norm_LTU);
        void build_NTU(arma::vec n_NTU, arma::vec m_NTU, arma::mat A_NTU, arma::mat B_NTU, bool need_norm_NTU);
        void build_TU(arma::vec n_TU, arma::vec m_TU, arma::mat K_TU, bool need_norm_TU);

        arma::mat M(arma::mat a_xs, arma::mat b_ys);
        arma::mat M(arma::mat a_xs, arma::mat b_ys, arma::uvec* xs, arma::uvec* ys);
        arma::mat M(double a_xs, arma::mat b_ys, arma::uvec* xs, arma::uvec* ys);
        arma::mat M(arma::mat a_xs, double b_ys, arma::uvec* xs, arma::uvec* ys);

        arma::mat Mx0(arma::mat a_x);
        arma::mat M0y(arma::mat b_y);

        arma::vec marg_x_inv(arma::uvec* xs, arma::mat B_ys);
        arma::vec marg_y_inv(arma::uvec* ys, arma::mat A_xs);

        void trans();

    private:
        arma::mat aux_zeta;
        arma::mat aux_log_C;
        arma::mat aux_log_D;

        // member functions
        double marg_x_inv_fn(double z, const trame_zeroin_data& opt_data);
        double marg_y_inv_fn(double z, const trame_zeroin_data& opt_data);

        double zeroin_mmf(double ax, double bx, double (mmf::*f)(double x, const trame_zeroin_data& opt_data), const trame_zeroin_data& zeroin_data, double* tol_inp, int* max_iter_inp);
};

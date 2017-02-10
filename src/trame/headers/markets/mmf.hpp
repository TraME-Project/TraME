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
        void build_ETU(const arma::vec& n_ETU, const arma::vec& m_ETU, const arma::mat& C_ETU, const arma::mat& D_ETU, const arma::mat& kappa_ETU, bool need_norm_ETU);
        void build_LTU(const arma::vec& n_LTU, const arma::vec& m_LTU, const arma::mat& lambda_LTU, const arma::mat& K_LTU, bool need_norm_LTU);
        void build_NTU(const arma::vec& n_NTU, const arma::vec& m_NTU, const arma::mat& A_NTU, const arma::mat& B_NTU, bool need_norm_NTU);
        void build_TU(const arma::vec& n_TU, const arma::vec& m_TU, const arma::mat& K_TU, bool need_norm_TU);

        // M functions are const restricted because of calls to const market pointers in solvers
        arma::mat M(const arma::mat& a_xs, const arma::mat& b_ys) const;
        arma::mat M(const arma::mat& a_xs, const arma::mat& b_ys, arma::uvec* xs, arma::uvec* ys) const;
        arma::mat M(const double& a_xs, const arma::mat& b_ys, arma::uvec* xs, arma::uvec* ys) const;
        arma::mat M(const arma::mat& a_xs, const double& b_ys, arma::uvec* xs, arma::uvec* ys) const;

        arma::mat Mx0(const arma::mat& a_x) const;
        arma::mat M0y(const arma::mat& b_y) const;

        arma::vec marg_x_inv(const arma::mat& B_ys);
        arma::vec marg_x_inv(const arma::mat& B_ys, arma::uvec* xs);
        arma::vec marg_y_inv(const arma::mat& A_xs);
        arma::vec marg_y_inv(const arma::mat& A_xs, arma::uvec* ys);

        void trans();

    private:
        arma::mat aux_zeta;
        arma::mat aux_log_C;
        arma::mat aux_log_D;

        // member functions
        static double marg_x_inv_fn(double z, void* opt_data);
        static double marg_y_inv_fn(double z, void* opt_data);
};

struct trame_mmf_zeroin_data {
    bool coeff;

    int x_ind;
    int y_ind;

    arma::mat A_xs;
    arma::mat B_ys;

    mmf mmf_obj;
};

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
 * Cupids LP solver
 *
 * Keith O'Hara
 * 08/25/2016
 */

template<typename Ta>
bool cupids_lp(dse<Ta> market, bool xFirst, arma::mat& mu, arma::vec& mux0, arma::vec& mu0y, arma::mat& U_out, arma::mat& V_out)
{
    if (market.need_norm) {
        printf("CupidsLP does not yet allow for the case without unmatched agents.\n");
        return false;
    }
    //
    bool success = false;

    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;

    Ta arums_G = market.arums_G;
    Ta arums_H = market.arums_H;

    transfers trans_obj = market.trans_obj;

    arma::mat phi = trans_obj.phi;
    //
    arma::mat epsilon_iy, epsilon0_i, I_ix;
    arma::mat eta_xj, eta0_j, I_yj;

    int nbDraws_1 = build_disaggregate_epsilon(n,arums_G,epsilon_iy,epsilon0_i,I_ix);
    int nbDraws_2 = build_disaggregate_epsilon(m,arums_H,eta_xj,eta0_j,I_yj);

    eta_xj = eta_xj.t();

    epsilon0_i = arma::vectorise(epsilon0_i);
    eta0_j = arma::vectorise(eta0_j);

    I_yj = I_yj.t();
    //
    arma::vec n_i = arma::vectorise(I_ix * n) / (double) nbDraws_1;
    arma::vec m_j = arma::vectorise(m * I_yj) / (double) nbDraws_2;

    int nbI = n_i.n_elem;
    int nbJ = m_j.n_elem;
    //
    arma::mat A_11_lp = arma::kron(arma::ones(nbY,1), arma::eye(nbI,nbI));
    arma::mat A_12_lp = arma::zeros(nbI*nbY,nbJ);
    arma::mat A_13_lp = arma::kron(- arma::eye(nbY,nbY), I_ix);

    arma::mat A_1_lp = arma::join_rows(arma::join_rows(A_11_lp,A_12_lp),A_13_lp);

    arma::mat A_21_lp = arma::zeros(nbX*nbJ,nbI);
    arma::mat A_22_lp = arma::kron(arma::eye(nbJ,nbJ),arma::ones(nbX,1));
    arma::mat A_23_lp = arma::kron(I_yj.t(), arma::eye(nbX,nbX));

    arma::mat A_2_lp = arma::join_rows(arma::join_rows(A_21_lp,A_22_lp),A_23_lp);

    arma::mat A_lp = arma::join_cols(A_1_lp,A_2_lp);

    int k_lp = A_lp.n_rows;
    int n_lp = A_lp.n_cols;

    int jj;
    char* sense_lp = new char[k_lp];
    for (jj=0; jj<k_lp; jj++) {
        sense_lp[jj] = '>';
    }

    arma::vec lb_lp(epsilon0_i.n_elem + eta0_j.n_elem + nbX*nbY);
    lb_lp.rows(0,epsilon0_i.n_elem-1) = arma::vectorise(epsilon0_i);
    lb_lp.rows(epsilon0_i.n_elem,epsilon0_i.n_elem + eta0_j.n_elem - 1) = eta0_j;
    lb_lp.rows(epsilon0_i.n_elem + eta0_j.n_elem, lb_lp.n_rows - 1).fill(-arma::datum::inf);

    arma::vec rhs_lp(epsilon_iy.n_elem + eta_xj.n_elem);
    rhs_lp.rows(0,epsilon_iy.n_elem-1) = arma::vectorise(epsilon_iy);
    rhs_lp.rows(epsilon_iy.n_elem,rhs_lp.n_elem-1) = eta_xj + phi * I_yj;

    arma::vec obj_lp(nbI + nbJ + nbX*nbY);
    obj_lp.rows(0,nbI-1) = n_i;
    obj_lp.rows(nbI, nbI + nbJ-1) = m_j;
    obj_lp.rows(nbI + nbJ, obj_lp.n_elem-1).zeros();

    int modelSense = 0; // minimize
    
    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    bool lp_optimal = false;
    double val_lp = 0.0;
    //
    try {
        lp_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), A_lp.memptr(), modelSense, rhs_lp.memptr(), sense_lp, NULL, lb_lp.memptr(), NULL, NULL, val_lp, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

        if (lp_optimal) {
            U_out = arma::reshape(sol_mat(arma::span(nbI+nbJ,nbI+nbJ+nbX*nbY-1),0),nbX,nbY);
            V_out = phi - U_out;

            arma::mat mu_iy = arma::reshape(dual_mat(arma::span(0,nbI*nbY-1),0),nbI,nbY);
            mu = I_ix.t() * mu_iy;

            mux0 = n - arma::sum(mu,1);
            mu0y = m - arma::trans(arma::sum(mu,0));

            success = true;
        } else {
            std::cout << "Non-optimal value found during optimization" << std::endl;
        }
#if !defined(TRAME_USE_GUROBI_C)
    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
#endif
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //
    return success;
}

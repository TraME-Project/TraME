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
 * eap nash
 *
 * Keith O'Hara
 * 08/16/2016
 */

template<typename Ta>
bool eap_nash(dse<Ta> market, bool xFirst, double* tol_inp, arma::mat& mu, arma::vec& mux0, arma::vec& mu0y, arma::mat& u, arma::mat& v)
{
    transfers trans_obj = market.trans_obj;

    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;
    //
    arma::mat v_curr, v_next, v_err;

    if (xFirst) {
        arma::mat subdiff;
        v_curr = v_from_us(trans_obj,arma::zeros(nbX,1),NULL,&subdiff);
    } else {
        v_curr = arma::zeros(nbY,1);
    }
    //
    double tol;
    if (tol_inp) {
        tol = *tol_inp;
    } else {
        tol = 1E-12;
    }

    bool success = false;
    double err = 2*tol;
    int iter = 0;
    int max_iter = 10000;
    //
    while (err > tol && iter < max_iter) {
        iter++;

        v_next = update_v(trans_obj,v_curr,n,m,xFirst);
        v_err = arma::abs(v_next - v_curr);
        err = arma::as_scalar(arma::max(v_err));

        v_curr = v_next;
    }
    //
    v = v_curr;

    arma::mat subdiff;
    u = u_from_vs(trans_obj,v,NULL,&subdiff);

    arma::vec uv_vec = arma::join_cols(arma::vectorise(u),arma::vectorise(v)); 
    //
    arma::vec obj_lp = arma::vectorise(subdiff);

    arma::mat A_1_lp = arma::kron(arma::ones(1,nbY), arma::eye(nbX,nbX));
    arma::mat A_2_lp = arma::kron(arma::eye(nbY,nbY), arma::ones(1,nbX));
    arma::mat A_lp = arma::join_cols(A_1_lp, A_2_lp);

    arma::vec rhs_lp = arma::join_cols(n,m);

    int k_lp = A_lp.n_rows;
    int n_lp = A_lp.n_cols;

    int jj;
    char* sense_lp = new char[k_lp];
    for (jj=0; jj<k_lp; jj++) {
        if (uv_vec(jj) - 0 < tol) {
            sense_lp[jj] = '<';
        } else {
            sense_lp[jj] = '=';
        }
    }

    int modelSense = 1; // minimize

    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    bool LP_optimal = false;
    double val_lp = 0.0;
    //
    try {
        LP_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), A_lp.memptr(), modelSense, rhs_lp.memptr(), sense_lp, NULL, NULL, NULL, NULL, val_lp, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

        if (LP_optimal) {
            mu = arma::reshape(sol_mat.col(0),nbX,nbY);
            
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

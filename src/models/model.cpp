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
 * general model class
 * template specialization
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * This version:
 * 03/09/2017
 */

#include "trame.hpp"

namespace trame
{

template<>
void model<arums::logit,mmfs::tu>::build_market_TU(const arma::mat& theta)
{
    mfe_obj.build_TU(n,m,Phi_xy_theta(theta),NULL,need_norm);
}

template<>
bool model<arums::logit,mmfs::tu>::solve(arma::mat& mu_sol)
{
    bool res = mfe_obj.solve(mu_sol,NULL);
    //
    return res;
}

template<>
bool model<arums::logit,mmfs::tu>::solve(arma::mat& mu_sol, const char* solver)
{
    bool res = mfe_obj.solve(mu_sol,solver);
    //
    return res;
}

template<>
bool model<arums::logit,mmfs::tu>::solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver)
{
    bool res = mfe_obj.solve(mu_sol,U,V,solver);
    //
    return res;
}
/*
template<>
void model<arums::logit,transfers::tu>::dtheta_mu(const arma::mat& theta, const arma::mat* dtheta, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& dmu_out)
{
    build_market_TU(theta);
    //
    int range_params = theta.n_cols;
    
    arma::mat mu, U, V;
    solve(mu,U,V,NULL);
    
    arma::vec mu_x0 = mfe_obj.n - arma::sum(mu,1);
    arma::vec mu_0y = mfe_obj.m - arma::trans(arma::sum(mu,0));

    arma::mat dparams_Psi, dparams_G, dparams_H;
    dparam(dtheta,dparams_Psi,NULL,NULL);

    arma::mat du_Psi = mfe_obj.mmfs_obj.du_Psi(U,V);
    arma::mat dv_Psi = 1.0 - du_Psi;
    //
    arma::mat dtheta_psis = mfe_obj.mmfs_obj.dtheta_Psi(U,V,dparams_Psi);
    arma::vec mu_dthetapsi_vec = arma::vectorise(mu) % arma::vectorise(dtheta_psis);

    arma::cube mu_dthetapsi(mu_dthetapsi_vec.memptr(),nbX,nbY,range_params,false);

    arma::mat d_1 = cube_sum(mu_dthetapsi,0) / mfe_obj.sigma;
    arma::mat d_2 = cube_sum(mu_dthetapsi,1) / mfe_obj.sigma;

    arma::mat numer = arma::join_cols(d_1,d_2);
    //
    arma::mat Delta_11 = arma::diagmat(mu_x0 + arma::sum(mu % du_Psi,1));
    arma::mat Delta_22 = arma::diagmat(mu_0y + arma::trans(arma::sum(mu % dv_Psi,0)));

    arma::mat Delta_12 = mu % dv_Psi;
    arma::mat Delta_21 = arma::trans(mu % du_Psi);

    arma::mat Delta = arma::join_cols(arma::join_rows(Delta_11,Delta_12),arma::join_rows(Delta_21,Delta_22));
    //
    arma::mat dlogmu_singles = arma::solve(Delta,numer);
    arma::mat dlogmu_x0 = dlogmu_singles.rows(0,nbX-1);
    arma::mat dlogmu_0y = dlogmu_singles.rows(nbX,nbX+nbY-1);
}*/

template<>
bool model<arums::empirical,transfers::tu>::mme(const arma::mat& mu_hat, arma::mat& theta_hat, double* val_out, arma::mat* mu_out, arma::mat* U_out, arma::mat* V_out)
{
    bool success = false;
    //
    arma::mat kron_mat = Phi_xy();
    arma::mat kron_mat_2 = arma::reshape(kron_mat.t(),nbParams*nbX,nbY);

    arma::vec C_hat = arma::vectorise(arma::vectorise(mu_hat)*kron_mat);
    //
    arma::mat epsilon_iy, epsilon0_i, I_ix;
    arma::mat eta_xj, eta_0j, I_yj;

    int nbDraws_1 = build_disaggregate_epsilon(n,market_obj.arums_G,epsilon_iy,epsilon0_i,I_ix);
    int nbDraws_2 = build_disaggregate_epsilon(m,market_obj.arums_H,eta_xj,eta_0j,I_yj);

    epsilon0_i = arma::vectorise(epsilon0_i);

    eta_xj = eta_xj.t();
    eta_0j = arma::vectorise(eta_0j);
    I_yj = I_yj.t();
    
    arma::vec n_i = arma::vectorise(I_ix * n) / (double) nbDraws_1;
    arma::vec m_j = arma::vectorise(m.t() * I_yj) / (double) nbDraws_2;

    int nbI = n_i.n_elem;
    int nbJ = m_j.n_elem;

    arma::mat kron_data_mat = - arma::reshape(kron_mat_2*I_yj,nbParams,nbX*nbJ);
    /*
     * use batch allocation to construct the sparse constraint matrix (A)
     *
     * A_sp_t size: (nbI + nbJ + nbX*nbY + nbParams) x (nbI*nbY + nbJ*nbX)
     *
     * first block involves nbY blocks of nbI diagonal matrices, so nbI x (nbY*nbI) 
     *
     * second block begins on row nbI+1 and ends on row nbI+nbJ (inclusive); from column 
     * 1 to nbI*nbY nothing but zeros; then nbJ blocks of nbX-length row vectors of ones 
     *
     * third block begins on row nbI+nbJ+1 and ends on nbI+nbJ+nbX*nbY; 
     * there are nbX*nbY blocks of length nbDraws_1 from columns (1,nbI*nbY);
     * for columns (nbI+1)
     *
     * fourth block is filled from (nbI*nbY+1,nbI*nbY + nbJ*nbX) with kron_data_mat
     */
    int jj, kk, ll, count_val = 0;
    int num_nonzero_elem = nbI*nbY + nbJ*nbX + nbX*nbY*nbDraws_1 + nbY*nbDraws_2*nbX + nbParams*nbJ*nbX;

    arma::umat location_mat_1(2,num_nonzero_elem);
    arma::rowvec vals_mat_1(num_nonzero_elem);

    for (jj=0; jj < nbY; jj++) { // first block
        for (kk=0; kk < nbI; kk++) {
            location_mat_1(0,count_val) = kk; // rows
            location_mat_1(1,count_val) = kk + jj*nbI; // columns

            vals_mat_1(count_val) = 1.0;

            ++count_val;
        }
    }

    for (jj=0; jj < nbJ; jj++) { // second block
        for (kk=0; kk < nbX; kk++) {
            location_mat_1(0,count_val) = nbI + jj;
            location_mat_1(1,count_val) = nbI*nbY + kk + jj*nbX;

            vals_mat_1(count_val) = 1.0;

            ++count_val;
        }
    }

    for (jj=0; jj < nbX*nbY; jj++) { // third block, part 1
        for (kk=0; kk < nbDraws_1; kk++) {
            location_mat_1(0,count_val) = nbI + nbJ + jj;
            location_mat_1(1,count_val) = kk + jj*nbDraws_1;

            vals_mat_1(count_val) = -1.0;

            ++count_val;
        }
    }

    for (jj=0; jj < nbY; jj++) { // third block, part 2
        for (kk=0; kk < nbDraws_2; kk++) {
            for (ll=0; ll < nbX; ll++) {
                location_mat_1(0,count_val) = nbI + nbJ + jj*nbX + ll;
                location_mat_1(1,count_val) = nbI*nbY + jj*(nbDraws_2*nbX) + kk*nbX + ll;

                vals_mat_1(count_val) = 1.0;

                ++count_val;
            }
        }
    }

    for (jj=0; jj < nbParams; jj++) { // fourth block, data
        for (kk=0; kk < nbJ*nbX; kk++) {
            location_mat_1(0,count_val) = nbI + nbJ + nbX*nbY + jj;
            location_mat_1(1,count_val) = nbI*nbY + kk;

            vals_mat_1(count_val) = kron_data_mat(jj,kk);

            ++count_val;
        }
    }
    //
    // now proceed to solve the LP problem
    arma::sp_mat A_sp_t(location_mat_1,vals_mat_1); // transpose of A

    int k_lp = A_sp_t.n_cols; // cols as we're working with the transpose
    int n_lp = A_sp_t.n_rows; // rows as we're working with the transpose

    const arma::uword* row_vals = &(*A_sp_t.row_indices);
    const arma::uword* col_vals = &(*A_sp_t.col_ptrs);

    int* vind_lp = new int[num_nonzero_elem];
    int* vbeg_lp = new int[k_lp+1];
    double* vval_lp = new double[num_nonzero_elem];

    for (jj=0; jj<num_nonzero_elem; jj++) {
        vind_lp[jj] = row_vals[jj];
        vval_lp[jj] = A_sp_t.values[jj];
    }

    for (jj=0; jj<k_lp+1; jj++) {
        vbeg_lp[jj] = col_vals[jj];
    }
    //
    char* sense_lp = new char[k_lp];
    for (jj=0; jj<k_lp; jj++) {
        sense_lp[jj] = '>';
    }

    arma::vec lb_lp(epsilon0_i.n_elem + eta_0j.n_elem + nbX*nbY + nbParams);
    lb_lp.rows(0,epsilon0_i.n_elem-1) = arma::vectorise(epsilon0_i);
    lb_lp.rows(epsilon0_i.n_elem,epsilon0_i.n_elem + eta_0j.n_elem - 1) = eta_0j;
    lb_lp.rows(epsilon0_i.n_elem + eta_0j.n_elem, lb_lp.n_rows - 1).fill(-arma::datum::inf);

    arma::vec rhs_lp(epsilon_iy.n_elem + eta_xj.n_elem);
    rhs_lp.rows(0,epsilon_iy.n_elem-1) = arma::vectorise(epsilon_iy);
    rhs_lp.rows(epsilon_iy.n_elem,rhs_lp.n_elem-1) = arma::vectorise(eta_xj);

    arma::vec obj_lp(nbI + nbJ + nbX*nbY + C_hat.n_elem);
    obj_lp.rows(0,nbI-1) = n_i;
    obj_lp.rows(nbI, nbI + nbJ - 1) = m_j;
    obj_lp.rows(nbI + nbJ, nbI + nbJ + nbX*nbY - 1).zeros();
    obj_lp.rows(nbI + nbJ + nbX*nbY, obj_lp.n_elem-1) = - C_hat;

    int modelSense = 0; // minimize
    
    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    bool lp_optimal = false;
    double val_lp = 0.0;
    //
    try {
        lp_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), num_nonzero_elem, vbeg_lp, vind_lp, vval_lp, modelSense, rhs_lp.memptr(), sense_lp, NULL, lb_lp.memptr(), NULL, NULL, val_lp, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));
        
        if (lp_optimal) {
            arma::mat mu_iy = arma::reshape(dual_mat(arma::span(0,nbI*nbY-1),0),nbI,nbY);
            arma::mat mu = I_ix.t() * mu_iy;

            theta_hat = sol_mat(arma::span(nbI+nbJ+nbX*nbY,nbI+nbJ+nbX*nbY+nbParams-1),0);
            //
            // package up solution
            if (mu_out) {
                *mu_out = mu;
            }

            if (U_out && V_out) {
                *U_out = arma::reshape(sol_mat(arma::span(nbI+nbJ,nbI+nbJ+nbX*nbY-1),0),nbX,nbY);
                *V_out = arma::reshape(kron_mat*theta_hat,nbX,nbY) - *U_out;
            } else if (U_out) {
                *U_out = arma::reshape(sol_mat(arma::span(nbI+nbJ,nbI+nbJ+nbX*nbY-1),0),nbX,nbY);
            } else if (V_out) {
                *V_out = arma::reshape(kron_mat*theta_hat,nbX,nbY) - arma::reshape(sol_mat(arma::span(nbI+nbJ,nbI+nbJ+nbX*nbY-1),0),nbX,nbY);
            }

            if (val_out) {
                *val_out = val_lp;
            }
            //
            success = true;
        } else {
            std::cout << "Non-optimal value found during optimization" << std::endl;
        }
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //
    return success;
}

}

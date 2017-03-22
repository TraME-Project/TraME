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
 * Cupids LP solver
 *
 * Keith O'Hara
 * 08/25/2016
 *
 * This version:
 * 03/22/2017
 */

// internal cupids_lp

template<typename Tg, typename Th, typename Tm>
bool cupids_lp_int(const dse<Tg,Th,Tm>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out)
{
    bool success = false;
    //
    if (market.need_norm) {
        printf("CupidsLP does not yet allow for the case without unmatched agents.\n");
        return false;
    }
    //
    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;

    Tg arums_G = market.arums_G;
    Th arums_H = market.arums_H;

    Tm trans_obj = market.trans_obj;

    arma::mat phi = trans_obj.phi;
    //
    arma::mat epsilon_iy, epsilon_0i, I_ix;
    arma::mat eta_xj, eta_0j, I_yj;

    int nbDraws_1 = build_disaggregate_epsilon(n,arums_G,epsilon_iy,epsilon_0i,I_ix);
    int nbDraws_2 = build_disaggregate_epsilon(m,arums_H,eta_xj,eta_0j,I_yj);

    eta_xj = eta_xj.t();

    epsilon_0i = arma::vectorise(epsilon_0i);
    eta_0j = arma::vectorise(eta_0j);

    I_yj = I_yj.t();
    //
    arma::vec n_i = arma::vectorise(I_ix * n) / (double) nbDraws_1;
    arma::vec m_j = arma::vectorise(m.t() * I_yj) / (double) nbDraws_2;

    int nbI = n_i.n_elem;
    int nbJ = m_j.n_elem;
    //
    // use batch allocation to construct the sparse constraint matrix (A)
    int jj, kk, ll, count_val = 0;
    int num_non_zero = nbI*nbY + nbJ*nbX + nbDraws_1*(nbX*nbY) + nbX*nbDraws_2*nbY;

    arma::umat location_mat(2,num_non_zero);
    arma::rowvec vals_mat(num_non_zero);

    // upper diagonal blocks
    for (jj=0; jj<nbY; jj++) {
        for (kk=0; kk<nbI; kk++) {
            location_mat(0,count_val) = kk; // rows
            location_mat(1,count_val) = kk + jj*nbI; // columns

            vals_mat(count_val) = 1;

            ++count_val;
        }
    }

    // mid diagonal block (mid-right)
    for (jj=0; jj<nbJ; jj++) {
        for (kk=0; kk<nbX; kk++) {
            location_mat(0,count_val) = nbI + jj; // rows
            location_mat(1,count_val) = nbI*nbY + jj*nbX + kk; // columns

            vals_mat(count_val) = 1;

            ++count_val;
        }
    }

    // lower blocks
    for (jj=0; jj<(nbX*nbY); jj++) {
        for (kk=0; kk<nbDraws_1; kk++) {
            location_mat(0,count_val) = nbI + nbJ + jj; // rows
            location_mat(1,count_val) = jj*nbDraws_1 + kk; // columns

            vals_mat(count_val) = -1;

            ++count_val;
        }
    }

    for (jj=0; jj<nbY; jj++) {
        for (kk=0; kk<nbDraws_2; kk++) {
            for (ll=0; ll<nbX; ll++) {
                location_mat(0,count_val) = nbI + nbJ + jj*nbX + ll; // rows
                location_mat(1,count_val) = nbI*nbY + kk*nbX + jj*nbX*nbDraws_2 + ll; // columns

                vals_mat(count_val) = 1;

                ++count_val;
            }
        }
    }
    //
    // now proceed to solve the LP problem
    arma::sp_mat A_sp_t(location_mat,vals_mat); // transpose of A
    
    int k_lp = A_sp_t.n_cols; // cols as we're working with the transpose
    int n_lp = A_sp_t.n_rows; // rows as we're working with the transpose

    const arma::uword* row_vals = &(*A_sp_t.row_indices);
    const arma::uword* col_vals = &(*A_sp_t.col_ptrs);

    int* vind_lp = new int[num_non_zero];
    int* vbeg_lp = new int[k_lp+1];
    double* vval_lp = new double[num_non_zero];

    for (jj=0; jj<num_non_zero; jj++) {
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

    arma::vec lb_lp(epsilon_0i.n_elem + eta_0j.n_elem + nbX*nbY);
    lb_lp.rows(0,epsilon_0i.n_elem-1) = arma::vectorise(epsilon_0i);
    lb_lp.rows(epsilon_0i.n_elem,epsilon_0i.n_elem + eta_0j.n_elem - 1) = eta_0j;
    lb_lp.rows(epsilon_0i.n_elem + eta_0j.n_elem, lb_lp.n_rows - 1).fill(-arma::datum::inf);

    arma::vec rhs_lp(epsilon_iy.n_elem + eta_xj.n_elem);
    rhs_lp.rows(0,epsilon_iy.n_elem-1) = arma::vectorise(epsilon_iy);
    rhs_lp.rows(epsilon_iy.n_elem,rhs_lp.n_elem-1) = arma::vectorise(eta_xj + phi * I_yj);

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
        lp_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), num_non_zero, vbeg_lp, vind_lp, vval_lp, modelSense, rhs_lp.memptr(), sense_lp, NULL, lb_lp.memptr(), NULL, NULL, val_lp, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));
        
        if (lp_optimal) {
            arma::mat mu_iy = arma::reshape(dual_mat(arma::span(0,nbI*nbY-1),0),nbI,nbY);
            arma::mat mu = I_ix.t() * mu_iy;
            //
            // package up solution
            if (mu_out) {
                *mu_out = mu;
            }

            if (U_out && V_out) {
                *U_out = arma::reshape(sol_mat(arma::span(nbI+nbJ,nbI+nbJ+nbX*nbY-1),0),nbX,nbY);
                *V_out = phi - *U_out;
            } else if (U_out) {
                *U_out = arma::reshape(sol_mat(arma::span(nbI+nbJ,nbI+nbJ+nbX*nbY-1),0),nbX,nbY);
            } else if (V_out) {
                *V_out = phi - arma::reshape(sol_mat(arma::span(nbI+nbJ,nbI+nbJ+nbX*nbY-1),0),nbX,nbY);
            }

            if (mu_x0_out) {
                *mu_x0_out = n - arma::sum(mu,1);
            }
            if (mu_0y_out) {
                *mu_0y_out = m - arma::trans(arma::sum(mu,0));
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

// wrappers

template<typename Tg, typename Th, typename Tm>
bool cupids_lp(const dse<Tg,Th,Tm>& market, arma::mat& mu_out)
{
    bool res = cupids_lp_int(market,&mu_out,NULL,NULL,NULL,NULL);

    return res;
}

template<typename Tg, typename Th, typename Tm>
bool cupids_lp(const dse<Tg,Th,Tm>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    bool res = cupids_lp_int(market,&mu_out,NULL,NULL,&U_out,&V_out);
    
    return res;
}

template<typename Tg, typename Th, typename Tm>
bool cupids_lp(const dse<Tg,Th,Tm>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out)
{
    bool res = cupids_lp_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out);

    return res;
}

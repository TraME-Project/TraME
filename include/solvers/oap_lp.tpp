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
 * Optimal Assignment Problem (OAP) LP solver for TU case
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 09/23/2017
 */

// internal oap_lp

template<typename Tg, typename Th, typename Tt>
bool
oap_lp_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::vec* u_out, arma::vec* v_out, const bool* x_first_inp, double* val_out, arma::mat* residuals_out)
{
    printf("oap_lp only works for TU transfers.\n");
    return false;
}

template<typename Tg, typename Th>
bool
oap_lp_int(const dse<Tg,Th,transfers::tu>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::vec* u_out, arma::vec* v_out, const bool* x_first_inp, double* val_out, arma::mat* residuals_out)
{
    bool success = false;
    
    //

    const int nbX = market.nbX;
    const int nbY = market.nbY;
    
    // build constraint matrix

    const int num_non_zero_lp = nbX*nbY*2;

    arma::umat location_mat_lp(2,num_non_zero_lp);
    arma::rowvec vals_mat_lp = arma::ones(1,num_non_zero_lp);

    int count_val = 0;

    for (int kk=0; kk < nbY; kk++) {
        for (int jj=0; jj < nbX; jj++) {
            location_mat_lp(0,count_val) = jj + kk*nbX;
            location_mat_lp(1,count_val) = jj;
            ++count_val;
        }
        for (int jj=0; jj < nbX; jj++) {
            location_mat_lp(0,count_val) = jj + kk*nbX;
            location_mat_lp(1,count_val) = kk + nbX;
            ++count_val;
        }
    }

    arma::sp_mat A_lp_t(location_mat_lp,vals_mat_lp); // this is the transpose of the constraint matrix

    const int k_lp = A_lp_t.n_cols; // n_cols as we are working with the transpose of A
    const int n_lp = A_lp_t.n_rows; // n_rows as we are working with the transpose of A

    int* vind_lp = uword_to_int(A_lp_t.row_indices,num_non_zero_lp); // index of what row each non-zero value belongs to
    int* vbeg_lp = uword_to_int(A_lp_t.col_ptrs,k_lp+1);    // index of how many non-zero values are in each column

    double* vval_lp = new double[num_non_zero_lp];
    std::memcpy(vval_lp, A_lp_t.values, num_non_zero_lp * sizeof(double));

    //

    arma::vec obj_lp = arma::vectorise(market.trans_obj.phi);

    arma::vec rhs_lp = arma::join_cols(market.n,market.m);

    //

    char* sense_lp = new char[k_lp];
    std::memset(sense_lp, '<', k_lp * sizeof (char));

    //

    bool lp_optimal = false;
    int modelSense = 1; // maximize
    double objval, val_lp = 0.0;

    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    try {
        lp_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), num_non_zero_lp, vbeg_lp, vind_lp, vval_lp, modelSense, rhs_lp.memptr(), sense_lp, nullptr, nullptr, nullptr, nullptr, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

        if (lp_optimal) {
            const arma::mat mu = arma::reshape(sol_mat.col(0),nbX,nbY);

            if (mu_out) {
                *mu_out = mu;
            }

            if (mu_x0_out) {
                *mu_x0_out = market.n - arma::sum(mu,1);
            }
            if (mu_0y_out) {
                *mu_0y_out = market.m - arma::trans(arma::sum(mu,0));
            }

            // u_0 = dual_mat.col(0).rows(0,nbX-1);
            // v_0 = dual_mat.col(0).rows(nbX,nbX+nbY-1);

            val_lp = objval;

            if (val_out) {
                *val_out = objval;
            }
        } else {
            std::cout << "Non-optimal value found during optimization" << std::endl;
        }
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

    delete[] vbeg_lp;
    delete[] vind_lp;
    delete[] vval_lp;
    delete[] sense_lp;

    //
    // step2
    
    if ( lp_optimal && (u_out || v_out || residuals_out) ) {
        arma::vec obj_bis = arma::join_cols(market.n,-market.m);

        const bool x_first = (x_first_inp) ? *x_first_inp : true;

        if (!x_first) {
            obj_bis *= -1.0;
        }

        // build constraint matrix

        const int num_non_zero_bis = nbX*nbY*2 + nbX + nbY;

        arma::umat location_mat_bis(2,num_non_zero_bis);
        arma::rowvec vals_mat_bis(num_non_zero_bis);

        location_mat_bis.cols(0,num_non_zero_lp-1) = location_mat_lp;
        vals_mat_bis.cols(0,num_non_zero_lp-1) = vals_mat_lp;

        location_mat_bis.swap_rows(0,1);

        count_val = num_non_zero_lp;

        for (int kk=0; kk < nbX + nbY; kk++) {
            location_mat_bis(0,count_val) = kk;
            location_mat_bis(1,count_val) = nbX*nbY;

            vals_mat_bis(count_val) = rhs_lp(kk);

            ++count_val;
        }

        arma::sp_mat A_bis_t(location_mat_bis,vals_mat_bis); // this is the transpose of A_sp

        const int k_bis = A_bis_t.n_cols; // n_cols as we are working with the transpose of A
        const int n_bis = A_bis_t.n_rows; // n_rows as we are working with the transpose of A

        int* vind_bis = uword_to_int(A_bis_t.row_indices,num_non_zero_bis); // index of what row each non-zero value belongs to
        int* vbeg_bis = uword_to_int(A_bis_t.col_ptrs,k_bis+1);    // index of how many non-zero values are in each column
    
        double* vval_bis = new double[num_non_zero_bis];
        std::memcpy(vval_bis, A_bis_t.values, num_non_zero_bis * sizeof(double));

        //

        arma::vec rhs_bis(obj_lp.n_rows+1);
        rhs_bis.rows(0,rhs_bis.n_rows-2) = obj_lp;
        rhs_bis(rhs_bis.n_rows-1) = val_lp;

        //

        char* sense_bis = new char[k_bis];
        std::memset(sense_bis, '>', (k_bis-1) * sizeof (char));
    
        sense_bis[k_bis-1] = '=';

        //

        int modelSense_bis = 1; // maximize
        double objval_bis;

        arma::mat sol_mat_bis(n_bis, 2);
        arma::mat dual_mat_bis(k_bis, 2);
        
        try {
            lp_optimal = generic_LP(k_bis, n_bis, obj_bis.memptr(), num_non_zero_bis, vbeg_bis, vind_bis, vval_bis, modelSense_bis, rhs_bis.memptr(), sense_bis, nullptr, nullptr, nullptr, nullptr, objval_bis, sol_mat_bis.colptr(0), sol_mat_bis.colptr(1), dual_mat_bis.colptr(0), dual_mat_bis.colptr(1));

            if (lp_optimal) {
                const arma::mat u = sol_mat_bis.col(0).rows(0,nbX-1);
                const arma::mat v = sol_mat_bis.col(0).rows(nbX,nbX+nbY-1);

                if (u_out) {
                    *u_out = u;
                }
                if (v_out) {
                    *v_out = v;
                }

                if (residuals_out) {
                    const arma::mat u_Psi = arma::repmat(u,1,nbY);     // Keith: check use of byrow here
                    const arma::mat v_Psi = arma::repmat(v.t(),nbX,1);

                    *residuals_out = market.trans_obj.Psi(u_Psi,v_Psi);
                }
            } else {
                std::cout << "Non-optimal value found during optimization" << std::endl;
            }
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }
        //

        delete[] vbeg_bis;
        delete[] vind_bis;
        delete[] vval_bis;
        delete[] sense_bis;

    }

    //

    success = lp_optimal;

    return success;
}

// wrappers

template<typename Tg, typename Th, typename Tt>
bool
oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out)
{
    return oap_lp_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp)
{
    return oap_lp_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,&x_first_inp,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& residuals_out)
{
    return oap_lp_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr,&residuals_out);
}

template<typename Tg, typename Th, typename Tt>
bool
oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp, arma::mat& residuals_out)
{
    return oap_lp_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,&x_first_inp,nullptr,&residuals_out);
}

template<typename Tg, typename Th, typename Tt>
bool
oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& u_out, arma::vec& v_out)
{
    return oap_lp_int(market,&mu_out,nullptr,nullptr,&u_out,&v_out,nullptr,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::vec& u_out, arma::vec& v_out, const bool x_first_inp, double& val_out, arma::mat& residuals_out)
{
    return oap_lp_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&u_out,&v_out,&x_first_inp,&val_out,&residuals_out);
}

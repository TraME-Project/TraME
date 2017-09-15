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
 * Optimal Assignment Problem (OAP) LP solver for TU case
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 07/26/2017
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
    
    //

    arma::vec obj_lp = arma::vectorise(market.trans_obj.phi);

    arma::mat A_1_lp = arma::kron(arma::ones(1,nbY),arma::eye(nbX,nbX));
    arma::mat A_2_lp = arma::kron(arma::eye(nbY,nbY),arma::ones(1,nbX));
    arma::mat A_lp = arma::join_cols(A_1_lp,A_2_lp);

    arma::vec rhs_lp = arma::join_cols(market.n,market.m);

    const int k_lp = A_lp.n_rows;
    const int n_lp = A_lp.n_cols;

    char* sense_lp = new char[k_lp];
    for (int jj=0; jj<k_lp; jj++) {
        sense_lp[jj] = '<';
    }
    //
    bool LP_optimal = false;
    int modelSense = 1; // maximize
    double objval, val_lp = 0.0;

    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    try {
        LP_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), A_lp.memptr(), modelSense, rhs_lp.memptr(), sense_lp, nullptr, nullptr, nullptr, nullptr, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

        if (LP_optimal) {
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
    /*
     * step 2
     */
    if ( LP_optimal && (u_out || v_out || residuals_out) ) {
        arma::vec obj_bis = arma::join_cols(market.n,-market.m);

        const bool x_first = (x_first_inp) ? *x_first_inp : true;

        if (!x_first) {
            obj_bis *= -1.0;
        }
        //
        arma::mat A_bis = arma::join_cols(A_lp.t(),rhs_lp.t());

        arma::vec rhs_bis(obj_lp.n_rows+1);
        rhs_bis.rows(0,rhs_bis.n_rows-2) = obj_lp;
        rhs_bis(rhs_bis.n_rows-1) = val_lp;

        const int k_bis = A_bis.n_rows;
        const int n_bis = A_bis.n_cols;

        char* sense_bis = new char[k_bis];
        for (int jj=0; jj<k_bis-1; jj++) {
            sense_bis[jj] = '>';
        }
        sense_bis[k_bis-1] = '=';

        int modelSense_bis = 1; // maximize
        double objval_bis;

        arma::mat sol_mat_bis(n_bis, 2);
        arma::mat dual_mat_bis(k_bis, 2);
        //
        try {
            LP_optimal = generic_LP(k_bis, n_bis, obj_bis.memptr(), A_bis.memptr(), modelSense_bis, rhs_bis.memptr(), sense_bis, nullptr, nullptr, nullptr, nullptr, objval_bis, sol_mat_bis.colptr(0), sol_mat_bis.colptr(1), dual_mat_bis.colptr(0), dual_mat_bis.colptr(1));

            if (LP_optimal) {
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
    }
    //
    success = LP_optimal;

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
    return oap_lp_int(market,&mu_out,nullptr,nullptr,&u_out,&v_out,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
oap_lp(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::vec& u_out, arma::vec& v_out, const bool x_first_inp, double& val_out, arma::mat& residuals_out)
{
    return oap_lp_int(market,&mu_out,&x_first_inp,&mu_x0_out,&mu_0y_out,&u_out,&v_out,&val_out,&residuals_out);
}

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
 * OAP LP for TU
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/15/2017
 */

// internal oap_lp

template<typename Ta>
bool oap_lp_int(const dse<Ta,tu>& market, arma::mat* mu_out, const bool* x_first_inp, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::vec* u_out, arma::vec* v_out, double* val_out, arma::mat* residuals_out)
{
    bool success = false;
    //
    // warnings
    if (!market.TU) {
        printf("oap_lp only works for TU transfers.\n");
        return false;
    }
    //
    arma::mat phi = market.trans_obj.phi;

    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;
    //
    arma::vec obj_lp = arma::vectorise(phi);

    arma::mat A_1_lp = arma::kron(arma::ones(1,nbY),arma::eye(nbX,nbX));
    arma::mat A_2_lp = arma::kron(arma::eye(nbY,nbY),arma::ones(1,nbX));
    arma::mat A_lp = arma::join_cols(A_1_lp,A_2_lp);

    arma::vec rhs_lp = arma::join_cols(n,m);

    int k_lp = A_lp.n_rows;
    int n_lp = A_lp.n_cols;

    int jj;
    char* sense_lp = new char[k_lp];
    for (jj=0; jj<k_lp; jj++) {
        sense_lp[jj] = '<';
    }
    
    bool LP_optimal = false;
    int modelSense = 1; // maximize
    double objval; 
    
    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);
    //
    double val_lp = 0.0;
    arma::mat mu, u_0, v_0;

    try {
        LP_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), A_lp.memptr(), modelSense, rhs_lp.memptr(), sense_lp, NULL, NULL, NULL, NULL, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

        if (LP_optimal) {
            mu = arma::reshape(sol_mat.col(0),nbX,nbY);

            if (mu_out) {
                *mu_out = mu;
            }

            if (mu_x0_out) {
                *mu_x0_out = n - arma::sum(mu,1);
            }
            if (mu_0y_out) {
                *mu_0y_out = m - arma::trans(arma::sum(mu,0));
            }

            u_0 = dual_mat.col(0).rows(0,nbX-1);
            v_0 = dual_mat.col(0).rows(nbX,nbX+nbY-1);

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
        arma::vec obj_bis = arma::join_cols(n,-m);
        
        bool x_first;
        if (x_first_inp) {
            x_first = *x_first_inp;
        } else {
            x_first = true;
        }

        if (!x_first) {
            obj_bis *= -1.0;
        }
        //
        arma::mat A_bis = arma::join_cols(A_lp.t(),rhs_lp.t());

        arma::vec rhs_bis(obj_lp.n_rows+1);
        rhs_bis.rows(0,rhs_bis.n_rows-2) = obj_lp;
        rhs_bis(rhs_bis.n_rows-1) = val_lp;

        int k_bis = A_bis.n_rows;
        int n_bis = A_bis.n_cols;

        char* sense_bis = new char[k_bis];
        for (jj=0; jj<k_bis-1; jj++) {
            sense_bis[jj] = '>';
        }
        sense_bis[k_bis-1] = '=';
        
        int modelSense_bis = 1; // maximize
        double objval_bis;
        
        arma::mat sol_mat_bis(n_bis, 2);
        arma::mat dual_mat_bis(k_bis, 2);
        //
        try {
            LP_optimal = generic_LP(k_bis, n_bis, obj_bis.memptr(), A_bis.memptr(), modelSense_bis, rhs_bis.memptr(), sense_bis, NULL, NULL, NULL, NULL, objval_bis, sol_mat_bis.colptr(0), sol_mat_bis.colptr(1), dual_mat_bis.colptr(0), dual_mat_bis.colptr(1));

            if (LP_optimal) {
                arma::mat u = sol_mat_bis.col(0).rows(0,nbX-1);
                arma::mat v = sol_mat_bis.col(0).rows(nbX,nbX+nbY-1);

                if (u_out) {
                    *u_out = u;
                }
                if (v_out) {
                    *v_out = v;
                }

                if (residuals_out) {
                    arma::mat u_Psi = arma::repmat(u,1,nbY);
                    arma::mat v_Psi = arma::repmat(v.t(),nbX,1);

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

template<typename Ta, typename Tm>
bool oap_lp(const dse<Ta,Tm>& market, arma::mat& mu_out)
{
    bool res = oap_lp_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    
    return res;
}

template<typename Ta, typename Tm>
bool oap_lp(const dse<Ta,Tm>& market, arma::mat& mu_out, arma::mat& residuals_out)
{
    bool res = oap_lp_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&residuals_out);
    
    return res;
}

template<typename Ta, typename Tm>
bool oap_lp(const dse<Ta,Tm>& market, arma::mat& mu_out, const bool& x_first_inp, arma::mat& residuals_out)
{
    bool res = oap_lp_int(market,&mu_out,&x_first_inp,NULL,NULL,NULL,NULL,NULL,&residuals_out);
    
    return res;
}

template<typename Ta, typename Tm>
bool oap_lp(const dse<Ta,Tm>& market, arma::mat& mu_out, arma::vec& u_out, arma::vec& v_out)
{
    bool res = oap_lp_int(market,&mu_out,NULL,NULL,NULL,&u_out,&v_out,NULL,NULL);
    
    return res;
}

template<typename Ta, typename Tm>
bool oap_lp(const dse<Ta,Tm>& market, arma::mat& mu_out, const bool& x_first_inp, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::vec& u_out, arma::vec& v_out, double& val_out, arma::mat& residuals_out)
{
    bool res = oap_lp_int(market,&mu_out,&x_first_inp,&mu_x0_out,&mu_0y_out,&u_out,&v_out,&val_out,&residuals_out);
    
    return res;
}

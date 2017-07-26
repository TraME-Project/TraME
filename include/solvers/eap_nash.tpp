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
 * Equilibrium Assignment Problem (EAP) - Nash
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 07/26/2017
 */

// internal eap_nash

template<typename Tg, typename Th, typename Tt>
bool
eap_nash_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* u_out, arma::mat* v_out, const bool* x_first_inp, const double* tol_inp, const int* max_iter_inp)
{
    bool success = false;
    //
    const int nbX = market.nbX;
    const int nbY = market.nbY;

    const bool x_first = (x_first_inp) ? *x_first_inp : true;
    const double tol = (tol_inp) ? *tol_inp : 1E-12;
    const int max_iter = (max_iter_inp) ? *max_iter_inp : 10000;
    //
    arma::mat v_curr, v_next;

    if (x_first) {
        v_curr = v_from_us(market.trans_obj,arma::zeros(nbX,1),NULL,NULL);
    } else {
        v_curr = arma::zeros(nbY,1);
    }

    int iter = 0;
    double err = 2*tol;
    
    while (err > tol && iter < max_iter) {
        iter++;

        v_next = update_v(market.trans_obj,v_curr,market.n,market.m,x_first);
        err = elem_max(arma::abs(v_next - v_curr));

        v_curr = v_next;
    }
    //
    arma::mat subdiff;
    arma::mat u = u_from_vs(market.trans_obj,v_curr,NULL,&subdiff);

    arma::vec uv_vec = arma::join_cols(arma::vectorise(u),arma::vectorise(v_curr));

    if (u_out) {
        *u_out = u;
    }
    if (v_out) {
        *v_out = v_curr;
    }
    //
    arma::vec obj_lp = arma::vectorise(subdiff);

    arma::mat A_1_lp = arma::kron(arma::ones(1,nbY), arma::eye(nbX,nbX));
    arma::mat A_2_lp = arma::kron(arma::eye(nbY,nbY), arma::ones(1,nbX));
    arma::mat A_lp = arma::join_cols(A_1_lp, A_2_lp);

    arma::vec rhs_lp = arma::join_cols(market.n,market.m);

    const int k_lp = A_lp.n_rows;
    const int n_lp = A_lp.n_cols;

    char* sense_lp = new char[k_lp];
    for (int jj=0; jj<k_lp; jj++) {
        if (uv_vec(jj) - 0 < tol) {
            sense_lp[jj] = '<';
        } else {
            sense_lp[jj] = '=';
        }
    }

    int modelSense = 1; // maximize

    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    bool LP_optimal = false;
    double val_lp = 0.0;
    //
    try {
        LP_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), A_lp.memptr(), modelSense, rhs_lp.memptr(), sense_lp, NULL, NULL, NULL, NULL, val_lp, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

        if (LP_optimal) {
            arma::mat mu = arma::reshape(sol_mat.col(0),nbX,nbY);
            //
            if (mu_out) {
                *mu_out = mu;
            }
            
            if (mu_x0_out) {
                *mu_x0_out = market.n - arma::sum(mu,1);
            }
            if (mu_0y_out) {
                *mu_0y_out = market.m - arma::trans(arma::sum(mu,0));
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

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out)
{
    return eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp)
{
    return eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,&x_first_inp,NULL,NULL);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double tol_inp)
{
    return eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,NULL);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const int max_iter_inp)
{
    return eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp, const double tol_inp)
{
    return eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,&x_first_inp,&tol_inp,NULL);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp, const int max_iter_inp)
{
    return eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double tol_inp, const int max_iter_inp)
{
    return eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,&tol_inp,&max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp, const double tol_inp, const int max_iter_inp)
{
    return eap_nash_int(market,&mu_out,NULL,NULL,NULL,NULL,&x_first_inp,&tol_inp,&max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& u_out, arma::mat& v_out)
{
    return eap_nash_int(market,&mu_out,NULL,NULL,&u_out,&v_out,NULL,NULL,NULL);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::vec& u_out, arma::vec& v_out, const bool* x_first_inp, const double* tol_inp, const int* max_iter_inp)
{
    return eap_nash_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&u_out,&v_out,x_first_inp,tol_inp,max_iter_inp);
}

// internal functions

template<typename Tt>
arma::mat
u_from_vs(const Tt& trans_obj, const arma::mat& v, const double* tol_inp, arma::mat* subdiff)
{
    const arma::mat us = trans_obj.Ucal(v,NULL,NULL);
    arma::mat u  = elem_max(arma::max(us,1),0.0);
    //
    if (subdiff) {
        const double tol = (tol_inp) ? *tol_inp : 0.0;

        *subdiff = arma::zeros(trans_obj.nbX,trans_obj.nbY);
        subdiff->elem( arma::find(arma::abs(elem_sub(u,us)) <= tol) ).ones();
    }
    //
    return u;
}

template<typename Tt>
arma::mat
v_from_us(const Tt& trans_obj, const arma::mat& u, const double* tol_inp, arma::mat* subdiff)
{
    const arma::mat vs = trans_obj.Vcal(u,NULL,NULL);
    arma::mat v  = arma::trans(elem_max(arma::max(vs,0),0.0));
    //
    if (subdiff) {
        const double tol = (tol_inp) ? *tol_inp : 0.0;

        *subdiff = arma::zeros(trans_obj.nbY,trans_obj.nbX);
        subdiff->elem( arma::find(arma::abs(elem_sub(v,vs.t())) <= tol) ).ones();
        *subdiff = subdiff->t();
    }
    //
    return v;
}

template<typename Tt>
arma::mat
update_v(const Tt& trans_obj, const arma::mat& v, const arma::vec& n, const arma::vec& m, const bool x_first)
{
    const int nbX = trans_obj.nbX;
    const int nbY = trans_obj.nbY;
    
    arma::mat the_mat = arma::zeros(nbX,nbY);
    arma::mat v_updated = arma::zeros(nbY,1);
    //
    // LP setup
    int jj, y, x, yp;

    arma::vec obj_lp(nbX+1,1);
    obj_lp.rows(0,nbX-1) = n;
    obj_lp.row(nbX) = m(0);

    arma::vec rhs_lp = arma::zeros(nbX,1);
    arma::mat A_lp = arma::join_rows(arma::eye(nbX,nbX),arma::ones(nbX,1));
    arma::mat lb_lp = arma::zeros(nbX+1,1);

    const int k_lp = A_lp.n_rows;
    const int n_lp = A_lp.n_cols;

    char* sense_lp = new char[k_lp];
    for (jj=0; jj<k_lp; jj++) {
        sense_lp[jj] = '>';
    }

    bool LP_optimal = false;
    int modelSense = 0; // minimize
    double objval;

    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    double val_lp = 0.0;
    //
    arma::vec obj_bis = arma::zeros(nbX+1,1);
    obj_bis(nbX,0) = 1.0;

    arma::vec rhs_bis = arma::zeros(nbX+1,1);
    rhs_bis.rows(0,nbX-1) = rhs_lp;

    arma::mat A_bis(nbX+1,nbX+1);
    A_bis.rows(0,nbX-1) = A_lp;
    A_bis.row(nbX) = obj_lp.t();

    arma::mat lb_bis = arma::zeros(nbX+1,1);

    const int k_bis = A_bis.n_rows;
    const int n_bis = A_bis.n_cols;

    char* sense_bis = new char[k_bis];
    for (jj=0; jj<k_bis-1; jj++) {
        sense_bis[jj] = '>';
    }
    sense_bis[k_bis-1] = '=';

    int modelSense_bis = (x_first) ? 0 : 1;

    arma::mat sol_mat_bis(n_bis, 2);
    arma::mat dual_mat_bis(k_bis, 2);

    double val_bis = 0.0;
    //
    double temp_u;
    arma::mat u;

    for (y=0; y<nbY; y++) {
        for (x=0; x<nbX; x++) {
            for (yp=0; yp<nbY; yp++) {
                if (yp==y) {
                    the_mat(x,yp) = trans_obj.Vcal(0.0,x,y);
                } else {
                    temp_u = trans_obj.Ucal(v(yp),x,yp);
                    the_mat(x,yp) = trans_obj.Vcal(temp_u,x,y);
                }
            } // end of yp loop
        } // end of x loop
        //

        lb_lp.rows(0,nbX-1) = - arma::min(the_mat,1);
        obj_lp(nbX,0) = m(y);

        try {
            LP_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), A_lp.memptr(), modelSense, rhs_lp.memptr(), sense_lp, NULL, lb_lp.memptr(), NULL, NULL, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

            if (LP_optimal) {
                //u0  = sol_mat.col(0).rows(0,nbX-1);
                //v0y = sol_mat(nbX,0);

                val_lp = objval;
            } else {
                std::cout << "Non-optimal value found during optimization" << std::endl;
            }
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }

        //
        // second LP
        //

        A_bis(nbX,nbX) = m(y);
        rhs_bis(nbX) = val_lp;

        try {
            LP_optimal = generic_LP(k_bis, n_bis, obj_bis.memptr(), A_bis.memptr(), modelSense_bis, rhs_bis.memptr(), sense_bis, NULL, lb_lp.memptr(), NULL, NULL, val_bis, sol_mat_bis.colptr(0), sol_mat_bis.colptr(1), dual_mat_bis.colptr(0), dual_mat_bis.colptr(1));

            if (LP_optimal) {
                u = sol_mat.col(0).rows(0,nbX-1);
                v_updated(y,0) = sol_mat(nbX,0);
            } else {
                std::cout << "Non-optimal value found during optimization" << std::endl;
            }
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }
    }
    //
    return v_updated;
}

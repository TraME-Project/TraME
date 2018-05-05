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
 * Equilibrium Assignment Problem (EAP) - Nash
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/04/2018
 */

// internal eap_nash

template<typename Tg, typename Th, typename Tt>
bool
eap_nash_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* u_out, arma::mat* v_out,
             const bool x_first, const double err_tol, const uint_t max_iter)
{
    bool success = false;

    //

    const uint_t nbX = market.nbX;
    const uint_t nbY = market.nbY;

    //

    arma::mat v_curr;

    if (x_first) {
        v_curr = v_from_us(market.transfers_obj,arma::zeros(nbX,1),nullptr,nullptr);
    } else {
        v_curr = arma::zeros(nbY,1); // Keith: should this be u_from_vs?
    }

    uint_t iter = 0;
    double err = 2*err_tol;
    arma::mat v_next;

    while (err > err_tol && iter < max_iter)
    {
        iter++;

        v_next = update_v(market.transfers_obj,v_curr,market.n,market.m,x_first);
        err = elem_max(arma::abs(v_next - v_curr));

        v_curr = std::move(v_next);
    }

    //

    arma::mat subdiff;
    arma::mat u = u_from_vs(market.transfers_obj,v_curr,nullptr,&subdiff);

    arma::vec uv_vec = arma::join_cols(arma::vectorise(u),arma::vectorise(v_curr));

    if (u_out) {
        *u_out = u;
    }
    if (v_out) {
        *v_out = v_curr;
    }

    // build constraint matrix

    const uint_t num_non_zero_lp = nbX*nbY*2;

    arma::umat location_mat_lp(2,num_non_zero_lp);
    arma::rowvec vals_mat_lp = arma::ones(1,num_non_zero_lp);

    uint_t count_val = 0;

    for (uint_t kk=0; kk < nbY; kk++) {
        for (uint_t jj=0; jj < nbX; jj++) {
            location_mat_lp(0,count_val) = jj + kk*nbX;
            location_mat_lp(1,count_val) = jj;
            ++count_val;
        }
        for (uint_t jj=0; jj < nbX; jj++) {
            location_mat_lp(0,count_val) = jj + kk*nbX;
            location_mat_lp(1,count_val) = kk + nbX;
            ++count_val;
        }
    }

    arma::sp_mat A_lp_t(location_mat_lp,vals_mat_lp); // this is the transpose of the constraint matrix

    const uint_t k_lp = A_lp_t.n_cols; // n_cols as we are working with the transpose of A
    const uint_t n_lp = A_lp_t.n_rows; // n_rows as we are working with the transpose of A

    int* vind_lp = uword_to_int(A_lp_t.row_indices,num_non_zero_lp); // index of what row each non-zero value belongs to
    int* vbeg_lp = uword_to_int(A_lp_t.col_ptrs,k_lp+1);    // index of how many non-zero values are in each column

    double* vval_lp = new double[num_non_zero_lp];
    std::memcpy(vval_lp, A_lp_t.values, num_non_zero_lp * sizeof(double));

    //

    arma::vec obj_lp = std::move(arma::vectorise(subdiff));

    arma::vec rhs_lp = arma::join_cols(market.n,market.m);

    //

    char* sense_lp = new char[k_lp];
    for (uint_t jj=0; jj<k_lp; jj++)
    {
        if (uv_vec(jj) - 0 < err_tol) {
            sense_lp[jj] = '<';
        } else {
            sense_lp[jj] = '=';
        }
    }

    //

    bool lp_optimal = false;
    int modelSense = 1; // maximize
    double val_lp = 0.0;

    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    try {
        lp_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), num_non_zero_lp, vbeg_lp, vind_lp, vval_lp, modelSense, rhs_lp.memptr(), sense_lp, nullptr, nullptr, nullptr, nullptr, val_lp, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

        if (lp_optimal) {

            // arma::mat mu = arma::reshape(sol_mat.col(0),nbX,nbY);
            arma::mat mu(sol_mat.colptr(0),nbX,nbY,false,true);
            
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
    return eap_nash_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp)
{
    return eap_nash_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,x_first_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double err_tol_inp, const uint_t max_iter_inp)
{
    return eap_nash_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,true,err_tol_inp,max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const bool x_first_inp, const double err_tol_inp, const uint_t max_iter_inp)
{
    return eap_nash_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,x_first_inp,err_tol_inp,max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& u_out, arma::mat& v_out)
{
    return eap_nash_int(market,&mu_out,nullptr,nullptr,&u_out,&v_out);
}

template<typename Tg, typename Th, typename Tt>
bool
eap_nash(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& u_out, arma::mat& v_out,
         const bool x_first_inp, const double err_tol_inp, const uint_t max_iter_inp)
{
    return eap_nash_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&u_out,&v_out,x_first_inp,err_tol_inp,max_iter_inp);
}

// internal functions

template<typename Tt>
arma::mat
u_from_vs(const Tt& transfers_obj, const arma::mat& v, const double* tol_inp, arma::mat* subdiff)
{
    const arma::mat us = transfers_obj.Ucal(v,nullptr,nullptr);
    arma::mat u  = elem_max(arma::max(us,1),0.0);
    //
    if (subdiff) {
        const double tol = (tol_inp) ? *tol_inp : 0.0;

        *subdiff = arma::zeros(transfers_obj.nbX,transfers_obj.nbY);
        subdiff->elem( arma::find(arma::abs(elem_sub(u,us)) <= tol) ).ones();
    }
    //
    return u;
}

template<typename Tt>
arma::mat
v_from_us(const Tt& transfers_obj, const arma::mat& u, const double* tol_inp, arma::mat* subdiff)
{
    const arma::mat vs = transfers_obj.Vcal(u,nullptr,nullptr);
    arma::mat v  = arma::trans(elem_max(arma::max(vs,0),0.0));
    //
    if (subdiff) {
        const double tol = (tol_inp) ? *tol_inp : 0.0;

        *subdiff = arma::zeros(transfers_obj.nbY,transfers_obj.nbX);
        subdiff->elem( arma::find(arma::abs(elem_sub(v,vs.t())) <= tol) ).ones();
        *subdiff = subdiff->t();
    }
    //
    return v;
}

template<typename Tt>
arma::mat
update_v(const Tt& transfers_obj, const arma::mat& v, const arma::vec& n, const arma::vec& m, const bool x_first)
{
    const uint_t nbX = transfers_obj.nbX;
    const uint_t nbY = transfers_obj.nbY;

    arma::mat the_mat = arma::zeros(nbX,nbY);
    arma::mat v_updated = arma::zeros(nbY,1);

    //
    // LP setup

    const uint_t num_non_zero_lp = 2*nbX;

    arma::umat location_mat_lp(2,num_non_zero_lp);
    arma::rowvec vals_mat_lp = arma::ones(1,num_non_zero_lp);

    uint_t count_val = 0;

    for (uint_t kk=0; kk < nbX; kk++) {
        location_mat_lp(0,count_val) = kk;
        location_mat_lp(1,count_val) = kk;
        ++count_val;

        location_mat_lp(0,count_val) = nbX;
        location_mat_lp(1,count_val) = kk;

        ++count_val;
    }

    arma::sp_mat A_lp_t(location_mat_lp,vals_mat_lp); // this is the transpose of the constraint matrix

    const uint_t k_lp = A_lp_t.n_cols; // n_cols as we are working with the transpose of A
    const uint_t n_lp = A_lp_t.n_rows; // n_rows as we are working with the transpose of A

    int* vind_lp = uword_to_int(A_lp_t.row_indices,num_non_zero_lp); // index of what row each non-zero value belongs to
    int* vbeg_lp = uword_to_int(A_lp_t.col_ptrs,k_lp+1);    // index of how many non-zero values are in each column

    double* vval_lp = new double[num_non_zero_lp];
    std::memcpy(vval_lp, A_lp_t.values, num_non_zero_lp * sizeof(double));

    //

    arma::vec obj_lp(nbX+1,1);
    obj_lp.rows(0,nbX-1) = n;
    obj_lp.row(nbX) = m(0);

    arma::vec rhs_lp = arma::zeros(nbX,1);
    arma::mat lb_lp = arma::zeros(nbX+1,1);

    //

    char* sense_lp = new char[k_lp];
    std::memset(sense_lp, '>', k_lp * sizeof (char));

    //

    bool lp_optimal = false;
    int modelSense_lp = 0; // minimize
    double val_lp = 0.0;

    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    //
    // BIS setup

    const uint_t num_non_zero_bis = 2*nbX + nbX + 1;

    arma::umat location_mat_bis(2,num_non_zero_bis);
    arma::rowvec vals_mat_bis = arma::ones(1,num_non_zero_bis);

    location_mat_bis.cols(0,num_non_zero_lp-1) = location_mat_lp;

    vals_mat_bis.cols(0,num_non_zero_lp-1) = vals_mat_lp;
    vals_mat_bis.cols(num_non_zero_lp,num_non_zero_bis-2) = n.t();
    vals_mat_bis.col(num_non_zero_bis-1) = m(0);

    count_val = num_non_zero_lp;

    for (uint_t kk=0; kk < nbX + 1; kk++) {
        location_mat_bis(0,count_val) = kk;
        location_mat_bis(1,count_val) = nbX;
        ++count_val;
    }

    arma::sp_mat A_bis_t(location_mat_bis,vals_mat_bis); // this is the transpose of the constraint matrix

    const uint_t k_bis = A_bis_t.n_cols; // n_cols as we are working with the transpose of A
    const uint_t n_bis = A_bis_t.n_rows; // n_rows as we are working with the transpose of A

    int* vind_bis = uword_to_int(A_bis_t.row_indices,num_non_zero_bis); // index of what row each non-zero value belongs to
    int* vbeg_bis = uword_to_int(A_bis_t.col_ptrs,k_bis+1);    // index of how many non-zero values are in each column

    double* vval_bis = new double[num_non_zero_bis];
    std::memcpy(vval_bis, A_bis_t.values, num_non_zero_bis * sizeof(double));

    //

    arma::vec obj_bis = arma::zeros(nbX+1,1);
    obj_bis(nbX,0) = 1.0;

    arma::vec rhs_bis = arma::zeros(nbX+1,1);
    rhs_bis.rows(0,nbX-1) = rhs_lp;

    //

    char* sense_bis = new char[k_bis];
    std::memset(sense_bis, '>', (k_bis-1) * sizeof (char));
    
    sense_bis[k_bis-1] = '=';

    //

    int modelSense_bis = (x_first) ? 0 : 1;

    arma::mat sol_mat_bis(n_bis, 2);
    arma::mat dual_mat_bis(k_bis, 2);

    double val_bis = 0.0;

    // begin loop

    for (uint_t y=0; y<nbY; y++) {
        for (uint_t x=0; x<nbX; x++) {
            for (uint_t yp=0; yp<nbY; yp++) {
                if (yp==y) {
                    the_mat(x,yp) = transfers_obj.Vcal(0.0,x,y);
                } else {
                    double temp_u = transfers_obj.Ucal(v(yp),x,yp);
                    the_mat(x,yp) = transfers_obj.Vcal(temp_u,x,y);
                }
            } // end of yp loop
        } // end of x loop

        lb_lp.rows(0,nbX-1) = - arma::min(the_mat,1);
        obj_lp(nbX,0) = m(y);

        try {
            lp_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), num_non_zero_lp, vbeg_lp, vind_lp, vval_lp, modelSense_lp, rhs_lp.memptr(), sense_lp, nullptr, lb_lp.memptr(), nullptr, nullptr, val_lp, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

            if (!lp_optimal) {
                std::cout << "Non-optimal value found during optimization" << std::endl;
            }
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }

        //
        // second LP

        vval_bis[num_non_zero_bis-1] = m(y);
        rhs_bis(nbX) = val_lp;

        try {
            lp_optimal = generic_LP(k_bis, n_bis, obj_bis.memptr(), num_non_zero_bis, vbeg_bis, vind_bis, vval_bis, modelSense_bis, rhs_bis.memptr(), sense_bis, nullptr, lb_lp.memptr(), nullptr, nullptr, val_bis, sol_mat_bis.colptr(0), sol_mat_bis.colptr(1), dual_mat_bis.colptr(0), dual_mat_bis.colptr(1));

            if (lp_optimal) {
                v_updated(y,0) = sol_mat(nbX,0);
            } else {
                std::cout << "Non-optimal value found during optimization" << std::endl;
            }
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }
    }

    //

    delete[] vbeg_lp;
    delete[] vind_lp;
    delete[] vval_lp;
    delete[] sense_lp;

    delete[] vbeg_bis;
    delete[] vind_bis;
    delete[] vval_bis;
    delete[] sense_bis;

    //

    return v_updated;
}

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
 * empirical additive random utility model (ARUM) class
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 07/03/2017
 */

#include "trame.hpp"

//
// build functions

trame::arums::empirical::empirical(const int nbX_inp, const int nbY_inp)
{
    this->build(nbX_inp, nbY_inp);
}

trame::arums::empirical::empirical(const int nbX_inp, const int nbY_inp, const arma::cube& atoms_inp, const bool x_homogeneous_inp, const bool outside_option_inp)
{
    this->build(nbX_inp, nbY_inp, atoms_inp, x_homogeneous_inp, outside_option_inp);
}

void
trame::arums::empirical::build(const int nbX_inp, const int nbY_inp)
{
    nbX = nbX_inp;
    nbY = nbY_inp;
}

void
trame::arums::empirical::build(const int nbX_inp, const int nbY_inp, const arma::cube& atoms_inp, const bool x_homogeneous_inp, const bool outside_option_inp)
{
    nbX = nbX_inp;
    nbY = nbY_inp;

    atoms = atoms_inp;

    dim_params = atoms_inp.n_elem;
    aux_nbDraws = atoms.n_rows;

    x_homogeneous = x_homogeneous_inp;
    outside_option = outside_option_inp;

    nb_options = (outside_option_inp) ? nbY_inp + 1 : nbY_inp;

    presolve_LP_Gstar();
    presolve_LP_Gbar();
}

//
// indirect utility

double
trame::arums::empirical::G(const arma::vec& n)
{
    return this->G(n,U,mu_sol);
}

double
trame::arums::empirical::G(const arma::vec& n, const arma::mat& U_inp, arma::mat& mu_out)
const
{
    double val=0.0, val_x;
    mu_out.set_size(nbX,nbY);
    //
    arma::mat mu_x_temp;

    for (int i=0; i<nbX; i++) {
        val_x = Gx(U_inp.row(i).t(), mu_x_temp, i);

        val += n(i)*val_x;
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double
trame::arums::empirical::Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, const int x)
const
{
    mu_x_out.set_size(nbY,1);
    //
    arma::mat U_xs = (outside_option) ? arma::join_cols(arma::vectorise(U_x_inp),arma::zeros(1,1)) : U_x_inp;
    const arma::mat Utilde = (x_homogeneous) ? arma::ones(aux_nbDraws,1) * U_xs.t() + atoms.slice(0) : arma::ones(aux_nbDraws,1) * U_xs.t() + atoms.slice(x);
    //
    const arma::vec argmaxs = arma::max(Utilde,1);       // take max over dim = 1
    const arma::uvec argmax_inds = which_max(Utilde,1);

    double val_x = arma::accu(argmaxs) / (double)(aux_nbDraws);
    //
    arma::uvec temp_find;

    for (int tt=0; tt < nbY; tt++) {
        temp_find = arma::find(argmax_inds==tt);
        mu_x_out(tt,0) = (double)(temp_find.n_elem)/(double)(aux_nbDraws);
    }
    //
    return val_x;
}

//
// Fenchel transform of G

double
trame::arums::empirical::Gstar(const arma::vec& n)
{
    return this->Gstar(n,mu_sol,U_sol);
}

double
trame::arums::empirical::Gstar(const arma::vec& n, const arma::mat& mu_inp, arma::mat& U_out)
const
{
    if (!TRAME_PRESOLVED_GSTAR) {
        printf("TraME: Gstar cannot be called before first running presolve.\n");
    }
    //
    double val = 0.0, val_x = 0.0;
    U_out.set_size(nbX,nbY);
    //
    arma::mat U_x;

    for (int i=0; i < nbX; i++) {
        val_x = Gstarx((mu_inp.row(i).t())/n(i),U_x,i);

        val += n(i)*val_x;
        U_out.row(i) = arma::trans(U_x);
    }
    //
    return val;
}

double
trame::arums::empirical::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x)
const
{
    if (!TRAME_PRESOLVED_GSTAR) {
        printf("TraME: Gstarx cannot be called before first running presolve.\n");
    }
    //
    const arma::mat Phi = (x_homogeneous) ? atoms.slice(0) : atoms.slice(x);
    //
    arma::vec p = arma::ones(aux_nbDraws,1)/aux_nbDraws;
    arma::mat q;

    if (outside_option) {
        arma::mat temp_q(1,1);
        temp_q(0,0) = 1 - arma::accu(mu_x_inp);

        q = arma::join_cols(arma::vectorise(mu_x_inp),temp_q);
    } else {
        q = arma::vectorise(mu_x_inp);
    }
    //
    arma::vec obj_lp = arma::vectorise(Phi);

    arma::vec rhs_lp = arma::join_cols(p,q);

    char* sense_lp = new char[k_Gstar];
    for (int jj=0; jj < k_Gstar; jj++) {
        sense_lp[jj] = '=';
    }

    bool LP_optimal;
    int modelSense = 1; // maximize
    double objval, val_x = 0.0;

    arma::mat sol_mat(n_Gstar, 2);
    arma::mat dual_mat(k_Gstar, 2);
    try {
        LP_optimal = generic_LP(k_Gstar, n_Gstar, obj_lp.memptr(), numnz_Gstar, vbeg_Gstar, vind_Gstar, vval_Gstar, modelSense, rhs_lp.memptr(), sense_lp, NULL, NULL, NULL, NULL, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));
        //
        if (LP_optimal) {
            arma::mat u = dual_mat.col(0).rows(0,aux_nbDraws-1);

            if (outside_option) {
                const arma::mat U_x_temp = dual_mat.col(0).rows(aux_nbDraws,aux_nbDraws+nbY);
                U_x_out = - U_x_temp.rows(0,nbY-1) + arma::as_scalar(U_x_temp.row(nbY));
            } else {
                const arma::mat U_x_temp = dual_mat.col(0).rows(aux_nbDraws,aux_nbDraws+nbY-1);
                U_x_out = - U_x_temp - arma::accu(p % u);
            }
            //
            val_x = -objval;
        } else {
            std::cout << "Non-optimal value found during optimization" << std::endl;
        }
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //
    delete[] sense_lp;
    //
    return val_x;
}

//
// Gbar is used by DARUM

double
trame::arums::empirical::Gbar(const arma::mat& Ubar, const arma::mat& mubar, const arma::vec& n, arma::mat& U_out, arma::mat& mu_out)
const
{
    if (!TRAME_PRESOLVED_GBAR) {
        printf("TraME: Gbar cannot be called before first running presolve.\n");
    }
    //
    double val = 0.0;

    U_out.set_size(nbX,nbY);
    mu_out.set_size(nbX,nbY);
    arma::mat U_x_temp, mu_x_temp;
    //
    for (int i=0; i<nbX; i++) {
        const double val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),U_x_temp,mu_x_temp,i);

        val += n(i)*val_temp;
        U_out.row(i) = arma::trans(U_x_temp);
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double
trame::arums::empirical::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, int x)
const
{
    if (!TRAME_PRESOLVED_GBAR) {
        printf("TraME: Gbarx cannot be called before first running presolve.\n");
    }
    //
    arma::mat Phi, U_x_temp;

    if (!outside_option) {
        printf("Gbarx not implemented for empirical with outside_option = false\n");
        return 0;
    }

    if (x_homogeneous) {
        Phi = atoms.slice(0);
    } else {
        Phi = atoms.slice(x);
    }
    //
    arma::vec obj_lp_1 = mubar_x;
    arma::vec obj_lp_2 = - arma::ones(aux_nbDraws,1)/aux_nbDraws;
    arma::vec obj_lp   = arma::join_cols(obj_lp_1,obj_lp_2);

    arma::vec rhs_lp_1 = Ubar_x;
    arma::vec rhs_lp_2 = arma::vectorise(-Phi);
    arma::vec rhs_lp = arma::join_cols(rhs_lp_1,rhs_lp_2);

    char* sense_lp = new char[k_Gbar];
    for (int jj=0; jj<k_Gbar; jj++) {
        sense_lp[jj] = '<';
    }

    bool LP_optimal;
    int modelSense = 1; // maximize
    double objval;

    arma::mat sol_mat(n_Gbar,2);
    arma::mat dual_mat(k_Gbar,2);

    double val_x = 0.0;

    try {
        LP_optimal = generic_LP(k_Gbar, n_Gbar, obj_lp.memptr(), numnz_Gbar, vbeg_Gbar, vind_Gbar, vval_Gbar, modelSense, rhs_lp.memptr(), sense_lp, NULL, NULL, NULL, NULL, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));
        //
        if (LP_optimal) {
            U_x_out = sol_mat.col(0).rows(0,nbY-1);
            arma::vec delta_mu_x = dual_mat.col(0).rows(0,nbY-1);
            mu_x_out = mubar_x - delta_mu_x;
            //
            val_x = arma::accu(mubar_x % Ubar_x) - objval;
        } else {
            std::cout << "Non-optimal value found during optimization" << std::endl;
        }
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //
    delete[] sense_lp;
    //
    return val_x;
}

//
// presolve functions for Gstar and Gbar

void
trame::arums::empirical::presolve_LP_Gstar()
{
    /*
     * Here we build and store the linear constraint matrix ('A') used in Gstarx.
     *
     * (Note: we actually construct A-transpose because Armadillo sparse matrices
     * are stored in compressed sparse column (CSC) format, whereas  Gurobi 
     * requires inputs to be in compressed sparse row (CSR) format.)
     * 
     * We use batch allocation because this is *much* faster than first 
     * constructing a sparse matrix and then (ex-post) inserting values.
     */

    int jj = 0, count_val=0;

    arma::umat location_mat(2,aux_nbDraws*nb_options*2);
    arma::rowvec vals_mat(aux_nbDraws*nb_options*2);

    for (int kk=0; kk < nb_options; kk++) {
        for (jj=0; jj < aux_nbDraws; jj++) {
            location_mat(0,count_val) = jj + kk*aux_nbDraws;
            location_mat(1,count_val) = jj;
            ++count_val;
        }
        for (jj=0; jj < aux_nbDraws; jj++) {
            location_mat(0,count_val) = jj + kk*aux_nbDraws;
            location_mat(1,count_val) = kk + aux_nbDraws;
            ++count_val;
        }
    }

    vals_mat.fill(1);

    arma::sp_mat A_sp_Gstar_t(location_mat,vals_mat); // this is the transpose of A_sp_Gstar

    k_Gstar = A_sp_Gstar_t.n_cols; // cols as we're working with the transpose
    n_Gstar = A_sp_Gstar_t.n_rows; // rows as we're working with the transpose

    numnz_Gstar = aux_nbDraws*nb_options*2;

    const arma::uword* row_vals = &(*A_sp_Gstar_t.row_indices);
    const arma::uword* col_vals = &(*A_sp_Gstar_t.col_ptrs);

    vind_Gstar = new int[numnz_Gstar];    // index of what row each non-zero value belongs to
    vbeg_Gstar = new int[k_Gstar+1];      // index of how many non-zero values are in each column
    vval_Gstar = new double[numnz_Gstar]; // vals

    for (jj=0; jj<numnz_Gstar; jj++) {
        vind_Gstar[jj] = row_vals[jj];
        vval_Gstar[jj] = A_sp_Gstar_t.values[jj];
    }

    for (jj=0; jj<k_Gstar+1; jj++) {
        vbeg_Gstar[jj] = col_vals[jj];
    }
    //
    TRAME_PRESOLVED_GSTAR = true;
}

void
trame::arums::empirical::presolve_LP_Gbar()
{
    /*
     * Here we build and store the linear constraint matrix ('A') used in Gbarx.
     *
     * (Note: we actually construct A-transpose because Armadillo sparse matrices
     * are stored in compressed sparse column (CSC) format, whereas  Gurobi 
     * requires inputs to be in compressed sparse row (CSR) format.)
     * 
     * We use batch allocation because this is *much* faster than first 
     * constructing a sparse matrix and then (ex-post) inserting values.
     */

    int jj = 0, count_val=0;

    arma::umat location_mat_2(2,nbY + aux_nbDraws*nbY + aux_nbDraws*(nbY+1));
    arma::rowvec vals_mat_2(nbY + aux_nbDraws*nbY + aux_nbDraws*(nbY+1));

    for (jj=0; jj < nbY; jj++) { // top-left diagonal block
        location_mat_2(0,count_val) = jj;
        location_mat_2(1,count_val) = jj;

        vals_mat_2(count_val) = 1;

        ++count_val;
    }

    for (int kk=0; kk < (nbY+1); kk++) {
        if (kk < nbY) { // top section 
            for (jj=0; jj<aux_nbDraws; jj++) {
                location_mat_2(0,count_val) = kk;
                location_mat_2(1,count_val) = nbY + jj + kk*aux_nbDraws;

                vals_mat_2(count_val) = 1;

                ++count_val;
            }
        }

        for (jj=0; jj<aux_nbDraws; jj++) { // diagonal terms (nbY+1 number of blocks)
            location_mat_2(0,count_val) = nbY + jj;
            location_mat_2(1,count_val) = nbY + jj + kk*aux_nbDraws;

            vals_mat_2(count_val) = -1;

            ++count_val;
        }
    }

    arma::sp_mat A_sp_Gbar_t(location_mat_2,vals_mat_2);

    k_Gbar = A_sp_Gbar_t.n_cols; // cols as we're working with the transpose
    n_Gbar = A_sp_Gbar_t.n_rows; // rows as we're working with the transpose

    numnz_Gbar = vals_mat_2.n_elem;

    const arma::uword* row_vals_2 = &(*A_sp_Gbar_t.row_indices);
    const arma::uword* col_vals_2 = &(*A_sp_Gbar_t.col_ptrs);

    vind_Gbar = new int[numnz_Gbar];
    vbeg_Gbar = new int[k_Gbar+1];
    vval_Gbar = new double[numnz_Gbar];

    for (jj=0; jj<numnz_Gbar; jj++) {
        vind_Gbar[jj] = row_vals_2[jj];
        vval_Gbar[jj] = A_sp_Gbar_t.values[jj];
    }

    for (jj=0; jj<k_Gbar+1; jj++) {
        vbeg_Gbar[jj] = col_vals_2[jj];
    }
    //
    TRAME_PRESOLVED_GBAR = true;
}

/*
void empirical::presolve_LP()
{
    int jj;
    //
    // Gstar
    arma::sp_mat A_sp_Gstar = arma::join_cols(kron_sp(arma::ones(1,nb_options),arma::speye(aux_nbDraws,aux_nbDraws)),kron_sp(arma::speye(nb_options,nb_options),arma::ones(1,aux_nbDraws)));

    k_Gstar = A_sp_Gstar.n_rows;
    n_Gstar = A_sp_Gstar.n_cols;

    arma::sp_mat A_sp_Gstar_t = arma::trans(A_sp_Gstar); // need to transpose to get data into CSR format (not CSC)

    numnz_Gstar = nonzeros(A_sp_Gstar).n_elem;

    const arma::uword* row_vals = &(*A_sp_Gstar_t.row_indices);
    const arma::uword* col_vals = &(*A_sp_Gstar_t.col_ptrs);

    vind_Gstar = new int[numnz_Gstar];
    vbeg_Gstar = new int[k_Gstar+1];
    vval_Gstar = new double[numnz_Gstar];

    for(jj=0; jj<numnz_Gstar; jj++){
        vind_Gstar[jj] = row_vals[jj];
        vval_Gstar[jj] = A_sp_Gstar_t.values[jj];
    }

    for(jj=0; jj<k_Gstar+1; jj++){
        vbeg_Gstar[jj] = col_vals[jj];
    }
    //
    // Gbar

    arma::sp_mat A1 = arma::speye(nbY,nbY);
    arma::sp_mat A2(nbY,aux_nbDraws);
    arma::sp_mat A3 = kron_sp(arma::speye(nbY,nbY),arma::ones(aux_nbDraws,1));
    arma::sp_mat A4 = - kron_sp(arma::ones(nbY,1),arma::speye(aux_nbDraws,aux_nbDraws));
    arma::sp_mat A5(aux_nbDraws,nbY);
    arma::sp_mat A6 = -arma::speye(aux_nbDraws,aux_nbDraws);

    arma::sp_mat A_sp_Gbar = arma::join_cols(arma::join_rows(A1,A2), arma::join_cols(arma::join_rows(A3,A4),arma::join_rows(A5,A6)));

    k_Gbar = A_sp_Gbar.n_rows;
    n_Gbar = A_sp_Gbar.n_cols;

    arma::sp_mat A_sp_Gbar_t = arma::trans(A_sp_Gbar); // need to transpose to get data into CSR format (not CSC)

    numnz_Gbar = nonzeros(A_sp_Gbar).n_elem;

    const arma::uword* row_vals_2 = &(*A_sp_Gbar_t.row_indices);
    const arma::uword* col_vals_2 = &(*A_sp_Gbar_t.col_ptrs);

    vind_Gbar = new int[numnz_Gbar];
    vbeg_Gbar = new int[k_Gbar+1];
    vval_Gbar = new double[numnz_Gbar];

    for(jj=0; jj<numnz_Gbar; jj++){
        vind_Gbar[jj] = row_vals_2[jj];
        vval_Gbar[jj] = A_sp_Gbar_t.values[jj];
    }

    for(jj=0; jj<k_Gbar+1; jj++){
        vbeg_Gbar[jj] = col_vals_2[jj];
    }
}*/

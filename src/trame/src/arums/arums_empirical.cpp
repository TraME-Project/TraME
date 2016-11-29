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
 * empirical class
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 11/28/2016
 */

#include "trame.hpp"

trame::empirical::empirical(int nbX_inp, int nbY_inp, arma::cube atoms_inp, bool xHomogenous_inp, bool outsideOption_inp)
{
    this->build(nbX_inp, nbY_inp, atoms_inp, xHomogenous_inp, outsideOption_inp);
}

void trame::empirical::build(int nbX_inp, int nbY_inp, arma::cube atoms_inp, bool xHomogenous_inp, bool outsideOption_inp)
{
    nbX = nbX_inp;
    nbY = nbY_inp;

    atoms = atoms_inp;

    nbParams = atoms_inp.n_elem;
    aux_nbDraws = atoms.n_rows;

    xHomogenous = xHomogenous_inp;
    outsideOption = outsideOption_inp;
}

double trame::empirical::G(arma::vec n)
{
    double val = this->G(n,U,mu_sol);
    //
    return val;
}

double trame::empirical::G(arma::vec n, const arma::mat& U_inp, arma::mat& mu_out)
{
    int i;
    double val=0.0, val_x;

    mu_out.set_size(nbX,nbY);
    arma::mat mu_x_temp;
    //
    for (i=0; i<nbX; i++) {
        val_x = Gx(U_inp.row(i).t(), mu_x_temp, i);

        val += n(i)*val_x;
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double trame::empirical::Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, int x)
{
    arma::mat U_xs, Utilde;

    if (outsideOption) {
        U_xs = arma::join_cols(arma::vectorise(U_x_inp),arma::zeros(1,1));
    } else {
        U_xs = U_x_inp;
    }

    if (xHomogenous) {
        Utilde = arma::ones(aux_nbDraws,1) * U_xs.t() + atoms.slice(0);
    } else {
        Utilde = arma::ones(aux_nbDraws,1) * U_xs.t() + atoms.slice(x);
    }
    //
    int tt;
    arma::vec argmaxs = arma::max(Utilde,1);
    arma::uvec argmax_inds = which_max(Utilde, 1);

    double thesum = 0.0;
    for (tt=0; tt < aux_nbDraws; tt++) {
        thesum += argmaxs(tt,0);
    }

    double val_x = thesum/(double)(aux_nbDraws);
    //
    mu_x_out.set_size(nbY,1);
    arma::uvec temp_find;

    for (tt=0; tt<nbY; tt++) {
        temp_find = arma::find(argmax_inds==tt);
        mu_x_out(tt,0) = (double)(temp_find.n_elem)/(double)(aux_nbDraws);
    }
    //
    return val_x;
}

double trame::empirical::Gstar(arma::vec n)
{
    double val = this->Gstar(n,mu_sol,U_sol);
    //
    return val;
}

double trame::empirical::Gstar(arma::vec n, const arma::mat& mu_inp, arma::mat& U_out)
{
    if (!TRAME_PRESOLVED_GSTAR) {
        presolve_LP_Gstar();
    }
    //
    double val = 0.0, val_x = 0.0;

    U_out.set_size(nbX,nbY);
    arma::mat U_x;
    //
    for (int i=0; i<nbX; i++) {
        val_x = Gstarx((mu_inp.row(i).t())/n(i),U_x,i);

        val += n(i)*val_x;
        U_out.row(i) = arma::trans(U_x);
    }
    //
    return val;
}

double trame::empirical::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x)
{
    if (!TRAME_PRESOLVED_GSTAR) {
        presolve_LP_Gstar();
    }
    //
    int jj;
    double val_x = 0.0;
    arma::mat Phi, U_x_temp;

    if (xHomogenous) {
        Phi = atoms.slice(0);
    } else {
        Phi = atoms.slice(x);
    }
    //
    arma::vec p = arma::ones(aux_nbDraws,1)/aux_nbDraws;
    arma::mat q;

    if (outsideOption) {
        arma::mat temp_q(1,1);
        temp_q(0,0) = 1 - arma::accu(mu_x_inp);
        q = arma::join_cols(arma::vectorise(mu_x_inp),temp_q);
    } else {
        q = arma::vectorise(mu_x_inp);
    }
    //
    arma::vec obj_grbi = arma::vectorise(Phi);

    arma::vec rhs_grbi = arma::join_cols(p,q);

    char* sense_grbi = new char[k_Gstar];
    for (jj=0; jj<k_Gstar; jj++) {
        sense_grbi[jj] = '=';
    }

    bool LP_optimal;
    int modelSense = 1; // maximize
    double objval;

    arma::mat sol_mat(n_Gstar, 2);
    arma::mat dual_mat(k_Gstar, 2);
    try {
        //LP_optimal = generic_LP((int) A_grbi.n_rows, (int) A_grbi.n_cols, obj_grbi.memptr(), A_grbi.memptr(), modelSense, rhs_grbi.memptr(), sense_grbi, NULL, NULL, NULL, NULL, objval, sol_mat, dual_mat);
        //LP_optimal = generic_LP(k_Gstar, n_Gstar, obj_grbi.memptr(), numnz_Gstar, vbeg_Gstar, vind_Gstar, vval_Gstar, modelSense, rhs_grbi.memptr(), sense_grbi, NULL, NULL, NULL, NULL, objval, sol_mat, dual_mat);
        LP_optimal = generic_LP(k_Gstar, n_Gstar, obj_grbi.memptr(), numnz_Gstar, vbeg_Gstar, vind_Gstar, vval_Gstar, modelSense, rhs_grbi.memptr(), sense_grbi, NULL, NULL, NULL, NULL, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));
        //
        if (LP_optimal) {
            arma::mat u = dual_mat.col(0).rows(0,aux_nbDraws-1);
            if (outsideOption) {
                U_x_temp = dual_mat.col(0).rows(aux_nbDraws,aux_nbDraws+nbY);
                U_x_out = - U_x_temp.rows(0,nbY-1) + arma::as_scalar(U_x_temp.row(nbY));
            } else {
                U_x_temp = dual_mat.col(0).rows(aux_nbDraws,aux_nbDraws+nbY-1);
                U_x_out = -U_x_temp - arma::accu(p % u);
            }
            //
            val_x = -objval;
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
    delete[] sense_grbi;
    //
    return val_x;
}

double trame::empirical::Gbar(const arma::mat& Ubar, const arma::mat& mubar, const arma::vec& n, arma::mat& U_out, arma::mat& mu_out)
{
    if (!TRAME_PRESOLVED_GBAR) {
        presolve_LP_Gbar();
    }
    //
    double val=0.0, val_temp;

    U_out.set_size(nbX,nbY);
    mu_out.set_size(nbX,nbY);
    arma::mat U_x_temp, mu_x_temp;
    //
    for (int i=0; i<nbX; i++) {
        val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),U_x_temp,mu_x_temp,i);

        val += n(i)*val_temp;
        U_out.row(i) = arma::trans(U_x_temp);
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double trame::empirical::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, int x)
{
    if (!TRAME_PRESOLVED_GBAR) {
        presolve_LP_Gbar();
    }
    //
    int jj;
    double val_x=0.0;
    arma::mat Phi, U_x_temp;

    if (!outsideOption) {
        printf("Gbarx not implemented for empirical with outsideOption = false\n");
        return 0;
    }

    if (xHomogenous) {
        Phi = atoms.slice(0);
    } else {
        Phi = atoms.slice(x);
    }
    //
    arma::vec obj_grbi_1 = mubar_x;
    arma::vec obj_grbi_2 = - arma::ones(aux_nbDraws,1)/aux_nbDraws;
    arma::vec obj_grbi   = arma::join_cols(obj_grbi_1,obj_grbi_2);

    arma::vec rhs_grbi_1 = Ubar_x;
    arma::vec rhs_grbi_2 = arma::vectorise(-Phi);
    arma::vec rhs_grbi = arma::join_cols(rhs_grbi_1,rhs_grbi_2);

    char* sense_grbi = new char[k_Gbar];
    for (jj=0; jj<k_Gbar; jj++) {
        sense_grbi[jj] = '<';
    }

    bool LP_optimal;
    int modelSense = 1; // maximize
    double objval;

    arma::mat sol_mat(n_Gbar,2);
    arma::mat dual_mat(k_Gbar,2);

    try {
        //LP_optimal = generic_LP((int) A_grbi.n_rows, (int) A_grbi.n_cols, obj_grbi.memptr(), A_grbi.memptr(), modelSense, rhs_grbi.memptr(), sense_grbi, NULL, NULL, NULL, NULL, objval, sol_mat, dual_mat);
        //LP_optimal = generic_LP(k_Gbar, n_Gbar, obj_grbi.memptr(), numnz_Gbar, vbeg_Gbar, vind_Gbar, vval_Gbar, modelSense, rhs_grbi.memptr(), sense_grbi, NULL, NULL, NULL, NULL, objval, sol_mat, dual_mat);
        LP_optimal = generic_LP(k_Gbar, n_Gbar, obj_grbi.memptr(), numnz_Gbar, vbeg_Gbar, vind_Gbar, vval_Gbar, modelSense, rhs_grbi.memptr(), sense_grbi, NULL, NULL, NULL, NULL, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));
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
#if !defined(TRAME_USE_GUROBI_C)
    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
#endif
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //
    delete[] sense_grbi;
    //
    return val_x;
}

/*
 * presolve functions for Gstar and Gbar
 */

void trame::empirical::presolve_LP_Gstar()
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

    int jj, kk, count_val=0;

    arma::umat location_mat(2,aux_nbDraws*nbOptions*2);
    arma::rowvec vals_mat(aux_nbDraws*nbOptions*2);

    for (kk=0; kk<nbOptions; kk++) {
        for (jj=0; jj<aux_nbDraws; jj++) {
            location_mat(0,count_val) = jj + kk*aux_nbDraws;
            location_mat(1,count_val) = jj;
            ++count_val;
        }
        for (jj=0; jj<aux_nbDraws; jj++) {
            location_mat(0,count_val) = jj + kk*aux_nbDraws;
            location_mat(1,count_val) = kk + aux_nbDraws;
            ++count_val;
        }
    }

    vals_mat.fill(1);

    arma::sp_mat A_sp_Gstar_t(location_mat,vals_mat); // this is the transpose of A_sp_Gstar

    k_Gstar = A_sp_Gstar_t.n_cols; // cols as we're working with the transpose
    n_Gstar = A_sp_Gstar_t.n_rows; // rows as we're working with the transpose

    numnz_Gstar = aux_nbDraws*nbOptions*2;

    const arma::uword* row_vals = &(*A_sp_Gstar_t.row_indices);
    const arma::uword* col_vals = &(*A_sp_Gstar_t.col_ptrs);

    vind_Gstar = new int[numnz_Gstar];
    vbeg_Gstar = new int[k_Gstar+1];
    vval_Gstar = new double[numnz_Gstar];

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

void trame::empirical::presolve_LP_Gbar()
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

    int jj, kk, count_val=0;

    arma::umat location_mat_2(2,nbY + aux_nbDraws*nbY + aux_nbDraws*(nbY+1));
    arma::rowvec vals_mat_2(nbY + aux_nbDraws*nbY + aux_nbDraws*(nbY+1));

    for (jj=0; jj < nbY; jj++) { // top-left diagonal block
        location_mat_2(0,count_val) = jj;
        location_mat_2(1,count_val) = jj;

        vals_mat_2(count_val) = 1;

        ++count_val;
    }

    for (kk=0; kk < (nbY+1); kk++) {
        if (kk < nbY-1) { // top section 
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

    numnz_Gbar = nbY + aux_nbDraws*nbY*2 - aux_nbDraws;

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
    arma::sp_mat A_sp_Gstar = arma::join_cols(kron_sp(arma::ones(1,nbOptions),arma::speye(aux_nbDraws,aux_nbDraws)),kron_sp(arma::speye(nbOptions,nbOptions),arma::ones(1,aux_nbDraws)));

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

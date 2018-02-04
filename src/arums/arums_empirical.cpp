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
 * empirical additive random utility model (ARUM) class
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 02/04/2017
 */

#include "ancillary/ancillary.hpp"
#include "arums/arums.hpp"

//
// build functions

trame::arums::empirical::empirical(const uint_t nbX_inp, const uint_t nbY_inp)
{
    this->build(nbX_inp, nbY_inp);
}

trame::arums::empirical::empirical(const uint_t nbX_inp, const uint_t nbY_inp, const arma::cube& atoms_inp, const bool x_homogeneous_inp, const bool outside_option_inp)
{
    this->build(nbX_inp, nbY_inp, atoms_inp, x_homogeneous_inp, outside_option_inp);
}

void
trame::arums::empirical::build(const uint_t nbX_inp, const uint_t nbY_inp)
{
    nbX = nbX_inp;
    nbY = nbY_inp;
}

void
trame::arums::empirical::build(const uint_t nbX_inp, const uint_t nbY_inp, const arma::cube& atoms_inp, const bool x_homogeneous_inp, const bool outside_option_inp)
{
    nbX = nbX_inp;
    nbY = nbY_inp;

    atoms = atoms_inp;

    dim_params = atoms_inp.n_elem;
    aux_n_draws = atoms.n_rows;

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

    for (uint_t i=0; i < nbX; i++) 
    {
        val_x = Gx(U_inp.row(i).t(), mu_x_temp, i);

        val += n(i)*val_x;
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }

    //

    return val;
}

double
trame::arums::empirical::Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, const uint_t x)
const
{
    mu_x_out.set_size(nbY,1);
    
    arma::mat U_xs = (outside_option) ? arma::join_cols(arma::vectorise(U_x_inp),arma::zeros(1,1)) : U_x_inp;
    const arma::mat Utilde = (x_homogeneous) ? arma::ones(aux_n_draws,1) * U_xs.t() + atoms.slice(0) : arma::ones(aux_n_draws,1) * U_xs.t() + atoms.slice(x);

    //
    
    const arma::vec argmaxs = arma::max(Utilde,1);       // take max over dim = 1
    const arma::uvec argmax_inds = which_max(Utilde,1);

    double val_x = arma::accu(argmaxs) / static_cast<double>(aux_n_draws);

    //

    arma::uvec temp_find;

    for (uint_t tt=0; tt < nbY; tt++)
    {
        temp_find = arma::find(argmax_inds==tt);
        mu_x_out(tt,0) = static_cast<double>(temp_find.n_elem)/static_cast<double>(aux_n_draws);
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
        presolve_LP_Gstar();
    }

    double val = 0.0, val_x = 0.0;
    U_out.set_size(nbX,nbY);

    //

    arma::mat U_x;

    for (uint_t i=0; i < nbX; i++)
    {
        val_x = Gstarx((mu_inp.row(i).t())/n(i),U_x,i);

        val += n(i)*val_x;
        U_out.row(i) = arma::trans(U_x);
    }

    //

    return val;
}

double
trame::arums::empirical::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, const uint_t x)
const
{
    if (!TRAME_PRESOLVED_GSTAR) {
        presolve_LP_Gstar();
    }

    //

    const arma::mat Phi = (x_homogeneous) ? atoms.slice(0) : atoms.slice(x);
    
    arma::vec p = arma::ones(aux_n_draws,1)/ static_cast<double>(aux_n_draws);
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

    //

    char* sense_lp = new char[k_Gstar];
    std::memset(sense_lp, '=', k_Gstar * sizeof (char));

    //

    bool lp_optimal;
    int modelSense = 1; // maximize
    double objval, val_x = 0.0;

    arma::mat sol_mat(n_Gstar, 2);
    arma::mat dual_mat(k_Gstar, 2);
    
    try {
        lp_optimal = generic_LP(k_Gstar, n_Gstar, obj_lp.memptr(), num_non_zero_Gstar, vbeg_Gstar, vind_Gstar, vval_Gstar, modelSense, rhs_lp.memptr(), sense_lp, nullptr, nullptr, nullptr, nullptr, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));
        
        //

        if (lp_optimal)
        {
            arma::mat u = dual_mat.col(0).rows(0,aux_n_draws-1);

            if (outside_option) {
                const arma::mat U_x_temp = dual_mat.col(0).rows(aux_n_draws,aux_n_draws+nbY);
                U_x_out = - U_x_temp.rows(0,nbY-1) + arma::as_scalar(U_x_temp.row(nbY));
            } else {
                const arma::mat U_x_temp = dual_mat.col(0).rows(aux_n_draws,aux_n_draws+nbY-1);
                U_x_out = - U_x_temp - arma::accu(p % u);
            }
            
            val_x = -objval;
        } else {
            std::cout << "Non-optimal value found during optimization" << std::endl;
        }
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }

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
        presolve_LP_Gbar();
    }

    //

    double val = 0.0;

    U_out.set_size(nbX,nbY);
    mu_out.set_size(nbX,nbY);
    arma::mat U_x_temp, mu_x_temp;

    //
    
    for (uint_t i=0; i < nbX; i++)
    {
        const double val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),U_x_temp,mu_x_temp,i);

        val += n(i)*val_temp;
        U_out.row(i) = arma::trans(U_x_temp);
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }

    //

    return val;
}

double
trame::arums::empirical::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, const uint_t x)
const
{
    if (!TRAME_PRESOLVED_GBAR) {
        presolve_LP_Gbar();
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
    arma::vec obj_lp_2 = - arma::ones(aux_n_draws,1)/aux_n_draws;
    arma::vec obj_lp   = arma::join_cols(obj_lp_1,obj_lp_2);

    arma::vec rhs_lp_1 = Ubar_x;
    arma::vec rhs_lp_2 = arma::vectorise(-Phi);
    arma::vec rhs_lp = arma::join_cols(rhs_lp_1,rhs_lp_2);

    //

    char* sense_lp = new char[k_Gbar];
    std::memset(sense_lp, '<', k_Gbar * sizeof (char));

    //

    bool lp_optimal;
    int modelSense = 1; // maximize
    double objval, val_x = 0.0;;

    arma::mat sol_mat(n_Gbar,2);
    arma::mat dual_mat(k_Gbar,2);

    try {
        lp_optimal = generic_LP(k_Gbar, n_Gbar, obj_lp.memptr(), num_non_zero_Gbar, vbeg_Gbar, vind_Gbar, vval_Gbar, modelSense, rhs_lp.memptr(), sense_lp, nullptr, nullptr, nullptr, nullptr, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));
        //
        if (lp_optimal)
        {
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
// D2G

arma::mat
trame::arums::empirical::D2G(const arma::vec& n, const bool x_first)
const
{   
    arma::mat ret;

    printf("D2G not yet defined for the empirical case.\n");

    return ret;
}

void
trame::arums::empirical::D2G(arma::mat &H, const arma::vec& n, const bool x_first)
const
{   
    printf("D2G not yet defined for the empirical case.\n");
}

arma::mat
trame::arums::empirical::D2G(const arma::vec& n, const arma::mat& U_inp, const bool x_first)
const
{   
    arma::mat ret;

    printf("D2G not yet defined for the empirical case.\n");

    return ret;
}

void
trame::arums::empirical::D2G(arma::mat &H, const arma::vec& n, const arma::mat& U_inp, const bool x_first)
const
{   
    printf("D2G not yet defined for the empirical case.\n");
}

//
// D2Gstar

arma::mat
trame::arums::empirical::D2Gstar(const arma::vec& n, const bool x_first)
const
{   
    arma::mat ret;

    printf("D2Gstar not yet defined for the empirical case.\n");

    return ret;
}

void
trame::arums::empirical::D2Gstar(arma::mat &H, const arma::vec& n, const bool x_first)
const
{   
    printf("D2Gstar not yet defined for the empirical case.\n");
}

arma::mat
trame::arums::empirical::D2Gstar(const arma::vec& n, const arma::mat& mu_inp, const bool x_first)
const
{   
    arma::mat ret;

    printf("D2Gstar not yet defined for the empirical case.\n");

    return ret;
}

void
trame::arums::empirical::D2Gstar(arma::mat &H, const arma::vec& n, const arma::mat& mu_inp, const bool x_first)
const
{   
    printf("D2Gstar not yet defined for the empirical case.\n");
}

//
// dparams gradient

arma::mat
trame::arums::empirical::dparams_NablaGstar(const arma::vec& n, const arma::mat* dparams_inp, const bool x_first)
const
{
    arma::mat ret;
    
    printf("dparams_NablaGstar not yet defined for the empirical case.\n");

    return ret;
}

void
trame::arums::empirical::dparams_NablaGstar(arma::mat &ret, const arma::vec& n, const arma::mat* dparams_inp, const bool x_first)
const
{
    printf("dparams_NablaGstar not yet defined for the empirical case.\n");
}

arma::mat
trame::arums::empirical::dparams_NablaGstar(const arma::vec& n, const arma::mat& mu_inp, const arma::mat* dparams_inp, const bool x_first)
const
{
    arma::mat ret;
    
    printf("dparams_NablaGstar not yet defined for the empirical case.\n");

    return ret;
}

void
trame::arums::empirical::dparams_NablaGstar(arma::mat &ret, const arma::vec& n, const arma::mat& mu_inp, const arma::mat* dparams_inp, const bool x_first)
const
{
    printf("dparams_NablaGstar not yet defined for the empirical case.\n");
}

//
// presolve functions for Gstar and Gbar

void
trame::arums::empirical::presolve_LP_Gstar()
const
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

    num_non_zero_Gstar = aux_n_draws*nb_options*2;

    arma::umat location_mat(2,num_non_zero_Gstar);
    arma::rowvec vals_mat = arma::ones(1,num_non_zero_Gstar);

    uint_t count_val = 0;

    for (uint_t kk=0; kk < nb_options; kk++)
    {
        for (uint_t jj=0; jj < aux_n_draws; jj++)
        {
            location_mat(0,count_val) = jj + kk*aux_n_draws;
            location_mat(1,count_val) = jj;
            ++count_val;
        }
        for (uint_t jj=0; jj < aux_n_draws; jj++)
        {
            location_mat(0,count_val) = jj + kk*aux_n_draws;
            location_mat(1,count_val) = kk + aux_n_draws;
            ++count_val;
        }
    }

    //

    arma::sp_mat A_lp_t(location_mat,vals_mat); // this is the transpose of the constraint matrix

    k_Gstar = A_lp_t.n_cols; // n_cols as we are working with the transpose of A
    n_Gstar = A_lp_t.n_rows; // n_rows as we are working with the transpose of A

    vind_Gstar = uword_to_int(A_lp_t.row_indices,num_non_zero_Gstar); // index of what row each non-zero value belongs to
    vbeg_Gstar = uword_to_int(A_lp_t.col_ptrs,k_Gstar+1);    // index of how many non-zero values are in each column

    vval_Gstar = new double[num_non_zero_Gstar];
    std::memcpy(vval_Gstar, A_lp_t.values, num_non_zero_Gstar * sizeof(double));

    //

    TRAME_PRESOLVED_GSTAR = true;
}

void
trame::arums::empirical::presolve_LP_Gbar()
const
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

    num_non_zero_Gbar = nbY + aux_n_draws*nbY + aux_n_draws*(nbY+1);

    arma::umat location_mat_2(2,num_non_zero_Gbar);
    arma::rowvec vals_mat_2(num_non_zero_Gbar);

    int count_val = 0;

    for (uint_t jj=0; jj < nbY; jj++)
    {   // top-left diagonal block
        location_mat_2(0,count_val) = jj;
        location_mat_2(1,count_val) = jj;

        vals_mat_2(count_val) = 1;

        ++count_val;
    }

    for (uint_t kk=0; kk < (nbY+1); kk++)
    {
        if (kk < nbY) 
        {   // top section 
            for (uint_t jj=0; jj < aux_n_draws; jj++)
            {
                location_mat_2(0,count_val) = kk;
                location_mat_2(1,count_val) = nbY + jj + kk*aux_n_draws;

                vals_mat_2(count_val) = 1;

                ++count_val;
            }
        }

        for (uint_t jj=0; jj < aux_n_draws; jj++)
        {   // diagonal terms (nbY+1 number of blocks)
            location_mat_2(0,count_val) = nbY + jj;
            location_mat_2(1,count_val) = nbY + jj + kk*aux_n_draws;

            vals_mat_2(count_val) = -1;

            ++count_val;
        }
    }

    //

    arma::sp_mat A_lp_t(location_mat_2,vals_mat_2); // this is the transpose of the constraint matrix

    k_Gbar = A_lp_t.n_cols; // n_cols as we are working with the transpose of A
    n_Gbar = A_lp_t.n_rows; // n_rows as we are working with the transpose of A

    vind_Gbar = uword_to_int(A_lp_t.row_indices,num_non_zero_Gbar); // index of what row each non-zero value belongs to
    vbeg_Gbar = uword_to_int(A_lp_t.col_ptrs,k_Gbar+1);    // index of how many non-zero values are in each column

    vval_Gbar = new double[num_non_zero_Gbar];
    std::memcpy(vval_Gbar, A_lp_t.values, num_non_zero_Gbar * sizeof(double));

    //

    TRAME_PRESOLVED_GBAR = true;
}

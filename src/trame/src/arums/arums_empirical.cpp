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
 */

#include "trame.hpp"

void empirical::build(int nbX_inp, int nbY_inp, arma::cube atoms_inp, bool xHomogenous_inp, bool outsideOption_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
    
    atoms = atoms_inp;
    
    nbParams = atoms_inp.n_elem;
    aux_nbDraws = atoms.n_rows;
    
    xHomogenous = xHomogenous_inp;
    outsideOption = outsideOption_inp;
}

void empirical::presolve_LP_Gstar()
{   
    /*
     * Here we build and store the 'A' matrix that get passed to 
     * a linear programming solver in Gstarx. For this
     * we use batch allocation of the elements of sparse matrices. 
     * This is *much* faster than first allocating the sparse 
     * matrix A_sp_* then consecutively adding elements.
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

void empirical::presolve_LP_Gbar()
{   
    /*
     * Here we build and store the 'A' matrix that get passed to 
     * a linear programming solver in Gbarx. For this
     * we use batch allocation of the elements of sparse matrices. 
     * This is *much* faster than first allocating the sparse 
     * matrix A_sp_* then consecutively adding elements.
     */
    
    int jj, kk, count_val=0;
     
    arma::umat location_mat_2(2,nbY + aux_nbDraws*nbOptions*2 - aux_nbDraws);
    arma::rowvec vals_mat_2(nbY + aux_nbDraws*nbOptions*2 - aux_nbDraws);
    
    for (jj=0; jj<nbY; jj++) {
        location_mat_2(0,count_val) = jj;
        location_mat_2(1,count_val) = jj;
        
        vals_mat_2(count_val) = 1;
        
        ++count_val;
    }
    
    for (kk=0; kk<nbOptions; kk++) {
        if (kk < nbOptions-1) {
            for (jj=0; jj<aux_nbDraws; jj++) {
                location_mat_2(0,count_val) = kk;
                location_mat_2(1,count_val) = nbY + jj + kk*aux_nbDraws;
                
                vals_mat_2(count_val) = 1;
                
                ++count_val;
            }
        }
        
        for (jj=0; jj<aux_nbDraws; jj++) { // diagonal terms
            location_mat_2(0,count_val) = nbY + jj;
            location_mat_2(1,count_val) = nbY + jj + kk*aux_nbDraws;
            
            vals_mat_2(count_val) = -1;
            
            ++count_val;
        }
    }
    
    arma::sp_mat A_sp_Gbar_t(location_mat_2,vals_mat_2);
    
    k_Gbar = A_sp_Gbar_t.n_cols; // cols as we're working with the transpose
    n_Gbar = A_sp_Gbar_t.n_rows; // rows as we're working with the transpose
    
    numnz_Gbar = nbY + aux_nbDraws*nbOptions*2 - aux_nbDraws;
    
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

double empirical::G(arma::vec n)
{   
    int i;
    double val=0.0, val_x;
    
    mu_sol.set_size(nbX,nbY);
    arma::mat mu_x_temp;
    //
    for (i=0; i<nbX; i++) {
        val_x = Gx(U.row(i).t(),mu_x_temp,i);
        
        val += n(i)*val_x;
        mu_sol.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double empirical::Gx(arma::mat Ux, arma::mat& mu_x_inp, int x)
{   
    arma::mat Uxs, Utilde;
    
    if (outsideOption) {
        Uxs = arma::join_cols(arma::vectorise(Ux),arma::zeros(1,1));
    } else {
        Uxs = Ux;
    }
    
    if (xHomogenous) {
        Utilde = arma::ones(aux_nbDraws,1) * Uxs.t() + atoms.slice(0);
    } else {
        Utilde = arma::ones(aux_nbDraws,1) * Uxs.t() + atoms.slice(x);
    }
    //
    int tt;
    arma::vec argmaxs = arma::max(Utilde,1);
    arma::uvec argmax_inds = which_max(&Utilde, 1);
    
    double thesum = 0.0;
    for (tt=0; tt < aux_nbDraws; tt++) {
        thesum += argmaxs(tt,0);
    }

    double valx = thesum/(double)(aux_nbDraws);
    //
    mu_x_inp.set_size(nbY,1);
    arma::uvec temp_find;
    
    for (tt=0; tt<nbY; tt++) {
        temp_find = arma::find(argmax_inds==tt);
        mu_x_inp(tt,0) = (double)(temp_find.n_elem)/(double)(aux_nbDraws);
    }
    //
    return valx;
}

double empirical::Gstar(arma::vec n)
{   
    int i;
    double val=0.0, val_temp;
    
    if (TRAME_PRESOLVED_GSTAR!=true) {
        presolve_LP_Gstar();
    }
    
    U_sol.set_size(nbX,nbY);
    arma::mat Ux_temp;
    //
    for (i=0; i<nbX; i++) {
        val_temp = Gstarx((mu_sol.row(i).t())/n(i),Ux_temp,i);
        
        val += n(i)*val_temp;
        U_sol.row(i) = arma::trans(Ux_temp);
    }
    //
    return val;
}

double empirical::Gstar(arma::mat& U_inp, arma::vec n)
{   
    int i;
    double val=0.0, val_temp;
    
    if (TRAME_PRESOLVED_GSTAR!=true) {
        presolve_LP_Gstar();
    }
    
    U_inp.set_size(nbX,nbY);
    arma::mat Ux_temp;
    //
    for (i=0; i<nbX; i++) {
        val_temp = Gstarx((mu_sol.row(i).t())/n(i),Ux_temp,i);
        
        val += n(i)*val_temp;
        U_inp.row(i) = arma::trans(Ux_temp);
    }
    //
    return val;
}

double empirical::Gstarx(arma::mat mu_x, arma::mat& Ux_inp, int x)
{   
    int jj;
    double valx=0.0;
    arma::mat Phi, Ux_temp;
    
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
        temp_q(0,0) = 1 - arma::accu(mu_x);
        q = arma::join_cols(arma::vectorise(mu_x),temp_q);
    } else {
        q = arma::vectorise(mu_x);
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
                Ux_temp = dual_mat.col(0).rows(aux_nbDraws,aux_nbDraws+nbY);
                Ux_inp = - Ux_temp.rows(0,nbY-1) + arma::as_scalar(Ux_temp.row(nbY));
            } else { 
                Ux_temp = dual_mat.col(0).rows(aux_nbDraws,aux_nbDraws+nbY-1);
                Ux_inp = -Ux_temp - arma::accu(p % u);
            }
            //
            valx = -objval;
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
    return valx;
}

double empirical::Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_inp, arma::mat& mu_inp)
{   
    int i;
    double val=0.0, val_temp;
    
    if (!TRAME_PRESOLVED_GBAR) {
        presolve_LP_Gbar();
    }
    
    U_inp.set_size(nbX,nbY);
    mu_inp.set_size(nbX,nbY);
    arma::mat Ux_temp, mu_x_temp;
    //
    for (i=0; i<nbX; i++) {
        val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),Ux_temp,mu_x_temp,i);
        
        val += n(i)*val_temp;
        U_inp.row(i) = arma::trans(Ux_temp);
        mu_inp.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double empirical::Gbarx(arma::vec Ubarx, arma::vec mubarx, arma::mat& Ux_inp, arma::mat& mu_x_inp, int x)
{   
    int jj;
    double valx=0.0;
    arma::mat Phi, Ux_temp;
    
    if (!outsideOption) {
        printf("Gbarx not implemented for empirical with outsideOption = false\n");
        return 0; //Gbarx not implemented for empirical with outsideOption = false
    }
    
    if (xHomogenous) {
        Phi = atoms.slice(0);
    } else {
        Phi = atoms.slice(x);
    }
    //
    arma::vec obj_grbi_1 = mubarx;
    arma::vec obj_grbi_2 = - arma::ones(aux_nbDraws,1)/aux_nbDraws;
    arma::vec obj_grbi   = arma::join_cols(obj_grbi_1,obj_grbi_2);
 
    arma::vec rhs_grbi_1 = Ubarx;
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
            Ux_inp = sol_mat.col(0).rows(0,nbY-1);
            arma::vec delta_mu_x = dual_mat.col(0).rows(0,nbY-1);
            mu_x_inp = mubarx - delta_mu_x;
            //
            valx = arma::accu(mubarx % Ubarx) - objval;
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
    return valx;
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

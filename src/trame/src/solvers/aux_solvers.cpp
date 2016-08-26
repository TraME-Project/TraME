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
 * auxiliary functions for equilibrium solvers
 *
 * Keith O'Hara
 * 08/23/2016
 */

#include "trame.hpp"

arma::mat trame::u_from_vs(trame::transfers trans_obj, arma::mat v, double* tol_inp, arma::mat* subdiff)
{
    arma::mat us = trans_obj.Ucal(v,NULL,NULL);
    arma::mat u  = arma::max(arma::max(us,1),arma::zeros(us.n_rows,1));
    //
    double tol;
    if (tol_inp) {
        tol = *tol_inp;
    } else {
        tol = 0;
    }
    
    if (subdiff) {
        *subdiff = arma::zeros(trans_obj.nbX,trans_obj.nbY);
        subdiff->elem( arma::find(arma::abs(elem_sub(u,us)) <= tol) ).ones();
    }
    //
    return u;
}

arma::mat trame::v_from_us(trame::transfers trans_obj, arma::mat u, double* tol_inp, arma::mat* subdiff)
{
    arma::mat vs = trans_obj.Vcal(u,NULL,NULL);
    arma::mat v  = arma::trans(arma::max(arma::max(vs,0),arma::zeros(1,vs.n_cols)));
    //
    double tol;
    if (tol_inp) {
        tol = *tol_inp;
    } else {
        tol = 0;
    }

    if (subdiff) {
        *subdiff = arma::zeros(trans_obj.nbY,trans_obj.nbX);
        subdiff->elem( arma::find(arma::abs(elem_sub(v,vs.t())) <= tol) ).ones();
        *subdiff = subdiff->t();
    }
    //
    return v;
}

arma::mat trame::update_v(trame::transfers trans_obj, arma::mat v, arma::vec n, arma::vec m, bool xFirst)
{
    int nbX = trans_obj.nbX;
    int nbY = trans_obj.nbY;
    
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

    int k_lp = A_lp.n_rows;
    int n_lp = A_lp.n_cols;

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

    int k_bis = A_bis.n_rows;
    int n_bis = A_bis.n_cols;

    char* sense_bis = new char[k_bis];
    for (jj=0; jj<k_bis-1; jj++) {
        sense_bis[jj] = '>';
    }
    sense_bis[k_bis-1] = '=';

    int modelSense_bis = 0; // minimize
    if (!xFirst) {
        modelSense_bis = 1; // maximize
    }

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
#if !defined(TRAME_USE_GUROBI_C)
        } catch(GRBException e) {
            std::cout << "Error code = " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
#endif
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }
        /*
         * second LP
         */
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
#if !defined(TRAME_USE_GUROBI_C)
        } catch(GRBException e) {
            std::cout << "Error code = " << e.getErrorCode() << std::endl;
            std::cout << e.getMessage() << std::endl;
#endif
        } catch(...) {
            std::cout << "Exception during optimization" << std::endl;
        }
    }
    //
    return v_updated;
}

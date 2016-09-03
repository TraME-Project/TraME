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
 * Generic Linear/Quadratic Programming solver using Gurobi
 *
 * Keith O'Hara
 * 05/08/2016
 */

#include "trame.hpp"

#if !defined(TRAME_PREDEF_GUROBI_C)

    #include "gurobi_c++.h"

    bool trame::generic_LP(int k, int n, double *obj, double* A, int model_opt_sense, double* rhs, char* constr_sense, double* Q, double* lb, double* ub, double* start, double& objval, arma::mat& sol_mat, arma::mat& dual_mat)
    {
        // k: number of constraints ('rows')
        // n: number of variables ('columns')
        //
        // Initialize
        bool success = false;
        
        int i,j;
        //
        // environment
        GRBEnv* env = 0;
        env = new GRBEnv();
        
        GRBModel model = GRBModel(*env);
        model.getEnv().set(GRB_IntParam_OutputFlag, 0);
        //
        // Model variables and constraints
        GRBVar* vars = model.addVars(lb, ub, NULL, NULL, NULL, n);
        model.update();
        
        for (i=0; i < k; i++) {
            GRBLinExpr lhs = 0;
            for (j=0; j < n; j++) {
                if (A[i+j*k] != 0) {
                    lhs += A[i+j*k]*vars[j]; // need to switch order; was A[i*n+j]
                }
            }
            model.addConstr(lhs, constr_sense[i], rhs[i]);
        }
        //
        // Model objective
        if (Q==NULL) { // Linear programming problem
            GRBLinExpr LPobj = 0;
            
            for (j=0; j < n; j++) {
                LPobj += obj[j]*vars[j];
            }
            
            if (model_opt_sense==1) {
                model.setObjective(LPobj,GRB_MAXIMIZE);
            } else {
                model.setObjective(LPobj,GRB_MINIMIZE);
            }
        } else { // Quadratic programming problem
            GRBQuadExpr QPobj = 0;
            
            for (j=0; j < n; j++) {
                QPobj += obj[j]*vars[j];
            }
            for (i=0; i < n; i++) {
                for (j = 0; j < n; j++) {
                    if (Q[i+j*n] != 0) {
                        QPobj += Q[i+j*n]*vars[i]*vars[j];
                    }
                }
            }
            
            if (model_opt_sense==1) {
                model.setObjective(QPobj,GRB_MAXIMIZE);
            } else {
                model.setObjective(QPobj,GRB_MINIMIZE);
            }
        }
        
        model.update();
        //
        // Optimize and recover relevant solution objects
        model.optimize();
        
        GRBConstr* mycons = model.getConstrs();
        
        if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
            objval = model.get(GRB_DoubleAttr_ObjVal);
            for (i=0; i<n; i++) {
                sol_mat(i,0) = vars[i].get(GRB_DoubleAttr_X);
                sol_mat(i,1) = vars[i].get(GRB_DoubleAttr_RC);
            }
            for (j=0; j<k; j++) {
                dual_mat(j,0) = mycons[j].get(GRB_DoubleAttr_Pi);
                dual_mat(j,1) = mycons[j].get(GRB_DoubleAttr_Slack);
            }
            success = true;
        }
        //
        return success;
    }

#else

    extern "C" {
        #include "lp/generic_lp_c.h"
    }
    
    bool trame::generic_LP(int k, int n, double *obj, double* A, int model_opt_sense, double* rhs, char* constr_sense, double* Q, double* lb, double* ub, double* start, double& objval, double* sol_mat_X, double* sol_mat_RC, double* dual_mat_PI, double* dual_mat_SLACK)
    {
        // k: number of constraints ('rows')
        // n: number of variables ('columns')
        //
        // Initialize
        bool success = false;
        int solved;
        //
        // Call C-version of the solver
        solved = generic_LP_C_switch(k, n, obj, A, model_opt_sense, rhs, constr_sense, Q, lb, ub, 
                                     &objval, sol_mat_X, sol_mat_RC, dual_mat_PI, dual_mat_SLACK);
        //
        // Put solution matrices together
        if (solved == 1) {
            success = true;
        }
        //
        return success;
    }
    
    bool trame::generic_LP(int k, int n, double *obj, int numnz, int* vbeg, int* vind, double* vval, int model_opt_sense, double* rhs, char* constr_sense, double* Q, double* lb, double* ub, double* start, double& objval, double* sol_mat_X, double* sol_mat_RC, double* dual_mat_PI, double* dual_mat_SLACK)
    {
        // k: number of constraints ('rows')
        // n: number of variables ('columns')
        //
        // Initialize
        bool success = false;
        int solved;
        //
        // Call C-version of the solver
        solved = generic_LP_C_sparse(k, n, obj, numnz, vbeg, vind, vval, model_opt_sense, rhs, constr_sense, Q, lb, ub, 
                                     &objval, sol_mat_X, sol_mat_RC, dual_mat_PI, dual_mat_SLACK);
        //
        // Put solution matrices together
        if (solved == 1) {
            success = true;
        }
        //
        return success;
    }

#endif

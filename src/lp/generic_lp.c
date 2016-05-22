/*
 * Generic Linear/Quadratic Programming solver using Gurobi
 * Based on Gurobi's 'dense_c.c' example
 *
 * Keith O'Hara
 * 05/08/2016
 *
 * clang -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include generic_lp.c -c -o generic_lp.o
 */

#include "generic_lp.h"

int generic_LP_C(int rows, int cols, double* obj, double* A, int modelSense, double* rhs, char* sense, double* Q, double* lb, double* ub, double* objval, double* sol_mat_X, double* sol_mat_RC, double* dual_mat_PI, double* dual_mat_SLACK)
{
    int i, j, optimstatus;
    int error = 0;
    int success = 0;
    
    /* Create an empty model */
    
    GRBenv* env = NULL;
    
    error = GRBloadenv(&env, NULL);
    if (error) goto QUIT;
    
    GRBmodel *model = NULL;
    
    error = GRBnewmodel(env, &model, "trame_LP", cols, obj, lb, ub, NULL, NULL);
    if (error) goto QUIT;
    
    error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_OUTPUTFLAG, 0);
    if (error) goto QUIT;
    
    if(modelSense==0){ // minimize or maximize, resp.
        error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, 1);
    }else{
        error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, -1);
    }
    if (error) goto QUIT;

    error = GRBaddconstrs(model, rows, 0, NULL, NULL, NULL, sense, rhs, NULL);
    if (error) goto QUIT;
    
    /* Integrate new rows and columns */
    
    error = GRBupdatemodel(model);
    if (error) goto QUIT;
    
    /* Populate A matrix */
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            /* IMPORTANT: when using the generic_lp.hpp version, we switch
             * from:
             *
             * error = GRBchgcoeffs(model, 1, &i, &j, &A[i*rows+j]);
             *
             * to &A[i+j*cols] because of a weird issue where
             * the reference would run over rows rather than columns first
             */
            if (A[i*cols+j] != 0) {
                error = GRBchgcoeffs(model, 1, &i, &j, &A[i*cols+j]);
                if (error) goto QUIT;
            }
        }
    }
    
    /* Populate Q matrix */

    if (Q) {
        for (i = 0; i < cols; i++) {
            for (j = 0; j < cols; j++) {
                if (Q[i*cols+j] != 0) {
                    error = GRBaddqpterms(model, 1, &i, &j, &Q[i*cols+j]);
                    if (error) goto QUIT;
                }
            }
        }
    }
    
    error = GRBupdatemodel(model);
    if (error) goto QUIT;

    /* Optimize model */

    error = GRBoptimize(model);
    if (error) goto QUIT;

    /* Capture solution information */

    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
    if (error) goto QUIT;

    if (optimstatus == GRB_OPTIMAL) {

        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, objval);
        if (error) goto QUIT;

        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, cols, sol_mat_X);
        if (error) goto QUIT;
        
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_RC, 0, cols, sol_mat_RC);
        if (error) goto QUIT;
        
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, rows, dual_mat_PI);
        if (error) goto QUIT;
        
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_SLACK, 0, rows, dual_mat_SLACK);
        if (error) goto QUIT;

        success = 1;
    }

QUIT:

    /* Error reporting */

    if (error) {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }

    /* Free model */

    GRBfreemodel(model);
    GRBfreeenv(env);

    return success;
}

int generic_LP_C_switch(int rows, int cols, double* obj, double* A, int modelSense, double* rhs, char* sense, double* Q, double* lb, double* ub, double* objval, double* sol_mat_X, double* sol_mat_RC, double* dual_mat_PI, double* dual_mat_SLACK)
{
    int i, j, optimstatus;
    int error = 0;
    int success = 0;
    
    /* Create an empty model */
    
    GRBenv* env = NULL;
    
    error = GRBloadenv(&env, NULL);
    if (error) goto QUIT;
    
    GRBmodel *model = NULL;
    
    error = GRBnewmodel(env, &model, "trame_LP", cols, obj, lb, ub, NULL, NULL);
    if (error) goto QUIT;
    
    error = GRBsetintparam(GRBgetenv(model), GRB_INT_PAR_OUTPUTFLAG, 0);
    if (error) goto QUIT;
    
    if(modelSense==0){ // minimize or maximize, resp.
        error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, 1);
    }else{
        error = GRBsetintattr(model, GRB_INT_ATTR_MODELSENSE, -1);
    }
    if (error) goto QUIT;

    error = GRBaddconstrs(model, rows, 0, NULL, NULL, NULL, sense, rhs, NULL);
    if (error) goto QUIT;
    
    /* Integrate new rows and columns */
    
    error = GRBupdatemodel(model);
    if (error) goto QUIT;
    
    /* Populate A matrix */
    
    int act_row = 0;
    int act_col = 0;
    
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            /* IMPORTANT: when using the generic_lp.hpp version, we switch
             * from:
             *
             * error = GRBchgcoeffs(model, 1, &i, &j, &A[i*rows+j]);
             *
             * to &A[i+j*cols] because of a weird issue where
             * the reference would run over rows rather than columns first
             */
            if (A[i*cols+j] != 0) {
                error = GRBchgcoeffs(model, 1, &act_row, &act_col, &A[i*cols+j]);
                if (error) goto QUIT;
            }
            ++act_row;
            if(act_row>=rows){
                act_row=0;
                ++act_col;
            }
        }
    }
    
    /* Populate Q matrix */

    if (Q) {
        for (i = 0; i < cols; i++) {
            for (j = 0; j < cols; j++) {
                if (Q[i*cols+j] != 0) {
                    error = GRBaddqpterms(model, 1, &i, &j, &Q[i*cols+j]);
                    if (error) goto QUIT;
                }
            }
        }
    }
    
    error = GRBupdatemodel(model);
    if (error) goto QUIT;

    /* Optimize model */

    error = GRBoptimize(model);
    if (error) goto QUIT;

    /* Capture solution information */

    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, &optimstatus);
    if (error) goto QUIT;

    if (optimstatus == GRB_OPTIMAL) {

        error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, objval);
        if (error) goto QUIT;

        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, cols, sol_mat_X);
        if (error) goto QUIT;
        
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_RC, 0, cols, sol_mat_RC);
        if (error) goto QUIT;
        
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_PI, 0, rows, dual_mat_PI);
        if (error) goto QUIT;
        
        error = GRBgetdblattrarray(model, GRB_DBL_ATTR_SLACK, 0, rows, dual_mat_SLACK);
        if (error) goto QUIT;

        success = 1;
    }

QUIT:

    /* Error reporting */

    if (error) {
        printf("ERROR: %s\n", GRBgeterrormsg(env));
        exit(1);
    }

    /* Free model */

    GRBfreemodel(model);
    GRBfreeenv(env);

    return success;
}
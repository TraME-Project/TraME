/*
 * Generic Linear Programming solver using Gurobi
 */

#include "gurobi_c++.h"

void generic_LP(arma::vec *obj, arma::mat* A, arma::vec* rhs, char* sense, arma::mat* Q, arma::vec* lb, arma::vec* ub, arma::vec* start, double& objval, double& sol_vars, double& sol_PI, double& sol_RC)
{
    int i,j;
    int k = A->n_rows;   // number of constraints
    int n = obj->n_elem; // number of variables
    // environment
    GRBEnv* env = 0;
    env = new GRBEnv();
    
    GRBModel model = GRBModel(*env);
    model.getEnv().set(GRB_IntParam_OutputFlag, 0);
    //
    // Set bounds
    double** lb_grbi = new double*[n];
    double** ub_grbi = new double*[n];
    
    if(lb==NULL){
        *lb_grbi = NULL;
    }else{
        for(i=0; i<n; i++){
            lb_grbi[i] = &(*lb)[i];
        }
    }
    
    if(ub==NULL){
        *ub_grbi = NULL;
    }else{
        for(i=0; i<n; i++){
            ub_grbi[i] = &(*ub)[i];
        }
    }
    //
    // Model variables and constraints
    GRBVar* vars = model.addVars(*lb_grbi, *ub_grbi, NULL, NULL, NULL, n);
    model.update();
    
    for(i=0; i<k; i++){
        GRBLinExpr lhs = 0;
        for(j=0;j<n;j++){
            if((*A)(i,j) != 0){
                lhs += (*A)(i,j)*vars[j];
            }
        }
        model.addConstr(lhs, sense[i], (*rhs)(i));
    }
    //
    // Model objective
    GRBLinExpr LPobj = 0;
    
    for(j = 0; j < n; j++){
        LPobj += (*obj)(j)*vars[j];
    }
    
    model.setObjective(LPobj,GRB_MAXIMIZE);
    model.update();
    //
    // Optimize and recover relevant solution objects
    model.optimize();
    
    GRBConstr* mycons = model.getConstrs();
    
    double objvalP = 0;
    double* solution = new double[n];
    double* PI = new double[n];
    double* RC = new double[n];
    //double* slack = new double[n];
    
    if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL) {
        objvalP = model.get(GRB_DoubleAttr_ObjVal);
        for (i = 0; i < n; i++){
            solution[i] = vars[i].get(GRB_DoubleAttr_X);
            RC[i] = vars[i].get(GRB_DoubleAttr_RC);
            PI[i] = mycons[i].get(GRB_DoubleAttr_Pi);
            //slack[i] = mycons[i].get(GRB_DoubleAttr_Slack);
        }
        //success = true;
    }
    //
    objval = objvalP;
    sol_vars = *solution;
    sol_PI   = *PI;
    sol_RC   = *RC;
    //
    delete[] lb_grbi;
    delete[] ub_grbi;
    delete[] solution;
    delete[] PI;
    delete[] RC;
    
    //delete env;
}
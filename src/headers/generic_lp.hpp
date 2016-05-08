/*
 * Generic Linear Programming solver using Gurobi
 *
 * Keith O'Hara
 * 05/08/2016
 */

#include "gurobi_c++.h"

bool generic_LP(arma::vec *obj, arma::mat* A, int modelSense, arma::vec* rhs, char* sense, arma::mat* Q, arma::vec* lb, arma::vec* ub, arma::vec* start, double& objval, arma::mat& sol_mat)
{
    // Initialize
    bool success = false;
    
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
        for(i=0; i < n; i++){
            lb_grbi[i] = &(*lb)[i];
        }
    }
    
    if(ub==NULL){
        *ub_grbi = NULL;
    }else{
        for(i=0; i < n; i++){
            ub_grbi[i] = &(*ub)[i];
        }
    }
    //
    // Model variables and constraints
    GRBVar* vars = model.addVars(*lb_grbi, *ub_grbi, NULL, NULL, NULL, n);
    model.update();
    
    for(i=0; i < k; i++){
        GRBLinExpr lhs = 0;
        for(j=0; j < n; j++){
            if((*A)(i,j) != 0){
                lhs += (*A)(i,j)*vars[j];
            }
        }
        model.addConstr(lhs, sense[i], (*rhs)(i));
    }
    //
    // Model objective
    if(Q==NULL){ // Linear programming problem
        GRBLinExpr LPobj = 0;
        
        for(j = 0; j < n; j++){
            LPobj += (*obj)(j)*vars[j];
        }
        
        if(modelSense==1){
            model.setObjective(LPobj,GRB_MAXIMIZE);
        }else{
            model.setObjective(LPobj,GRB_MINIMIZE);
        }
    }else{ // Quadratic programming problem
        GRBQuadExpr QPobj = 0;
        
        for(j = 0; j < n; j++){
            QPobj += (*obj)(j)*vars[j];
        }
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                QPobj += (*A)(i,j)*vars[i]*vars[j];
            }
        }
        
        if(modelSense==1){
            model.setObjective(QPobj,GRB_MAXIMIZE);
        }else{
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
        for (i = 0; i < n; i++){
            sol_mat(i,0) = vars[i].get(GRB_DoubleAttr_X);
            sol_mat(i,1) = vars[i].get(GRB_DoubleAttr_RC);
            sol_mat(i,2) = mycons[i].get(GRB_DoubleAttr_Pi);
        }
        success = true;
    }
    //
    delete[] lb_grbi;
    delete[] ub_grbi;
    //
    return success;
}
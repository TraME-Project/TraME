/*
 * Generic Linear/Quadratic Programming solver using Gurobi
 *
 * Keith O'Hara
 * 05/08/2016
 */

#if defined(WIN32) || defined(_WIN32) || defined(TRAME_USE_GUROBI_C)
# define PREDEF_GUROBI_C
#endif

#if !defined(PREDEF_GUROBI_C)

    #include "gurobi_c++.h"

    bool generic_LP(arma::vec *obj, arma::mat* A, int modelSense, arma::vec* rhs, char* sense, arma::mat* Q, arma::vec* lb, arma::vec* ub, arma::vec* start, double& objval, arma::mat& sol_mat, arma::mat& dual_mat)
    {
        // Initialize
        bool success = false;
        
        int i,j;
        int k = A->n_rows;   // number of constraints
        int n = obj->n_elem; // number of variables
        //
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
                    if((*Q)(i,j) != 0){
                        QPobj += (*Q)(i,j)*vars[i]*vars[j];
                    }
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
            for(i = 0; i < n; i++){
                sol_mat(i,0) = vars[i].get(GRB_DoubleAttr_X);
                sol_mat(i,1) = vars[i].get(GRB_DoubleAttr_RC);
            }
            for(j = 0; j < k; j++){
                dual_mat(j,0) = mycons[j].get(GRB_DoubleAttr_Pi);
                dual_mat(j,1) = mycons[j].get(GRB_DoubleAttr_Slack);
            }
            success = true;
        }
        //
        delete[] lb_grbi;
        delete[] ub_grbi;
        //
        return success;
    }
#else
    #ifndef SWITCH_GRB_ROWCOL_ORDER
        #define SWITCH_GRB_ROWCOL_ORDER
    #endif

    extern "C" {
        #include "generic_lp.h"
    }
    
    bool generic_LP(arma::vec *obj, arma::mat* A, int modelSense, arma::vec* rhs, char* sense, arma::mat* Q, arma::vec* lb, arma::vec* ub, arma::vec* start, double& objval, arma::mat& sol_mat, arma::mat& dual_mat)
    {
        // Initialize
        bool success = false;
        int solved;
        
        int i,j;
        int k = A->n_rows;   // number of constraints ('rows')
        int n = obj->n_elem; // number of variables ('columns')
        //
        // Switch to double arrays
        double** obj_grbi = new double*[n];
        double** A_grbi   = new double*[k*n];
        double** rhs_grbi = new double*[k];
        double** Q_grbi   = new double*[n*n];
        double** lb_grbi  = new double*[n];
        double** ub_grbi  = new double*[n];
        
        for(i=0; i < n; i++){
            obj_grbi[i] = &(*obj)[i];
            //std::cout << *obj_grbi[i] << std::endl;
        }
        
        for(i=0; i < k; i++){
            for(j=0; j < n; j++){
                A_grbi[i*k+j] = &(*A)(i,j); // need to switch the order for C
            }
        }
        
        for(i=0; i < k; i++){
            rhs_grbi[i] = &(*rhs)[i];
            //std::cout << *rhs_grbi[i] << std::endl;
        }
        
        if(Q==NULL){
            *Q_grbi = NULL;
        }else{
            for(i=0; i < n; i++){
                for(j=0; j < n; j++){
                    Q_grbi[i*n+j] = &(*Q)(i,j);
                }
            }
        }
        
        if(lb==NULL){
            *lb_grbi = NULL;
        }else{
            for(i=0; i < n; i++){
                lb_grbi[i] = &(*lb)[i];
                //std::cout << *lb_grbi[i] << std::endl;
            }
        }
        
        if(ub==NULL){
            *ub_grbi = NULL;
        }else{
            for(i=0; i < n; i++){
                ub_grbi[i] = &(*ub)[i];
                //std::cout << *ub_grbi[i] << std::endl;
            }
        }
        
        double* sol_1_grbi  = new double[n];
        double* sol_2_grbi  = new double[n];
        double* dual_1_grbi = new double[k];
        double* dual_2_grbi = new double[k];
        //
        // Call C-version of the solver
        solved = generic_LP_C_switch(k, n, *obj_grbi, *A_grbi, modelSense, *rhs_grbi, sense, *Q_grbi, *lb_grbi,
                                     *ub_grbi, &objval, sol_1_grbi, sol_2_grbi, dual_1_grbi, dual_2_grbi);
        
        //
        // Put solution matrices together
        if (solved == 1) {
            for(i = 0; i < n; i++){
                sol_mat(i,0) = sol_1_grbi[i];
                sol_mat(i,1) = sol_2_grbi[i];
            }
            for(j = 0; j < k; j++){
                dual_mat(j,0) = dual_1_grbi[j];
                dual_mat(j,1) = dual_2_grbi[j];
            }
            success = true;
        }
        //
        delete[] obj_grbi;
        delete[] A_grbi;
        delete[] rhs_grbi;
        delete[] Q_grbi;
        delete[] lb_grbi;
        delete[] ub_grbi;
        
        delete[] sol_1_grbi;
        delete[] sol_2_grbi;
        delete[] dual_1_grbi;
        delete[] dual_2_grbi;
        //
        return success;
    }
#endif
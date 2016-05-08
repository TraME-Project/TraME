/*
 * generic_LP test 2: Quadratic programming
 * Based on Gurobi's 'dense_c++.cpp' example
 *
 * cd ~/Desktop/"Google Drive"/GitHub/TraME/src/tests
 * clang++ -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include generic_lp_2.cpp -o generic_lp_2.test -L/Library/gurobi650/mac64/lib -lgurobi_c++ -lgurobi65 -framework Accelerate
 */

#include "armadillo"

#include "../headers/generic_lp.hpp"

int main()
{
    arma::mat Q(3,3);
    Q  << 1.0 << 1.0 << 0.0 << arma::endr
       << 0.0 << 1.0 << 1.0 << arma::endr
       << 0.0 << 0.0 << 1.0 << arma::endr;
    
    arma::mat A(2,3);
    A  << 1.0 << 2.0 << 3.0 << arma::endr
       << 1.0 << 1.0 << 0.0 << arma::endr;
    
    arma::vec obj(3,1);
    
    obj << 1.0 << arma::endr
        << 1.0 << arma::endr
        << 0.0 << arma::endr;

    arma::vec rhs(2,1);
    
    rhs << 4.0 << arma::endr
        << 1.0 << arma::endr;
        
    char sense[] = {'>', '>'};
    
    arma::vec lb(3,1);
    
    lb  << 0.0 << arma::endr
        << 0.0 << arma::endr
        << 0.0 << arma::endr;
    
    arma::vec ub(3,1);
    ub  << 100.0 << arma::endr
        << 100.0 << arma::endr
        << 100.0 << arma::endr;
    
    bool LP_optimal;
    int modelSense = 0; // minimize
    double objval;
    
    arma::mat sol_mat(obj.n_elem,2);
    arma::mat dual_mat(A.n_rows,2);
    
    try {
        LP_optimal = generic_LP(&obj, &A, modelSense, &rhs, sense, &Q, &lb, &ub, NULL, objval, sol_mat, dual_mat);
        
        std::cout << "\nOptimal value: " << objval << ".\n" << std::endl;
        arma::cout << "Solution: [vars, RC] \n" << sol_mat << arma::endl;
        arma::cout << "Dual Variables: [Pi, Slack] \n" << dual_mat << arma::endl;
    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //  
    return 0;
}
/*
 * generic_LP test: Linear programming
 * Test the ifdef features
 *
 * Keith O'Hara
 * 05/17/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/lp/tests
 * clang++ -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include c_with_cpp_test_2.cpp -c -o c_with_cpp_test_2.o
 * clang++ -O2 -Wall -o cpptest_2.test ../generic_lp.o c_with_cpp_test_2.o -L/Library/gurobi650/mac64/lib -lgurobi65 -framework Accelerate
 */

#include "armadillo"

#define TRAMETESTLP

#include "../generic_lp.hpp"

int main()
{
    arma::mat A(3,3);
    A  << 3.0 << 4.0 << 2.0 << arma::endr
       << 2.0 << 1.0 << 2.0 << arma::endr
       << 1.0 << 3.0 << 2.0 << arma::endr;
    
    arma::vec obj(3,1);
    
    obj << 2.0 << arma::endr
        << 4.0 << arma::endr
        << 3.0 << arma::endr;

    arma::vec rhs(3,1);
    
    rhs << 60.0 << arma::endr
        << 40.0 << arma::endr
        << 80.0 << arma::endr;
        
    char sense[] = {'<', '<', '<'};
    
    arma::vec lb(3,1);
    
    lb  << 0.0 << arma::endr
        << 0.0 << arma::endr
        << 0.0 << arma::endr;
    
    arma::vec ub(3,1);
    ub  << 100.0 << arma::endr
        << 100.0 << arma::endr
        << 100.0 << arma::endr;
    
    bool LP_optimal;
    int modelSense = 1; // maximize
    double objval;
    
    arma::mat sol_mat(obj.n_elem,2);
    arma::mat dual_mat(A.n_rows,2);
    //
    LP_optimal = generic_LP(&obj, &A, modelSense, &rhs, sense, NULL, &lb, NULL, NULL, objval, sol_mat, dual_mat);
    //
    if(LP_optimal){
        std::cout << "\nOptimal value: " << objval << ".\n" << std::endl;
        arma::cout << "Solution: [vars, RC] \n" << sol_mat << arma::endl;
        arma::cout << "Dual Variables: [Pi, Slack] \n" << dual_mat << arma::endl;
    }
    
    return 0;
}
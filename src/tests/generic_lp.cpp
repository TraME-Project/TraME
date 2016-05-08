//
// generic_LP test
//
// cd ~/Desktop/"Google Drive"/GitHub/TraME/src/tests
// clang++ -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include generic_lp.cpp -o generic_lp.test -L/Library/gurobi650/mac64/lib -lgurobi_c++ -lgurobi65 -framework Accelerate
//

#include "armadillo"

#include "../headers/generic_lp.hpp"

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
    
    arma::mat sol_mat(rhs.n_elem,3);
    
    try {
        LP_optimal = generic_LP(&obj, &A, modelSense, &rhs, sense, NULL, &lb, &ub, NULL, objval, sol_mat);
        
        //std::cout << "\n" << std::endl;
        std::cout << "\nOptimal value: " << objval << ".\n" << std::endl;
        arma::cout << "Solution: [vars, Pi, RC] \n" << sol_mat << arma::endl;
    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //  
    return 0;
}
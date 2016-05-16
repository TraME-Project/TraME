/*
 * cd ~/Desktop/"Google Drive"/GitHub/TraME/src/tests
 * clang++ -O2 -Wall -std=c++11 -I/opt/local/include -I/Library/gurobi650/mac64/include arums_probit_test.cpp -o arums_probit.test -L/Library/gurobi650/mac64/lib -lgurobi_c++ -lgurobi65 -framework Accelerate
 */

#include "armadillo"

#include "../headers/trame_structs.hpp"

#include "../headers/arums_empirical.hpp"
#include "../headers/arums_probit.hpp"

int main()
{
    arma::mat U(2,3);
    U  << 1.6 << 3.2 << 1.1 << arma::endr
       << 2.9 << 1.0 << 3.1 << arma::endr;
    
    arma::mat mu(2,3);
    mu << 1.0 << 3.0 << 1.0 << arma::endr
       << 2.0 << 1.0 << 3.0 << arma::endr;
    
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;
    //
    int nbX = U.n_rows;
    int nbY = U.n_cols;
    
    arma::vec n = arma::sum(mu,1);
    //
    probit probits;
    probits.build(nbX,nbY,(bool) true);
    
    probits.rho = 0.5;
    //
    probits.unifCorrelCovMatrices();
    arma::cout << probits.Covar << arma::endl;
    
    int n_draws = 1000;
    empirical emp_obj;
    
    probits.simul(emp_obj, n_draws, (int) 777);
    
    emp_obj.U = U;

    emp_obj.G(n);

    arma::cout << "G emp -> mu: \n" << emp_obj.mu << arma::endl;
    
    arma::mat Ustar_emp;
    emp_obj.Gstar(n, Ustar_emp);
    
    arma::cout << "\\nabla G*(\\nabla G(U emp)): \n" << Ustar_emp << arma::endl;
    //
    return 0;
}
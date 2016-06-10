/*
 * probit class test
 *
 * Keith O'Hara
 * 05/17/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/tests
 * clang++ -O2 -Wall -std=c++11 -I/opt/local/include -I/Library/gurobi650/mac64/include arums_probit_test.cpp -o arums_probit.test -L/Library/gurobi650/mac64/lib -lgurobi_c++ -lgurobi65 -framework Accelerate
 *
 * gcc-mp-5 -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include ../lp/generic_lp.c -c -o ../lp/generic_lp.o
 * g++-mp-5 -O2 -Wall -std=c++11 -fopenmp -I/opt/local/include -I/Library/gurobi650/mac64/include arums_probit_test.cpp -c -o arums_probit_test.o
 * g++-mp-5 -O2 -Wall -fopenmp -o arums_probit.test ../lp/generic_lp.o arums_probit_test.o -L/Library/gurobi650/mac64/lib -lgurobi65 -framework Accelerate
 */

#ifndef __clang__
#define TRAME_USE_GUROBI_C
#endif

#include "armadillo"

#include "../headers/trame_structs.hpp"

#include "../headers/arums_empirical.hpp"
#include "../headers/arums_probit.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs
    arma::mat U(2,3);
    U  << 1.6 << 3.2 << 1.1 << arma::endr
       << 2.9 << 1.0 << 3.1 << arma::endr;
    
    arma::mat mu(2,3);
    mu << 1.0 << 3.0 << 1.0 << arma::endr
       << 2.0 << 1.0 << 3.0 << arma::endr;
    
    //
    // results
    printf("\n*===================   Start of Probit Test   ===================*\n");
    printf("\n");
    printf("Inputs: \n");
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;
    //
    // setup probit class object
    int nbX = U.n_rows;
    int nbY = U.n_cols;
    
    arma::vec n = arma::sum(mu,1);
    
    probit probits;
    probits.build(nbX,nbY,(bool) true);
    //
    // correlation matrices
    probits.rho = 0.5;
    
    probits.unifCorrelCovMatrices();
    arma::cout << "Correlation matrices: \n" << probits.Covar << arma::endl;
    //
    // empirical object:
    int n_draws = 1000;
    empirical emp_obj;
    
    probits.simul(emp_obj, n_draws, (int) 1777);
    
    emp_obj.U = U;
    emp_obj.mu = mu;
    //
    // first compute optimal assignment (mu)
    double G_sim_val = emp_obj.G(n);
    
    arma::cout << "G-sim(U): \n" << G_sim_val << arma::endl;
    arma::cout << "G-sim -> mu: \n" << emp_obj.mu_sol << arma::endl;
    //
    // solution to dual problem U*
    arma::mat U_star_sim;
    double Gstar_sim_val = emp_obj.Gstar(U_star_sim, n);
    
    arma::cout << "G*-sim(mu): \n" << Gstar_sim_val << arma::endl;
    arma::cout << "\\nabla G-sim*(\\nabla G-sim(U)): \n" << U_star_sim << arma::endl;
    //
    // Gbar
    arma::mat mu_bar(2,3);
    mu_bar.fill(2);
    
    arma::mat U_bar_temp, mu_bar_temp;

    double val_Gbar_sim = emp_obj.Gbar(U_star_sim, mu_bar, n, U_bar_temp, mu_bar_temp);
    
    arma::cout << "Gbar-sim val: \n" << val_Gbar_sim << arma::endl;
    //
    printf("\n*===================   End of Probit Test   ===================*\n");
    printf("\n");
    //
    end = std::chrono::system_clock::now();
        
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        
    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
    //
    return 0;
}
/*
 * 'none' class test
 *
 * Keith O'Hara
 * 05/17/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/tests
 *
 * clang++ -O2 -Wall -std=c++11 -I/opt/local/include -I/Library/gurobi650/mac64/include arums_none_test.cpp -o arums_none.test -L/Library/gurobi650/mac64/lib -lgurobi_c++ -lgurobi65 -framework Accelerate
 *
 * gcc-mp-5 -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include ../lp/generic_lp.c -c -o ../lp/generic_lp.o
 * g++-mp-5 -O2 -Wall -std=c++11 -fopenmp -I/opt/local/include -I/Library/gurobi650/mac64/include arums_none_test.cpp -c -o arums_none_test.o
 * g++-mp-5 -O2 -Wall -fopenmp -o arums_none.test ../lp/generic_lp.o arums_none_test.o -L/Library/gurobi650/mac64/lib -lgurobi65 -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs:
    arma::mat U(2,3);
    U  << 1.6 << 3.2 << 1.1 << arma::endr
       << 2.9 << 1.0 << 3.1 << arma::endr;
    
    arma::mat mu(2,3);
    mu << 1.0 << 3.0 << 1.0 << arma::endr
       << 2.0 << 1.0 << 3.0 << arma::endr;
    
    //
    // results
    printf("\n*===================   Start of None Test   ===================*\n");
    printf("\n");
    printf("Inputs: \n");
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;
    //
    // setup none class object
    int nbX = U.n_rows;
    int nbY = U.n_cols;
    
    arma::vec n = arma::sum(mu,1);
    //
    none none_obj;
    none_obj.nbX = nbX;
    none_obj.nbY = nbY;
    none_obj.nbParams = 0;
    
    none_obj.U = U;
    none_obj.mu = mu;
    //
    double val_G = none_obj.G(n);
    std::cout << "G val: \n" << val_G << std::endl;
    //
    arma::mat mubar(2,3);
    mubar.fill(2);

    arma::mat Ubar_temp, mubar_temp;
    double valGbar = none_obj.Gbar(U,mubar,n,Ubar_temp,mubar_temp);
    
    std::cout << "Gbar val: \n" << valGbar << std::endl;
    //
    printf("\n*===================   End of None Test   ===================*\n");
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
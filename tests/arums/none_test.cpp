/*
 * 'none' class test
 *
 * Keith O'Hara
 * 05/17/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/arums
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/usr/local/include/trame none_test.cpp -o none.test -L/usr/local/lib -ltrame -framework Accelerate
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
    
    int nbX = U.n_rows;
    int nbY = U.n_cols;
    
    arma::vec n = arma::sum(mu,1);
    //
    // results
    printf("\n*===================   Start of None Test   ===================*\n");
    printf("\n");
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;
    //
    // setup
    trame::arums::none none_obj(nbX,nbY);

    // first compute optimal assignment (mu)
    arma::mat mu_sol;
    //
    double G_val = none_obj.G(n,U,mu_sol);
    
    std::cout << "G(U) and G-sim(U): \n" << G_val << std::endl;
    arma::cout << "\nG -> mu: \n" << mu_sol << arma::endl;
    //
    arma::mat mubar(2,3);
    mubar.fill(2);

    arma::mat Ubar_temp, mubar_temp;
    double val_Gbar = none_obj.Gbar(U,mubar,n,Ubar_temp,mubar_temp);
    
    std::cout << "Gbar val: \n" << val_Gbar << std::endl;
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
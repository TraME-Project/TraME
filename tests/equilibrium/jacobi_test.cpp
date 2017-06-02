/*
 * jacobi solver test
 *
 * Keith O'Hara
 * 08/25/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/equilibrium
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I./../../include jacobi_test.cpp -o jacobi.test -L./../../ -ltrame -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs:
    arma::mat alpha(2,3);
    alpha << 1.6 << 3.2 << 1.1 << arma::endr 
          << 2.9 << 1.0 << 3.1 << arma::endr;

    arma::mat gamma(2,3);
    gamma << 1.6 << 3.2 << 1.1 << arma::endr 
          << 2.9 << 1.0 << 3.1 << arma::endr;

    arma::mat mu_hat(2,3);
    mu_hat << 1.0 << 3.0 << 1.0 << arma::endr 
           << 2.0 << 1.0 << 3.0 << arma::endr;

    arma::vec n = 1.2 * arma::sum(mu_hat,1);
    arma::vec m = 1.3 * arma::trans(arma::sum(mu_hat,0));

    int nbX = n.n_elem;
    int nbY = m.n_elem;

    arma::mat phi = alpha + gamma;
    //
    // results
    printf("\n*===================   Start of Jacobi Test   ===================*\n");
    printf("\n");
    //
    // TU
    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu> dse_obj_TU;

    trame::arums::logit logit_1(nbX,nbY), logit_2(nbY,nbX);

    dse_obj_TU.build(n,m,phi,logit_1,logit_2,false);
    //
    arma::vec mux0, mu0y;
    arma::mat mu_TU, U, V;
    //trame::jacobi(dse_obj_TU, true, NULL, NULL, NULL, mu_TU, mux0, mu0y, U, V);
    trame::jacobi(dse_obj_TU, mu_TU);

    arma::cout << "Solution of TU-logit problem using jacobi:\n" << mu_TU << arma::endl;
    //
    // NTU
    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::ntu> dse_obj_NTU;

    dse_obj_NTU.build(n,m,alpha,gamma,logit_1,logit_2,false);
    //
    arma::mat mu_NTU;
    //trame::jacobi(dse_obj_NTU, true, NULL, NULL, NULL, mu_NTU, mux0, mu0y, U, V);
    trame::jacobi(dse_obj_NTU, mu_NTU);

    arma::cout << "Solution of NTU-logit problem using jacobi:\n" << mu_NTU << arma::endl;
    //
    printf("\n*===================    End of Jacobi Test    ===================*\n");
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
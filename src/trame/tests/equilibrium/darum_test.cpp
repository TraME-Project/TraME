/*
 * DSE class test
 *
 * Keith O'Hara
 * 08/01/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/equilibrium
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include darum_test.cpp -o darum.test -L/opt/local/lib -ltrame -framework Accelerate
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
    //
    // results
    printf("\n*===================   Start of darum Test   ===================*\n");
    printf("\n");
    //
    // NTU
    trame::dse<trame::logit> dse_obj_NTU;

    trame::logit logit_1(nbX,nbY), logit_2(nbY,nbX);
    //logit_1.build(nbX,nbY,1.0,true);
    //logit_2.build(nbY,nbX,1.0,true);

    dse_obj_NTU.build_NTU(n,m,alpha,gamma,logit_1,logit_2,false);
    //
    arma::mat mu_NTU, U, V;
    //trame::darum(dse_obj_NTU, true, NULL, mu_NTU, mux0, mu0y, U, V);
    trame::darum(dse_obj_NTU, mu_NTU);

    std::cout << "Solution of NTU-logit problem using darum:\n" << mu_NTU << std::endl;
    //
    printf("\n*===================    End of DSE Test    ===================*\n");
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
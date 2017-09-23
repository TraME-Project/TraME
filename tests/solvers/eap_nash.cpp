/*
 * EAP Nash solver test
 *
 * Keith O'Hara
 * 08/01/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/equilibrium
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include eap_nash_test.cpp -o eap_nash.test -L/opt/local/lib -ltrame -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs:
    int nbX = 8;
    int nbY = 5;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha  = arma::randu(nbX,nbY);
    arma::mat gamma  = arma::randu(nbX,nbY);
    arma::mat lambda = 1 + arma::randu(nbX,nbY);
    arma::mat zeta   = 1 + arma::randu(nbX,nbY);

    arma::mat phi = alpha + gamma;
    //
    // results
    printf("\n*===================   Start of eap_nash Test   ===================*\n");
    printf("\n");
    //
    // TU
    arma::mat lambda_LTU = lambda/(lambda+zeta);
    arma::mat phi_LTU = (lambda%alpha + zeta%gamma) / (lambda+zeta);

    trame::dse<trame::arums::none,trame::arums::none,trame::transfers::ltu> dse_obj_LTU;
    dse_obj_LTU.build(n,m,lambda_LTU,phi_LTU,false);
    //
    arma::vec mux0, mu0y, u, v;
    arma::mat mu_LTU;
    trame::eap_nash(dse_obj_LTU,mu_LTU);

    std::cout << "Solution of TU-none problem using eap_nash:\n" << std::endl;
    arma::cout << "mu:\n" << mu_LTU << arma::endl;
    //arma::cout << "u:\n" << u << "\n v:\n" << v << arma::endl;
    //
    printf("\n*===================    End of eap_nash Test    ===================*\n");
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
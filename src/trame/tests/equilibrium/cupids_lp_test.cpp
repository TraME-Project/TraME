/*
 * EAP Nash solver test
 *
 * Keith O'Hara
 * 08/01/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/equilibrium
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include eap_nash_test.cpp -c -o eap_nash_test.o
 * g++-mp-5 -O2 -Wall -o eap_nash.test eap_nash_test.o -L/opt/local/lib -ltrame -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs:
    int nbX = 5;
    int nbY = 3;

    int nbDraws = 100;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha  = arma::randu(nbX,nbY);
    arma::mat gamma  = arma::randu(nbX,nbY);

    arma::mat phi = alpha + gamma;
    //
    // results
    printf("\n*===================   Start of cupids_lp Test   ===================*\n");
    printf("\n");
    //
    // TU
    trame::dse<trame::logit> dse_obj_TU;

    trame::logit logit_1, logit_2;
    logit_1.build(nbX,nbY,1.0,true);
    logit_2.build(nbY,nbX,1.0,true);

    trame::empirical logit_sim_1;
    trame::empirical logit_sim_2;

    logit_1.simul

    dse_obj_TU.build_TU(n,m,phi,logit_1,logit_2,false);
    //
    arma::vec mux0, mu0y;
    arma::mat mu_TU, U, V;

    //trame::cupids_lp(dse_obj_LTU, true, NULL, mu_TU, mux0, mu0y, U, V);

    //std::cout << "Solution of TU-none problem using eap_nash:\n" << std::endl;
    //arma::cout << "mu:\n" << mu_LTU << arma::endl;
    //arma::cout << "u:\n" << u << "\n v:\n" << v << arma::endl;
    //
    printf("\n*===================    End of cupids_lp Test    ===================*\n");
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
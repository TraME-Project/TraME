/*
 * Cupids-LP solver test
 *
 * Keith O'Hara
 * 08/01/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/equilibrium
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include cupids_lp_test.cpp -o cupids_lp.test -L/opt/local/lib -ltrame -framework Accelerate
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

    int nb_draws_1 = 1000;
    int nb_draws_2 = 1000;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha  = arma::randu(nbX,nbY);
    arma::mat gamma  = arma::randu(nbX,nbY);

    arma::mat phi = alpha + gamma;
    
    //
    // results

    printf("\n*===================   Start of cupids_lpTest   ===================*\n");
    printf("\n");
    
    //
    // build

    trame::dse<trame::arums::empirical,trame::arums::empirical,trame::transfers::tu> dse_obj_TU;

    trame::arums::logit logit_1, logit_2;
    logit_1.build(nbX,nbY,1.0,true);
    logit_2.build(nbY,nbX,1.0,true);

    trame::arums::empirical logit_sim_1 = logit_1.simul(nb_draws_1,123);
    trame::arums::empirical logit_sim_2 = logit_2.simul(nb_draws_2,321);

    dse_obj_TU.build(n,m,phi,logit_sim_1,logit_sim_2,false);
    
    //

    arma::vec mu_x0_out, mu_0y_out;
    arma::mat mu_TU, U_out, V_out;

    trame::cupids_lp(dse_obj_TU,mu_TU);
    trame::cupids_lp(dse_obj_TU,mu_TU,U_out,V_out);
    trame::cupids_lp(dse_obj_TU,mu_TU,mu_x0_out,mu_0y_out,U_out,V_out);

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

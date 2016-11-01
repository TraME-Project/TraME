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

    int nbDraws_1 = 1000;
    int nbDraws_2 = 1000;

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
    // TU logit
    trame::dse<trame::empirical> dse_obj_TU;

    trame::logit logit_1, logit_2;
    logit_1.build(nbX,nbY,1.0,true);
    logit_2.build(nbY,nbX,1.0,true);

    trame::empirical logit_sim_1 = logit_1.simul(&nbDraws_1,NULL);
    trame::empirical logit_sim_2 = logit_2.simul(&nbDraws_2,NULL);

    dse_obj_TU.build_TU(n,m,phi,logit_sim_1,logit_sim_2,false);
    //
    arma::vec mux0, mu0y;
    arma::mat mu_TU, U, V;

    trame::cupids_lp(dse_obj_TU, mu_TU);

    std::cout << "Solution of TU-logitSim problem using LP:\n" << std::endl;
    arma::cout << "mu:\n" << mu_TU << arma::endl;
    //
    // TU probit
    trame::dse<trame::empirical> dse_obj_TU_2;

    double rho = 0.3;

    trame::probit probit_1, probit_2;
    probit_1.build(nbX,nbY,rho,true);
    probit_2.build(nbY,nbX,rho,true);

    probit_1.unifCorrelCovMatrices();
    probit_2.unifCorrelCovMatrices();

    trame::empirical probit_sim_1 = probit_1.simul(&nbDraws_1,NULL);
    trame::empirical probit_sim_2 = probit_2.simul(&nbDraws_2,NULL);

    dse_obj_TU_2.build_TU(n,m,phi,probit_sim_1,probit_sim_2,false);
    //
    arma::mat mu_TU_2;

    trame::cupids_lp(dse_obj_TU_2, mu_TU_2);

    std::cout << "Solution of TU-probitSim problem using LP:\n" << std::endl;
    arma::cout << "mu:\n" << mu_TU_2 << arma::endl;
    //arma::cout << "u:\n" << u << "\n v:\n" << v << arma::endl;*/
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

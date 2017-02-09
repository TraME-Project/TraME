/*
 * arc_newton test
 *
 * Keith O'Hara
 * 02/07/2017
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/equilibrium
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include arc_newton_test.cpp -o arc_newton.test -L/opt/local/lib -ltrame -framework Accelerate
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
    double sigma = 1;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat phi = arma::randu(nbX,nbY);
    //
    // results
    printf("\n*===================   Start of arc_newton Test   ===================*\n");
    printf("\n");
    //
    // TU
    trame::logit logit_1(nbX,nbY), logit_2(nbY,nbX);
    
    trame::mfe<trame::mmf> mfe_obj_TU;
    trame::dse<trame::logit> dse_obj_TU;

    mfe_obj_TU.build_TU(n,m,phi,&sigma,false);
    dse_obj_TU.build_TU(n,m,phi,logit_1,logit_2,false);
    //
    //
    //double tol = 1E-06;
    //int max_iter = 5000;

    arma::mat mu_TU_1, mu_TU_2;
    trame::ipfp(mfe_obj_TU,mu_TU_1);

    //trame::max_welfare(dse_obj_TU,mu_TU_2);
    //dse_obj_TU.solve(mu_TU_2, (char*) "maxWelfare");
    std::cout << "begin arc_newton" << std::endl;
    trame::arc_newton(dse_obj_TU,mu_TU_2);

    arma::cout << "Solution of TU-logit problem using ipfp:\n" << mu_TU_1 << arma::endl;
    arma::cout << "Solution of TU-logit problem using arc_newton:\n" << mu_TU_2 << arma::endl;
    //
    printf("\n*===================    End of arc_newton Test    ===================*\n");
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
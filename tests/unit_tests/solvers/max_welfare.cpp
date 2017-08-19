/*
 * max_welfare test
 *
 * Keith O'Hara
 * 24/10/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/equilibrium
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include max_welfare_test.cpp -o max_welfare.test -L/opt/local/lib -ltrame -framework Accelerate
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

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat phi = arma::randu(nbX,nbY);

    //
    // results

    printf("\n*===================   Start of max_welfare Test   ===================*\n");
    printf("\n");

    //
    // build

    trame::arums::logit logit_1(nbX,nbY), logit_2(nbY,nbX);
    
    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu> dse_obj_TU;

    dse_obj_TU.build(n,m,phi,logit_1,logit_2,false);
    
    //

    double tol = 1E-06;
    int max_iter = 5000;

    arma::vec mux0, mu0y;
    arma::mat mu_TU, U, V;
    
    trame::max_welfare(dse_obj_TU,mu_TU);
    trame::max_welfare(dse_obj_TU,mu_TU,tol);
    trame::max_welfare(dse_obj_TU,mu_TU,max_iter);
    trame::max_welfare(dse_obj_TU,mu_TU,tol,max_iter);

    trame::max_welfare(dse_obj_TU,mu_TU,U,V);

    double val_out;

    trame::max_welfare(dse_obj_TU,mu_TU,mux0,mu0y,U,V,val_out,&tol,&max_iter);

    
    //
    printf("\n*===================    End of max_welfare Test    ===================*\n");
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
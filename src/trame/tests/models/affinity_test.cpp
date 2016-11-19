/*
 * affinity model class test
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include affinity_test.cpp -o affinity.test -L/opt/local/lib -ltrame -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs:
    int nbX = 100;
    int nbY = 100;

    arma::mat X_vals = arma::randu(nbX,10);
    arma::mat Y_vals = arma::randu(nbY,10);
    //
    // results
    printf("\n*===================   Start of affinity model Test   ===================*\n");
    printf("\n");
    //
    // TU
    trame::affinity aff_model;
    aff_model.build(X_vals,Y_vals);

    double lambda = 0.15;
    double val_hat;
    arma::mat mu_hat = arma::ones(nbX,nbY), theta_hat;

    aff_model.mme_regul(mu_hat,lambda,theta_hat,val_hat,NULL,NULL,NULL,NULL);

    //
    printf("\n*===================    End of affinity model Test    ===================*\n");
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
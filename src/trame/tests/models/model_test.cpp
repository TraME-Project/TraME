/*
 * general model class test
 *
 * Keith O'Hara
 * 11/27/2016
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include model_test.cpp -o model.test -L/opt/local/lib -ltrame -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs:
    int nbX = 5;
    int nbY = 4;
    int dX = 3;
    int dY = 2;

    arma::mat X_vals = arma::randu(nbX,dX);
    arma::mat Y_vals = arma::randu(nbY,dY);
    //
    // results
    printf("\n*===================   Start of general model Test   ===================*\n");
    printf("\n");
    //
    // TU
    trame::model<trame::logit> model_obj;
    model_obj.build(X_vals,Y_vals);

    //double val_hat;
    arma::mat mu_hat = arma::ones(nbX,nbY), theta_hat;

    /*aff_model.mme_regul(mu_hat,lambda,theta_hat,val_hat,NULL,NULL,NULL,NULL);*/

    //
    printf("\n*===================    End of general model Test    ===================*\n");
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
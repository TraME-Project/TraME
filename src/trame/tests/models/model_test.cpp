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
    int nbX = 6;
    int nbY = 5;
    int dX = 3;
    int dY = 3;

    arma::mat X_vals = arma::randu(nbX,dX);
    arma::mat Y_vals = arma::randu(nbY,dY);

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);
    //
    // results
    printf("\n*===================   Start of general model Test   ===================*\n");
    printf("\n");
    //
    // TU
    double noise_scale = 0.1, sigma = 1.0;

    arma::mat A = arma::ones(dX,dY);
    arma::mat phi = X_vals*A*Y_vals.t();
    trame::mfe<trame::mmf> mfe_obj_TU;
    mfe_obj_TU.build_TU(n,m,phi,&sigma,false);
    arma::mat noise = 1.0 + noise_scale*arma::randn(nbX,nbY);

    arma::mat mu_mfe;
    mfe_obj_TU.solve(mu_mfe);

    arma::mat mu_hat = mu_mfe % noise;
    //
    trame::model<trame::logit> TU_logit_model;
    TU_logit_model.build(X_vals,Y_vals,n,m);

    //double val_hat;
    arma::mat theta_hat;
    TU_logit_model.mme(mu_hat,theta_hat);

    arma::cout << "theta_hat mme1: \n" << theta_hat << arma::endl;

    TU_logit_model.mme_2(mu_hat,theta_hat);

    arma::cout << "theta_hat mme2: \n" << theta_hat << arma::endl;

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

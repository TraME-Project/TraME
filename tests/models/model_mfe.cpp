/*
 * general model class test
 *
 * Keith O'Hara
 * 11/27/2016
 *
 * This version:
 * 06/01/2017
 *
 * g++-mp-7 -O2 -Wall -std=c++11 -I/usr/local/include/trame model_test.cpp -o model.test -L/usr/local/lib -ltrame -framework Accelerate
 * g++-mp-7 -O2 -Wall -std=c++11 -I./../../include model.cpp -o model.test -L./../../ -ltrame -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    //
    // inputs:

    const int nbX = 80;
    const int nbY = 80;
    const int dX = 3;
    const int dY = 3;

    arma::mat X_vals = 0.5*arma::randn(nbX,dX);
    arma::mat Y_vals = 0.5*arma::randn(nbY,dY);

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    //
    // results

    printf("\n*===================   Start of general model Test   ===================*\n");
    printf("\n");

    //
    // build

    trame::model<trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu>> TU_logit_model;
    TU_logit_model.build(X_vals,Y_vals,n,m);

    // generate mu_hat 

    double noise_scale = 0.1; 

    // arma::mat A = arma::ones(dX,dY);
    arma::mat A = 1.0 + arma::randu(dX,dY);
    arma::mat Phi = 2.0*X_vals*A*Y_vals.t(); // KEITH: we multiply by 2 to account for sigma = 1 in the mme code (instead of sigma = 1/2)

    trame::mfe<trame::mmfs::geo> mfe_obj_TU;
    mfe_obj_TU.build(n,m,Phi);

    //

    arma::mat mu_mfe;
    mfe_obj_TU.solve(mu_mfe);

    arma::cout << "mu:\n " << mu_mfe(arma::span(0,9),arma::span(0,9));

    // arma::mat noise = 1.0 + noise_scale*arma::randn(nbX,nbY);
    // arma::mat mu_hat = mu_mfe % noise;
    arma::mat mu_hat = mu_mfe;

    //
    // Affinity model

    trame::model<trame::mfe<trame::mmfs::geo>> aff_model;
    aff_model.build(X_vals,Y_vals);

    // renormalize to have sum(row_i) = 1, assuming nbX_row_i = 1

    // mu_hat = arma::ones(nbX,nbY)/ static_cast<double>(nbX);
    // mu_hat = arma::eye(nbX,nbX);

    // mu_hat = mu_hat / static_cast<double>(nbX);
    // mu_hat = mu_hat / arma::accu(mu_hat);
    mu_hat = mu_hat / arma::repmat( arma::sum(mu_hat,1),1,nbY );

    arma::cout << "mu_hat:\n " << mu_hat(arma::span(0,9),arma::span(0,9));

    //

    double lambda = 0.01;
    double val_hat_1, val_hat_2;
    arma::mat theta_hat_aff_1, theta_hat_aff_2;

    aff_model.mme_regul(mu_hat,theta_hat_aff_1,lambda,&val_hat_1,nullptr,nullptr,nullptr,nullptr);
    aff_model.mme_woregul(mu_hat,theta_hat_aff_2,&val_hat_2,nullptr,nullptr,nullptr,nullptr,nullptr);

    std::cout << "\nobjective value with regularization:    " << val_hat_1 << std::endl;
    std::cout << "objective value without regularization: " << val_hat_2 << std::endl;

    arma::cout << "\n Real Affinity Matrix:\n" << A << arma::endl;

    arma::cout << "theta_hat with regularization:\n" << arma::reshape(theta_hat_aff_1,dX,dY) << arma::endl;
    arma::cout << "theta_hat without regularization:\n" << arma::reshape(theta_hat_aff_2,dX,dY) << arma::endl;

    //
    printf("\n*===================    End of general model Test    ===================*\n");
    printf("\n");
    //
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
    //
    return 0;
}

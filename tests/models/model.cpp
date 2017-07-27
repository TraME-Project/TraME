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
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs:
    int nbX = 5;
    int nbY = 4;
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
    double noise_scale = 0.1;

    arma::mat A = arma::ones(dX,dY);
    arma::mat phi = X_vals*A*Y_vals.t();

    trame::mfe<trame::mmfs::geo> mfe_obj_TU;
    // mfe_obj_TU.build(n,m,phi,&sigma,false);
    mfe_obj_TU.build(n,m,phi);
    arma::mat noise = 1.0 + noise_scale*arma::randn(nbX,nbY);

    arma::mat mu_mfe;
    mfe_obj_TU.solve(mu_mfe);

    arma::mat mu_hat = mu_mfe % noise;
    //
    trame::model<trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu>> TU_logit_model;
    TU_logit_model.build(X_vals,Y_vals,n,m);

    //double val_hat;
    arma::mat theta_hat_mme, theta_hat_mle;
    TU_logit_model.mme(mu_hat,theta_hat_mme,nullptr);

    arma::cout << "theta_hat mme: \n" << theta_hat_mme << arma::endl;

    TU_logit_model.mle(mu_hat,theta_hat_mle,nullptr);

    arma::cout << "theta_hat mle: \n" << theta_hat_mle << arma::endl;

    // arums_none
    trame::model< trame::dse<trame::arums::none, trame::arums::none, trame::transfers::tu> > TU_none_model;
    TU_none_model.build(X_vals,Y_vals,n,m);

    // arma::mat theta_0 = arma::zeros(dX*dY,1);
    TU_none_model.mme(mu_hat,theta_hat_mme,nullptr);

    arma::cout << "theta_hat mme for arums_none: \n" << theta_hat_mme << arma::endl;

    // aff_model.mme_regul(mu_hat,lambda,theta_hat,val_hat,nullptr,nullptr,nullptr,nullptr);

    //
    // Affinity model

    nbX = 100;
    nbY = 100;

    X_vals = arma::randu(nbX,10);
    Y_vals = arma::randu(nbY,10);

    trame::model<trame::mfe<trame::mmfs::geo>> aff_model;
    aff_model.build(X_vals,Y_vals);

    double lambda = 0.15;
    double val_hat_1, val_hat_2;
    mu_hat = arma::ones(nbX,nbY)/(nbX);
    arma::mat theta_hat_aff;

    aff_model.mme_regul(mu_hat,lambda,theta_hat_aff,val_hat_1,nullptr,nullptr,nullptr,nullptr);
    std::cout << "the_val with regularization: " << val_hat_1 << std::endl;
    aff_model.mme_woregul(mu_hat,theta_hat_aff,val_hat_2,nullptr,nullptr,nullptr,nullptr);
    std::cout << "the_val without regularization: " << val_hat_2 << std::endl;
    // arma::cout << "theta_hat: \n" << theta_hat << arma::endl;

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

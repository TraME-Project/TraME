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
    const int nbY = 72;
    const int dX = 3;
    const int dY = 3;

    arma::mat X_vals = arma::randu(nbX,dX);
    arma::mat Y_vals = arma::randu(nbY,dY);

    X_vals.load("xs.txt",arma::auto_detect);
    Y_vals.load("ys.txt",arma::auto_detect);

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

    arma::mat A = arma::ones(dX,dY);
    arma::mat phi = X_vals*A*Y_vals.t();

    trame::mfe<trame::mmfs::geo> mfe_obj_TU;
    mfe_obj_TU.build(n,m,phi);

    arma::mat noise = 1.0 + noise_scale*arma::randn(nbX,nbY);

    arma::mat mu_mfe;
    mfe_obj_TU.solve(mu_mfe);

    arma::mat mu_hat = mu_mfe % noise;

    mu_hat.load("mu_hat.txt",arma::auto_detect);

    arma::cout << mu_hat.rows(0,9) << arma::endl;  
    
    // MME and MLE

    arma::mat theta_hat_mme, theta_hat_mle;
    // TU_logit_model.mme(mu_hat,theta_hat_mme,nullptr);

    // arma::cout << "theta_hat mme: \n" << theta_hat_mme << arma::endl; 

    TU_logit_model.mle(mu_hat,theta_hat_mle,nullptr);

    arma::cout << "theta_hat mle: \n" << theta_hat_mle << arma::endl;

    //
    // sim model

    // const int n_draws = 1000;
    // const int sim_seed = 777;

    // trame::arums::logit logit_1(nbX,nbY), logit_2(nbY,nbX);

    // trame::arums::empirical logit_sim_1(nbX,nbY), logit_sim_2(nbY,nbX);
    // logit_1.simul(logit_sim_1, n_draws, sim_seed);
    // logit_2.simul(logit_sim_2, n_draws, sim_seed);

    // trame::model< trame::dse<trame::arums::empirical, trame::arums::empirical, trame::transfers::tu> > logit_sim_model;

    // logit_sim_model.build(X_vals,Y_vals,n,m);
    // logit_sim_model.market_obj.arums_G = logit_sim_1;
    // logit_sim_model.market_obj.arums_H = logit_sim_2;

    // arma::mat theta_hat_3;

    // logit_sim_model.mme(mu_hat,theta_hat_3,nullptr);
    
    // arma::cout << "logit_sim_model theta_hat:\n" << theta_hat_3 << arma::endl;

    // //
    // // arums_none model

    // trame::model< trame::dse<trame::arums::none, trame::arums::none, trame::transfers::tu> > TU_none_model;
    // TU_none_model.build(X_vals,Y_vals,n,m);

    // TU_none_model.mme(mu_hat,theta_hat_mme,nullptr);

    // arma::cout << "theta_hat mme for arums_none: \n" << theta_hat_mme << arma::endl;

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

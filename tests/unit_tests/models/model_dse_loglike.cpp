/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
  ##
  ##   This file is part of TraME.
  ##
  ##   TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * test loglikelihood function using data generated from R test code, where
 * seed=777, nbX=5, nbY=4, noiseScale=0.1, dX=3, dY=3
 *
 * Keith O'Hara
 * 11/27/2016
 *
 * This version:
 * 08/20/2017
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    //
    // inputs:

    const int nbX = 5;
    const int nbY = 4;
    const int dX = 3;
    const int dY = 3;

    arma::mat X_vals(nbX,dX);
    arma::mat Y_vals(nbY,dY);

    X_vals.load("xs_small.txt",arma::auto_detect);
    Y_vals.load("ys_small.txt",arma::auto_detect);

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    //
    // results

    printf("\n*===================   Start of model loglikelihood Test   ===================*\n");
    printf("\n");

    // 
    // build

    trame::model<trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu>> TU_logit_model; 
    TU_logit_model.build(X_vals,Y_vals,n,m);


    arma::mat mu_hat;
    mu_hat.load("mu_hat_small.txt",arma::auto_detect);

    // calculate loglikelihood

    arma::mat theta = arma::ones(dX*dY);

    TU_logit_model.model_to_market(theta);

    const double scale = std::max(arma::accu(n),arma::accu(m));

    const arma::vec mu_hat_x0 = n - arma::sum(mu_hat,1);
    const arma::vec mu_hat_0y = m - arma::trans(arma::sum(mu_hat,0));

    arma::cout << "mu_hat:\n" << mu_hat << arma::endl;
    arma::cout << "mu_hat_x0:\n" << mu_hat_x0 << arma::endl;
    arma::cout << "mu_hat_0y:\n" << mu_hat_0y << arma::endl;
    //
    // add optimization data
    trame::trame_model_mle_opt_data<trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu>> opt_data;

    opt_data.model_obj = TU_logit_model;
    opt_data.by_individual = true;
    opt_data.scale = scale;

    opt_data.mu_hat = mu_hat;
    opt_data.mu_hat_x0 = mu_hat_x0;
    opt_data.mu_hat_0y = mu_hat_0y;

    arma::vec grad_vec;

    double LL_val = TU_logit_model.log_likelihood(theta,&grad_vec,&opt_data);
    std::cout << "LL_val = " << LL_val << std::endl;
    arma::cout << "grad:\n" << grad_vec << arma::endl;

    //
    printf("\n*===================    End of model loglikelihood Test    ===================*\n");
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

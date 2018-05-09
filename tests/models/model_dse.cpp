/*################################################################################
  ##
  ##   Copyright (C) 2015-2018 the TraME Team:
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
 * DSE model class test using data generated from R test code, where
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

    bool small_ver = false;

    int nbX = 80;
    int nbY = 72;

    if (small_ver) {
        nbX = 5;
        nbY = 4;
    }

    const int dX = 3;
    const int dY = 3;

    arma::mat X_vals(nbX,dX);
    arma::mat Y_vals(nbY,dY);

    if (small_ver) {
        X_vals.load("xs_small.txt",arma::auto_detect);
        Y_vals.load("ys_small.txt",arma::auto_detect);
    } else {
        X_vals.load("xs.txt",arma::auto_detect);
        Y_vals.load("ys.txt",arma::auto_detect);
    }

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

    arma::mat mu_hat(nbX,nbY);

    if (small_ver) {
        mu_hat.load("mu_hat_small.txt",arma::auto_detect);
    } else {
        mu_hat.load("mu_hat.txt",arma::auto_detect);
    }

    // MME and MLE

    arma::mat theta_hat_mme, theta_hat_mle;

    TU_logit_model.mme(mu_hat,theta_hat_mme,nullptr);

    arma::cout << "theta_hat mme: \n" << theta_hat_mme << arma::endl; 

    TU_logit_model.mle(mu_hat,theta_hat_mle,nullptr);

    arma::cout << "theta_hat mle: \n" << theta_hat_mle << arma::endl;

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

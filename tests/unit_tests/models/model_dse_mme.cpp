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
 * test MME objective function using data generated from R test code, where
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

    printf("\n*===================   Start of model MME objfn Test   ===================*\n");
    printf("\n");

    // 
    // build

    trame::model<trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu>> TU_logit_model; 
    TU_logit_model.build(X_vals,Y_vals,n,m);


    arma::mat mu_hat;
    mu_hat.load("mu_hat_small.txt",arma::auto_detect);

    // calculate loglikelihood

    arma::mat theta = arma::ones(dX*dY);

    arma::mat dtheta_Psi;
    TU_logit_model.dtheta(nullptr,dtheta_Psi);

    TU_logit_model.model_to_market(theta);

    const arma::mat kron_term = dtheta_Psi;
    const arma::vec C_hat = arma::vectorise(arma::trans(arma::vectorise(mu_hat)) * kron_term);

    arma::cout << "mu_hat:\n" << mu_hat << arma::endl;
    arma::cout << "kron_term:\n" << kron_term << arma::endl;
    arma::cout << "C_hat:\n" << C_hat << arma::endl;

    //
    // add optimization data
        
    trame::trame_model_mme_opt_data<trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu>> opt_data;

    opt_data.market = TU_logit_model.market_obj;
    opt_data.dim_theta = TU_logit_model.dim_theta;
    opt_data.C_hat = C_hat;
    opt_data.kron_term = kron_term;
    //
    arma::vec sol_vec = arma::join_cols(arma::vectorise(kron_term * theta)/2.0,theta);

    arma::vec grad_vec;

    double mme_val = TU_logit_model.model_mme_opt_objfn(sol_vec,&grad_vec,&opt_data);
    std::cout << "mme_val = " << mme_val << std::endl;
    arma::cout << "grad:\n" << grad_vec << arma::endl;

    //
    printf("\n*===================    End of model MME objfn Test    ===================*\n");
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

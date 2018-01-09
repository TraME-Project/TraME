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
 * Cupids-LP solver test
 *
 * Keith O'Hara
 * 10/24/2016
 *
 * This version:
 * 08/18/2017
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

    int n_draws_1 = 1000;
    int n_draws_2 = 1000;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha  = arma::randu(nbX,nbY);
    arma::mat gamma  = arma::randu(nbX,nbY);

    arma::mat phi = alpha + gamma;
    
    //
    // results

    printf("\n*===================   Start of cupids_lpTest   ===================*\n");
    printf("\n");
    
    //
    // build

    trame::dse<trame::arums::empirical,trame::arums::empirical,trame::transfers::tu> dse_obj_TU;

    trame::arums::logit logit_1, logit_2;
    logit_1.build(nbX,nbY,1.0,true);
    logit_2.build(nbY,nbX,1.0,true);

    trame::arums::empirical logit_sim_1 = logit_1.simul(n_draws_1,123);
    trame::arums::empirical logit_sim_2 = logit_2.simul(n_draws_2,321);

    dse_obj_TU.build(n,m,phi,logit_sim_1,logit_sim_2,false);
    
    //

    double val_out;
    arma::vec mu_x0_out, mu_0y_out;
    arma::mat mu_TU, U_out, V_out;

    trame::cupids_lp(dse_obj_TU,mu_TU);
    trame::cupids_lp(dse_obj_TU,mu_TU,U_out,V_out);
    trame::cupids_lp(dse_obj_TU,mu_TU,mu_x0_out,mu_0y_out,U_out,V_out,val_out);

    //
    printf("\n*===================    End of cupids_lp Test    ===================*\n");
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

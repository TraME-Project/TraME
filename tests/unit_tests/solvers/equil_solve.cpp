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
 * ipfp test
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

    int nbX = 18;
    int nbY = 5;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha  = arma::randu(nbX,nbY);
    arma::mat gamma  = arma::randu(nbX,nbY);
    arma::mat lambda = 1 + arma::randu(nbX,nbY);
    arma::mat tau    = arma::randu(nbX,nbY);
    arma::mat zeta   = 1 + arma::randu(nbX,nbY);

    arma::mat phi = alpha + gamma;

    arma::mat lambda_LTU = lambda/(lambda+zeta);
    arma::mat phi_LTU = (lambda%alpha + zeta%gamma) / (lambda+zeta);

    //
    // results

    printf("\n*===================   Start of equil_solve Test   ===================*\n");
    printf("\n");
    
    //
    // build

    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::etu> dse_obj_ETU;
    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::ltu> dse_obj_LTU;
    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::ntu> dse_obj_NTU;
    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu> dse_obj_TU;

    trame::arums::logit logit_1(nbX,nbY), logit_2(nbY,nbX);

    dse_obj_ETU.build(n,m,alpha,gamma,tau,logit_1,logit_2,false);
    dse_obj_LTU.build(n,m,lambda_LTU,phi_LTU,logit_1,logit_2,false);
    dse_obj_NTU.build(n,m,alpha,gamma,logit_1,logit_2,false);
    dse_obj_TU.build(n,m,phi,logit_1,logit_2,false);
    
    //

    char slv = 'j';

    arma::mat mu, U, V;
    
    trame::equil_solve(dse_obj_ETU,mu);
    trame::equil_solve(dse_obj_ETU,mu,nullptr);
    trame::equil_solve(dse_obj_ETU,mu,&slv);

    trame::equil_solve(dse_obj_LTU,mu);
    trame::equil_solve(dse_obj_LTU,mu,nullptr);
    trame::equil_solve(dse_obj_LTU,mu,&slv);

    trame::equil_solve(dse_obj_NTU,mu);
    trame::equil_solve(dse_obj_NTU,mu,nullptr);
    trame::equil_solve(dse_obj_NTU,mu,&slv);

    trame::equil_solve(dse_obj_TU,mu);
    trame::equil_solve(dse_obj_TU,mu,nullptr);
    trame::equil_solve(dse_obj_TU,mu,&slv);

    // trame::jacobi(dse_obj_TU,mu_TU,U,V);

    //
    printf("\n*===================    End of equil_solve Test    ===================*\n");
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
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
 * EAP Nash solver test
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
    int nbX = 8;
    int nbY = 5;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha  = arma::randu(nbX,nbY);
    arma::mat gamma  = arma::randu(nbX,nbY);
    arma::mat lambda = 1 + arma::randu(nbX,nbY);
    arma::mat zeta   = 1 + arma::randu(nbX,nbY);

    arma::mat phi = alpha + gamma;

    //
    // results

    printf("\n*===================   Start of eap_nash Test   ===================*\n");
    printf("\n");

    //
    // build

    arma::mat lambda_LTU = lambda/(lambda+zeta);
    arma::mat phi_LTU = (lambda%alpha + zeta%gamma) / (lambda+zeta);

    trame::dse<trame::arums::none,trame::arums::none,trame::transfers::ltu> dse_obj_LTU;
    dse_obj_LTU.build(n,m,lambda_LTU,phi_LTU,false);

    // 

    double tol = 1E-06;
    int max_iter = 5000;
    bool x_first = true;

    arma::vec mux0, mu0y, u, v;
    arma::mat mu_LTU;

    trame::eap_nash(dse_obj_LTU,mu_LTU);

    trame::eap_nash(dse_obj_LTU,mu_LTU,true);
    trame::eap_nash(dse_obj_LTU,mu_LTU,false);

    trame::eap_nash(dse_obj_LTU,mu_LTU,tol);
    trame::eap_nash(dse_obj_LTU,mu_LTU,max_iter);

    trame::eap_nash(dse_obj_LTU,mu_LTU,tol,max_iter);
    trame::eap_nash(dse_obj_LTU,mu_LTU,true,tol,max_iter);

    trame::eap_nash(dse_obj_LTU,mu_LTU,u,v);

    trame::eap_nash(dse_obj_LTU,mu_LTU,mux0,mu0y,u,v,&x_first,&tol,&max_iter);
    

    //
    printf("\n*===================    End of eap_nash Test    ===================*\n");
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
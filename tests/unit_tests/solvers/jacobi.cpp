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
 * jacobi solver test
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

    arma::mat alpha(2,3);
    alpha << 1.6 << 3.2 << 1.1 << arma::endr
          << 2.9 << 1.0 << 3.1 << arma::endr;

    arma::mat gamma(2,3);
    gamma << 1.6 << 3.2 << 1.1 << arma::endr
          << 2.9 << 1.0 << 3.1 << arma::endr;

    arma::mat mu_hat(2,3);
    mu_hat << 1.0 << 3.0 << 1.0 << arma::endr
           << 2.0 << 1.0 << 3.0 << arma::endr;

    arma::vec n = 1.2 * arma::sum(mu_hat,1);
    arma::vec m = 1.3 * arma::trans(arma::sum(mu_hat,0));

    int nbX = n.n_elem;
    int nbY = m.n_elem;

    arma::mat phi = alpha + gamma;

    //
    // results

    printf("\n*===================   Start of Jacobi Test   ===================*\n");
    printf("\n");

    //
    // build

    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu> dse_obj_TU;

    trame::arums::logit logit_1(nbX,nbY), logit_2(nbY,nbX);

    dse_obj_TU.build(n,m,phi,logit_1,logit_2,false);

    //

    double tol = 1E-06;
    int max_iter = 5000;

    arma::vec mux0, mu0y;
    arma::mat mu_TU, U, V;

    trame::jacobi(dse_obj_TU,mu_TU);
    trame::jacobi(dse_obj_TU,mu_TU,tol,max_iter);

    trame::jacobi(dse_obj_TU,mu_TU,U,V);

    trame::jacobi(dse_obj_TU,mu_hat,mu_hat,mu_TU,mux0,mu0y,U,V,tol,max_iter);
    trame::jacobi_int(dse_obj_TU,nullptr,&mu_hat,&mu_TU,&mux0,&mu0y,&U,&V,tol,max_iter);


    //
    printf("\n*===================    End of Jacobi Test    ===================*\n");
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

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
 * max_welfare test
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

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat phi = arma::randu(nbX,nbY);

    //
    // results

    printf("\n*===================   Start of max_welfare Test   ===================*\n");
    printf("\n");

    //
    // build

    trame::arums::logit logit_1(nbX,nbY), logit_2(nbY,nbX);

    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu> dse_obj_TU;

    dse_obj_TU.build(n,m,phi,logit_1,logit_2,false);

    //

    double tol = 1E-06;
    int max_iter = 5000;

    arma::vec mux0, mu0y;
    arma::mat mu_TU, U, V;

    trame::max_welfare(dse_obj_TU,mu_TU);
    trame::max_welfare(dse_obj_TU,mu_TU,tol,max_iter);

    arma::cout << "max_welfare solution:\n" << mu_TU << arma::endl;

    //

    double val_out;
    trame::max_welfare(dse_obj_TU,mu_TU,U,V);
    trame::max_welfare(dse_obj_TU,mu_TU,mux0,mu0y,U,V,val_out,tol,max_iter);

    // check template specialization

    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::ltu> dse_obj_LTU;
    trame::max_welfare(dse_obj_LTU,mu_TU,tol,max_iter); // should print a warning and return false

    //
    printf("\n*===================    End of max_welfare Test    ===================*\n");
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

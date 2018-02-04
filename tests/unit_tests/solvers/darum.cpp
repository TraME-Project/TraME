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
 * Darum test
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

    //
    // results

    printf("\n*===================   Start of darum Test   ===================*\n");
    printf("\n");

    //
    // build

    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::ntu> dse_obj_NTU;

    trame::arums::logit logit_1(nbX,nbY), logit_2(nbY,nbX);

    dse_obj_NTU.build(n,m,alpha,gamma,logit_1,logit_2,false);

    // darum

    double tol = 1E-06;
    int max_iter = 5000;

    arma::mat mu_NTU;
    trame::darum(dse_obj_NTU,mu_NTU);
    trame::darum(dse_obj_NTU,mu_NTU,tol,max_iter);

    arma::cout << "darum solution:\n" << mu_NTU << arma::endl;

    //

    arma::vec mu_x0_out, mu_0y_out;
    arma::mat U_out, V_out;

    trame::darum(dse_obj_NTU,mu_NTU,U_out,V_out);
    trame::darum(dse_obj_NTU,mu_NTU,mu_x0_out,mu_0y_out,U_out,V_out,tol,max_iter);

    // check template specialization

    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::ltu> dse_obj_LTU;
    trame::darum(dse_obj_LTU,mu_NTU); // should print a warning and return false

    //
    printf("\n*===================    End of darum Test    ===================*\n");
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

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
    double sigma = 1;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha  = arma::randu(nbX,nbY);
    arma::mat gamma  = arma::randu(nbX,nbY);
    arma::mat lambda = 1 + arma::randu(nbX,nbY);
    arma::mat zeta   = 1 + arma::randu(nbX,nbY);

    arma::mat phi = alpha + gamma;

    arma::mat lambda_LTU = lambda/(lambda+zeta);
    arma::mat phi_LTU = (lambda%alpha + zeta%gamma) / (lambda+zeta);
    //
    // results
    printf("\n*===================   Start of ipfp Test   ===================*\n");
    printf("\n");
    //
    // TU, 
    trame::mfe<trame::mmfs::geo> mfe_obj_TU(sigma,false);
    mfe_obj_TU.build(n,m,phi);

    trame::mfe<trame::mmfs::cd> mfe_obj_LTU(sigma,false);
    mfe_obj_LTU.build(n,m,lambda_LTU,phi_LTU);

    trame::mfe<trame::mmfs::min> mfe_obj_NTU(sigma,false);
    mfe_obj_NTU.build(n,m,alpha,gamma);
    //
    //
    // double tol = 1E-06;
    // int max_iter = 5000;

    arma::mat mu_TU;
    trame::ipfp(mfe_obj_TU,mu_TU);
    // trame::ipfp(mfe_obj_TU,mu_TU,tol);
    // trame::ipfp(mfe_obj_TU,mu_TU,max_iter);
    // trame::ipfp(mfe_obj_TU,mu_TU,tol,max_iter);

    // mfe_obj_TU.solve(mu_TU);

    arma::cout << "Solution of mfe::geo (TU-logit) problem using ipfp:\n" << mu_TU << arma::endl;

    arma::mat mu_NTU;
    trame::ipfp(mfe_obj_NTU,mu_NTU);

    arma::cout << "Solution of mfe::min (NTU-logit) problem using ipfp:\n" << mu_NTU << arma::endl;

    arma::mat mu_LTU;
    trame::ipfp(mfe_obj_LTU,mu_LTU);

    arma::cout << "Solution of mfe::cd (LTU-logit) problem using ipfp:\n" << mu_LTU << arma::endl;
    //
    printf("\n*===================    End of ipfp Test    ===================*\n");
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
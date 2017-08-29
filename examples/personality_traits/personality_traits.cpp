

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

// g++-mp-7 -O3 -Wall -std=c++11 -fopenmp -I./../../include personality_traits.cpp -o personality_traits.test -framework Accelerate -L./../.. -ltrame

#include "trame.hpp"

int main()
{

    //
    // load data

    arma::mat Xvals, Yvals;
    Xvals.load("data/Xvals.txt",arma::auto_detect);
    Yvals.load("data/Yvals.txt",arma::auto_detect);
    
    const int nbX = Xvals.n_rows;
    const int nbY = Yvals.n_rows;

    const int dX = Xvals.n_cols;
    const int dY = Yvals.n_cols;

    arma::mat mu_hat = arma::eye(nbX,nbX);

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    //
    // estimate

    // trame::model< trame::mfe<trame::mmfs::geo> > aff_model;
    trame::model_affinity aff_model;

    aff_model.build(Xvals,Yvals,n,m);

    double lambda = 0.3, val_hat_1, val_hat_2;
    double err_tol = 1E-08;

    arma::mat theta_hat_aff;

    aff_model.mme_regul(mu_hat,theta_hat_aff,lambda,&val_hat_1,&err_tol,nullptr,nullptr,nullptr);
    std::cout << "obj. value with regularization: " << val_hat_1 << std::endl;

    arma::cout << "theta_hat:\n" << arma::reshape(theta_hat_aff,dX,dY) << arma::endl;

    aff_model.mme_woregul(mu_hat,theta_hat_aff,&val_hat_2,&err_tol,nullptr,nullptr,nullptr,nullptr);
    std::cout << "obj. value without regularization: " << val_hat_2 << std::endl;
    
    arma::cout << "theta_hat:\n" << arma::reshape(theta_hat_aff,dX,dY) << arma::endl;

    return 0;
}

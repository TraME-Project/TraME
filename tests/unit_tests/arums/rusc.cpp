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
 * logit class test
 *
 * Keith O'Hara
 * 05/17/2016
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
    arma::mat U(2,3);
    U  << 1.6 << 3.2 << 1.1 << arma::endr
       << 2.9 << 1.0 << 3.1 << arma::endr;
    
    arma::mat mu(2,3);
    mu << 1.0 << 3.0 << 1.0 << arma::endr
       << 2.0 << 1.0 << 3.0 << arma::endr;

    const int nbX = U.n_rows;
    const int nbY = U.n_cols;
    
    arma::vec n = arma::sum(mu,1) + 1.0;

    //
    // results
    printf("\n*===================   Start of arums::rusc test   ===================*\n");
    printf("\n");
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;

    // builds

    arma::mat zeta_temp(1,nbY+1);
    zeta_temp(0,0) = 0.1; zeta_temp(0,1) = 0.2; zeta_temp(0,2) = 0.3; zeta_temp(0,3) = 0.0;

    arma::mat zeta = arma::ones(nbX,1) * zeta_temp;

    arma::cout << "zeta:\n" << zeta << arma::endl;

    trame::arums::rusc arum_obj(zeta,true);
    trame::arums::rusc arum_obj2(nbX,nbY);

    arum_obj.build(zeta,true);
    arum_obj2.build(nbX,nbY);

    // G

    arma::mat mu_sol;
    arum_obj.U = U;
    arum_obj.mu = mu;

    arum_obj.G(n);
    arum_obj.G(n,U,mu_sol);

    // Gstar

    arma::mat U_sol;

    arum_obj.Gstar(n);
    arum_obj.Gstar(n,mu,U_sol);

    // Gbar

    arma::mat mu_bar(2,3);
    mu_bar.fill(2);
    
    arma::mat U_bar_out, mu_bar_out;
    
    arum_obj.Gbar(U,mu_bar,n,U_bar_out,mu_bar_out);

    // D2G

    arma::mat H;

    arum_obj.D2G(n,true);
    arum_obj.D2G(H,n,true);
    arum_obj.D2G(n,U,true);
    arum_obj.D2G(H,n,U,true);

    arum_obj.D2G(H,n,false);

    // D2Gstar

    arum_obj.D2Gstar(n,true);
    arum_obj.D2Gstar(H,n,true);
    arum_obj.D2Gstar(n,U,true);
    arum_obj.D2Gstar(H,n,U,true);

    arum_obj.D2Gstar(H,n,false);

    // dparams

    arma::mat nab_mat;

    arum_obj.dparams_NablaGstar(n,nullptr,true);
    arum_obj.dparams_NablaGstar(nab_mat,n,nullptr,true);
    arum_obj.dparams_NablaGstar(n,mu,nullptr,true);
    arum_obj.dparams_NablaGstar(nab_mat,n,mu,nullptr,true);

    arum_obj.dparams_NablaGstar(nab_mat,n,mu,nullptr,false);

    // simul
    
    trame::arums::empirical logit_sim(nbX,nbY);
    const int sim_seed = 1777, n_draws = 1000;

    logit_sim = arum_obj.simul();
    arum_obj.simul(logit_sim);
    arum_obj.simul(logit_sim, n_draws, sim_seed);

    //
    printf("\n*===================   End of arums::rusc test   ===================*\n");
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

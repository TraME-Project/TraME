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
    
    arma::vec n = arma::sum(mu,1);
    //
    // results
    printf("\n*===================   Start of arums::empirical test   ===================*\n");
    printf("\n");
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;

    // builds
    trame::arums::logit logits(nbX,nbY);

    trame::arums::empirical logit_sim(nbX,nbY);
    const int sim_seed = 1777, n_draws = 1000;

    logits.simul(logit_sim, n_draws, sim_seed);

    // G

    arma::mat mu_sol;
    logit_sim.U = U;
    logit_sim.mu = mu;

    logit_sim.G(n);
    logit_sim.G(n,U,mu_sol);

    // Gstar

    arma::mat U_sol;

    logit_sim.Gstar(n);
    logit_sim.Gstar(n,mu,U_sol);

    // logit_sim2.Gstar(n,mu,U_sol); // no outside option case

    // Gbar

    arma::mat mu_bar(2,3);
    mu_bar.fill(2);
    
    arma::mat U_bar_out, mu_bar_out;
    
    logit_sim.Gbar(U,mu_bar,n,U_bar_out,mu_bar_out);

    // D2G

    arma::mat H;

    logit_sim.D2G(n,true);
    logit_sim.D2G(H,n,true);
    logit_sim.D2G(n,U,true);
    logit_sim.D2G(H,n,U,true);

    logit_sim.D2G(H,n,false);

    // D2Gstar

    logit_sim.D2Gstar(n,true);
    logit_sim.D2Gstar(H,n,true);
    logit_sim.D2Gstar(n,U,true);
    logit_sim.D2Gstar(H,n,U,true);

    logit_sim.D2Gstar(H,n,false);

    // dparams

    arma::mat nab_mat;

    logit_sim.dparams_NablaGstar(n,nullptr,true);
    logit_sim.dparams_NablaGstar(nab_mat,n,nullptr,true);
    logit_sim.dparams_NablaGstar(n,mu,nullptr,true);
    logit_sim.dparams_NablaGstar(nab_mat,n,mu,nullptr,true);

    logit_sim.dparams_NablaGstar(nab_mat,n,mu,&nab_mat,true);
    logit_sim.dparams_NablaGstar(nab_mat,n,mu,nullptr,false);

    // logit_sim2.dparams_NablaGstar(nab_mat,n,mu,nullptr,true);
    // logit_sim2.dparams_NablaGstar(nab_mat,n,mu,nullptr,false);

    arma::mat nab_mat_2; // mat(0,0)
    logit_sim.dparams_NablaGstar(nab_mat,n,mu,&nab_mat_2,true);

    //
    printf("\n*===================   End of arums::empirical test   ===================*\n");
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

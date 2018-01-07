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
 * rsc class test
 *
 * Keith O'Hara
 * 05/17/2016
 *
 * This version:
 * 07/03/2017
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

    int nbX = U.n_rows;
    int nbY = U.n_cols;

    arma::vec n = arma::sum(mu,1) + 1.0;
    //
    // results
    printf("\n*===================   Start of RSC Test   ===================*\n");
    printf("\n");
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;
    //
    // setup
    arma::mat zeta_temp(1,nbY+1);
    zeta_temp(0,0) = 0.1; zeta_temp(0,1) = 0.2; zeta_temp(0,2) = 0.3; zeta_temp(0,3) = 0.0;

    arma::mat zeta = arma::ones(nbX,1) * zeta_temp;

    arma::cout << "zeta: \n" << zeta << arma::endl;
    //
    trame::arums::rsc rsc_obj;
    rsc_obj.build_beta(zeta,2.0,2.0);
    //
    // empirical object:
    int sim_seed = 1777, n_draws = 1000;
    trame::arums::empirical rsc_sim;
    rsc_obj.simul(rsc_sim, n_draws, sim_seed);
    //
    // first compute optimal assignment (mu)
    arma::mat mu_sol, mu_sol_sim;

    double G_val = rsc_obj.G(n,U,mu_sol);
    double G_sim_val = rsc_sim.G(n,U,mu_sol_sim);
    
    std::cout << "G(U) and G-sim(U): \n" << G_val << " and " << G_sim_val << std::endl;
    arma::cout << "\nG -> mu: \n" << mu_sol << "\nG-sim -> mu: \n" << mu_sol_sim << arma::endl;
    //
    // solution to dual problem U*
    arma::mat U_star, U_star_sim;
    
    double Gstar_val = rsc_obj.Gstar(n,mu_sol,U_star);
    double Gstar_sim_val = rsc_sim.Gstar(n,mu_sol,U_star_sim);
    
    std::cout << "G*(mu) and G*-sim(mu): \n" << Gstar_val << " and " << Gstar_sim_val << std::endl;
    arma::cout << "\n\\nabla G*(\\nabla G(U)): \n" << U_star << "\n\\nabla G-sim*(\\nabla G-sim(U)): \n" << U_star_sim << arma::endl;
    //
    // Gbar
    arma::mat mu_bar(2,3);
    mu_bar.fill(2);
    
    arma::mat U_bar_temp, mu_bar_temp;
    arma::mat U_bar_sim_temp, mu_bar_sim_temp;
    
    double val_Gbar = rsc_obj.Gbar(U,mu_bar,n,U_bar_temp,mu_bar_temp);
    double val_Gbar_sim = rsc_sim.Gbar(U,mu_bar,n,U_bar_sim_temp,mu_bar_sim_temp);

    std::cout << "Gbar val: \n" << val_Gbar << std::endl;
    std::cout << "Gbar-sim val: \n" << val_Gbar_sim << std::endl;

    arma::cout << "\nUbar: \n" << U_bar_temp << arma::endl;
    arma::cout << "\nUbar sim: \n" << U_bar_sim_temp << arma::endl;
    arma::cout << "mubar: \n" << mu_bar_temp << arma::endl;
    arma::cout << "mubar sim: \n" << mu_bar_temp << arma::endl;
    //
    // hessian objects
    arma::mat hess;

    rsc_obj.D2Gstar(hess,n,mu,true);

    arma::cout << "\nD2Gstar: \n" << hess << arma::endl;
    //
    arma::mat nablaGstar;
    arma::mat dtheta = arma::eye(rsc_obj.dim_params,rsc_obj.dim_params);

    rsc_obj.dparams_NablaGstar(nablaGstar,n,mu,&dtheta,true);

    arma::cout << "\nparams_NablaGstar: \n" << nablaGstar << arma::endl;
    //
    printf("\n*===================   End of RSC Test   ===================*\n");
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

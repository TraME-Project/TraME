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
 * probit class test
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
    // inputs
    arma::mat U(2,3);
    U  << 1.6 << 3.2 << 1.1 << arma::endr
       << 2.9 << 1.0 << 3.1 << arma::endr;
    
    arma::mat mu(2,3);
    mu << 1.0 << 3.0 << 1.0 << arma::endr
       << 2.0 << 1.0 << 3.0 << arma::endr;

    int nbX = U.n_rows;
    int nbY = U.n_cols;
    
    arma::vec n = arma::sum(mu,1);
    //
    // results
    printf("\n*===================   Start of Probit Test   ===================*\n");
    printf("\n");
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;
    //
    // setup
    double rho = 0.5;

    trame::arums::probit probits(nbX,nbY,rho,true);
    
    probits.unifCorrelCovMatrices();
    arma::cout << "Correlation matrices: \n" << probits.Covar << arma::endl;
    
    int sim_seed = 1777, n_draws = 1000;
    trame::arums::empirical emp_obj;
    
    probits.simul(emp_obj, n_draws, sim_seed);
    //
    // first compute optimal assignment (mu)
    arma::mat mu_sol, mu_sol_sim;

    double G_sim_val = emp_obj.G(n,U,mu_sol_sim);
    
    std::cout << "G-sim(U): \n" << G_sim_val << std::endl;
    arma::cout << "G-sim -> mu: \n" << mu_sol_sim << arma::endl;
    //
    // solution to dual problem U*
    arma::mat U_star_sim;

    double Gstar_sim_val = emp_obj.Gstar(n,mu_sol_sim,U_star_sim);
    
    arma::cout << "G*-sim(mu): \n" << Gstar_sim_val << arma::endl;
    arma::cout << "\\nabla G-sim*(\\nabla G-sim(U)): \n" << U_star_sim << arma::endl;
    //
    // Gbar
    arma::mat mu_bar(2,3);
    mu_bar.fill(2);
    
    arma::mat U_bar_temp, mu_bar_temp;

    double val_Gbar_sim = emp_obj.Gbar(U, mu_bar, n, U_bar_temp, mu_bar_temp);
    
    arma::cout << "Gbar-sim val: \n" << val_Gbar_sim << arma::endl;
    //
    printf("\n*===================   End of Probit Test   ===================*\n");
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

/*
 * logit class test
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

    const int nbX = U.n_rows;
    const int nbY = U.n_cols;
    
    arma::vec n = arma::sum(mu,1);
    //
    // results
    printf("\n*===================   Start of Logit Test   ===================*\n");
    printf("\n");
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;
    //
    // setup
    trame::arums::logit logits(nbX,nbY);

    trame::arums::empirical logit_sim(nbX,nbY);
    const int sim_seed = 1777, n_draws = 1000;
    logits.simul(logit_sim, n_draws, sim_seed);
    //
    // first compute optimal assignment (mu)
    arma::mat mu_sol, mu_sol_sim;

    double G_val = logits.G(n,U,mu_sol);
    double G_sim_val = logit_sim.G(n,U,mu_sol_sim);
    
    std::cout << "G(U) and G-sim(U): \n" << G_val << " and " << G_sim_val << std::endl;
    arma::cout << "\nG -> mu: \n" << mu_sol << "\nG-sim -> mu: \n" << mu_sol_sim << arma::endl;
    arma::cout << "mu normalized:\n" << trame::elem_div(mu_sol,n) << arma::endl;
    //
    // solution to dual problem U*
    arma::mat U_star, U_star_sim;
    
    double Gstar_val = logits.Gstar(n,mu_sol,U_star);
    double Gstar_sim_val = logit_sim.Gstar(n,mu_sol,U_star_sim);
    
    std::cout << "G*(mu) and G*-sim(mu): \n" << Gstar_val << " and " << Gstar_sim_val << std::endl;
    arma::cout << "\n\\nabla G*(\\nabla G(U)): \n" << U_star << "\n\\nabla G-sim*(\\nabla G-sim(U)): \n" << U_star_sim << arma::endl;
    //
    // Gbar
    arma::mat mu_bar(2,3);
    mu_bar.fill(2);
    
    arma::mat U_bar_temp, mu_bar_temp;
    arma::mat U_bar_sim_temp, mu_bar_sim_temp;
    
    double val_Gbar     = logits.Gbar(U,mu_bar,n,U_bar_temp,mu_bar_temp);
    double val_Gbar_sim = logit_sim.Gbar(U,mu_bar,n,U_bar_sim_temp,mu_bar_sim_temp);
    
    std::cout << "Gbar val: \n" << val_Gbar << "\nGbar-sim val: \n" << val_Gbar_sim << std::endl;
    //
    // Hessian tests
    // arma::mat H;
    // arma::mat Hstar;
    
    // logits.D2G(H, n, true);
    // logits.D2Gstar(Hstar, n, true);
    
    // arma::cout << "\nD2G: \n" << H << arma::endl;
    // arma::cout << "D2G*: \n" << Hstar << arma::endl;
    //
    printf("\n*===================   End of Logit Test   ===================*\n");
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
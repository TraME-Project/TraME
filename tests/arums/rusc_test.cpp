/*
 * RUSC class test
 *
 * Keith O'Hara
 * 09/30/2016
 *
 * cd ~/Desktop/SCM/GitHub/TraME/src/tests/arums
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include rusc_test.cpp -c -o rusc_test.o
 * g++-mp-5 -O2 -Wall -o rusc.test rusc_test.o -L/opt/local/lib -ltrame -framework Accelerate
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

    //
    // results
    printf("\n*===================   Start of rusc Test   ===================*\n");
    printf("\n");
    printf("Inputs: \n");
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;
    //
    // setup rusc class object
    int nbX = U.n_rows;
    int nbY = U.n_cols;

    arma::vec n = arma::sum(mu,1) + 1.0;

    arma::mat z_temp(1,nbY+1);
    z_temp(0,0) = 0.1; z_temp(0,1) = 0.2; z_temp(0,2) = 0.3; z_temp(0,3) = 0.0;

    arma::mat zeta = arma::ones(nbX,1) * z_temp;

    arma::cout << "zeta: \n" << zeta << arma::endl;
    //
    // rusc object
    trame::rusc rusc_obj;
    rusc_obj.U = U;
    rusc_obj.mu = mu;

    rusc_obj.build(zeta,true);
    //
    // empirical object:
    int sim_seed = 1777;
    int n_draws = 1000;
    trame::arums::empirical rusc_sim;
    
    rusc_obj.simul(rusc_sim, &n_draws, &sim_seed);
    
    rusc_sim.U = U;
    rusc_sim.mu = mu;
    //
    // first compute optimal assignment (mu)
    double G_val = rusc_obj.G(n);
    double G_sim_val = rusc_sim.G(n);

    std::cout << "G(U) and G-sim(U): \n" << G_val << " and " << G_sim_val << std::endl;

    arma::cout << "\nG -> mu: \n" << rusc_obj.mu_sol << arma::endl;
    arma::cout << "G-sim -> mu: \n" << rusc_sim.mu_sol << arma::endl;
    //
    // solution to dual problem U*
    double Gstar_val = rusc_obj.Gstar(n);
    double Gstar_sim_val = rusc_sim.Gstar(n);

    std::cout << "G*(mu) and G*-sim(mu): \n" << Gstar_val << " and " << Gstar_sim_val << std::endl;

    arma::cout << "\n\\nabla G*(\\nabla G(U)): \n" << rusc_obj.U_sol << arma::endl;
    arma::cout << "\\nabla G-sim*(\\nabla G-sim(U)): \n" << rusc_sim.U_sol << arma::endl;
    //
    // Gbar
    arma::mat mu_bar(2,3);
    mu_bar.fill(2);
    
    arma::mat U_bar_temp, mu_bar_temp;
    arma::mat U_bar_sim_temp, mu_bar_sim_temp;
    
    double val_Gbar = rusc_obj.Gbar(U,mu_bar,n,U_bar_temp,mu_bar_temp);
    double val_Gbar_sim = rusc_sim.Gbar(U,mu_bar,n,U_bar_sim_temp,mu_bar_sim_temp);

    std::cout << "Gbar val: \n" << val_Gbar << std::endl;
    std::cout << "Gbar-sim val: \n" << val_Gbar_sim << std::endl;

    arma::cout << "\nUbar: \n" << U_bar_temp << arma::endl;
    arma::cout << "\nUbar sim: \n" << U_bar_sim_temp << arma::endl;
    arma::cout << "mubar: \n" << mu_bar_temp << arma::endl;
    arma::cout << "mubar sim: \n" << mu_bar_temp << arma::endl;
    //
    // hessian objects
    /*arma::mat hess;

    rusc_obj.D2Gstar(hess,n,true);

    arma::cout << "\nD2Gstar: \n" << hess << arma::endl;
    //
    arma::mat nablaGstar;
    arma::mat dtheta = arma::eye(rusc_obj.nbParams,rusc_obj.nbParams);

    rusc_obj.dtheta_NablaGstar(nablaGstar,n,&dtheta,true);

    arma::cout << "\ndtheta_NablaGstar: \n" << nablaGstar << arma::endl;*/
    //
    printf("\n*===================   End of rusc Test   ===================*\n");
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

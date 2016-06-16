/*
 * RSC class test
 *
 * Keith O'Hara
 * 05/17/2016
 *
 * cd ~/Desktop/SCM/GitHub/TraME/src/tests
 *
 * g++-mp-5 -Wall -O2 -std=c++11 -I/opt/local/include -I/usr/local/include -I/Library/gurobi650/mac64/include arums_rsc_test.cpp -c -o arums_rsc_test.o
 * gcc-mp-5 -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include ../lp/generic_lp.c -c -o ../lp/generic_lp.o
 * gfortran-mp-5 -O2 ../prob/prob.f90  -c -o ../prob/prob.o
 * gfortran-mp-5 -O2 ../math/quadpack_double.f90  -c -o ../math/quadpack_double.o
 * gfortran-mp-5 -O2 ../prob/aux.f90  -c -o ../prob/aux.o
 * g++-mp-5 -o arums_rsc_test.test ../lp/generic_lp.o ../prob/prob.o ../math/quadpack_double.o ../prob/aux.o arums_rsc_test.o -L/Library/gurobi650/mac64/lib -L/usr/local/lib -lgurobi65 -lnlopt -lgfortran -framework Accelerate
 */

#ifndef __clang__
#define TRAME_USE_GUROBI_C
#endif

#include "armadillo"

#include "../headers/trame_structs.hpp"
#include "../headers/trame_aux.hpp"

#include "../nlopt/generic_nlopt.hpp"

#include "../headers/arums_empirical.hpp"
#include "../headers/arums_rsc.hpp"

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
    printf("\n*===================   Start of RSC Test   ===================*\n");
    printf("\n");
    printf("Inputs: \n");
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;
    //
    // setup RSC class object
    int nbX = U.n_rows;
    int nbY = U.n_cols;

    arma::vec n = arma::sum(mu,1) + 1.0;

    arma::mat z_temp(1,nbY+1);
    z_temp(0,0) = 0.1; z_temp(0,1) = 0.2; z_temp(0,2) = 0.3; z_temp(0,3) = 0.0;

    arma::mat zeta = arma::ones(nbX,1) * z_temp;

    arma::cout << "zeta: \n" << zeta << arma::endl;
    //
    // RSC object
    RSC rsc_obj;
    rsc_obj.U = U;
    rsc_obj.mu = mu;

    rsc_obj.build_beta(zeta,2.0,2.0);
    //
    // empirical object:
    int n_draws = 1000;
    empirical rsc_sim;
    
    rsc_obj.simul(rsc_sim, n_draws, (int) 1777);
    
    rsc_sim.U = U;
    rsc_sim.mu = mu;
    //
    // first compute optimal assignment (mu)
    double G_val = rsc_obj.G(n);
    double G_sim_val = rsc_sim.G(n);

    std::cout << "G(U) and G-sim(U): \n" << G_val << " and " << G_sim_val << std::endl;

    arma::cout << "\nG -> mu: \n" << rsc_obj.mu_sol << arma::endl;
    arma::cout << "G-sim -> mu: \n" << rsc_sim.mu_sol << arma::endl;
    //
    // solution to dual problem U*
    double Gstar_val = rsc_obj.Gstar(n);
    double Gstar_sim_val = rsc_sim.Gstar(n);

    std::cout << "G*(mu) and G*-sim(mu): \n" << Gstar_val << " and " << Gstar_sim_val << std::endl;

    arma::cout << "\n\\nabla G*(\\nabla G(U)): \n" << rsc_obj.U_sol << arma::endl;
    arma::cout << "\\nabla G-sim*(\\nabla G-sim(U)): \n" << rsc_sim.U_sol << arma::endl;
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
    arma::cout << "mubar: \n" << mu_bar_temp << arma::endl;
    //
    // hessian objects
    arma::mat hess;

    rsc_obj.D2Gstar(hess,n,true);

    arma::cout << "\nD2Gstar: \n" << hess << arma::endl;
    //
    arma::mat nablaGstar;
    arma::mat dtheta = arma::eye(rsc_obj.nbParams,rsc_obj.nbParams);

    rsc_obj.dtheta_NablaGstar(nablaGstar,n,&dtheta,true);

    arma::cout << "\ndtheta_NablaGstar: \n" << nablaGstar << arma::endl;
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

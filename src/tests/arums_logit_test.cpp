/*
 * logit class test
 *
 * Keith O'Hara
 * 05/17/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/tests
 * clang++ -O2 -Wall -std=c++11 -I/opt/local/include -I/Library/gurobi650/mac64/include arums_logit_test.cpp -o arums_logit.test -L/Library/gurobi650/mac64/lib -lgurobi_c++ -lgurobi65 -framework Accelerate
 * clang++ -O2 -Wall -static-libgcc -I/opt/local/include -I/Library/gurobi650/mac64/include arums_logit_test.cpp -o arums_logit.test -L/Library/gurobi650/mac64/lib -lgurobi65 -framework Accelerate
 *
 * gcc-mp-5 -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include ../lp/generic_lp.c -c -o ../lp/generic_lp.o
 * g++-mp-5 -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include arums_logit_test.cpp -c -o arums_logit_test.o
 * g++-mp-5 -O2 -Wall -o arums_logit.test ../lp/generic_lp.o arums_logit_test.o -L/Library/gurobi650/mac64/lib -lgurobi65 -framework Accelerate
 */

#define TRAME_USE_GUROBI_C

#include "armadillo"

#include "../headers/trame_structs.hpp"

#include "../headers/arums_empirical.hpp"
#include "../headers/arums_logit.hpp"

int main()
{
    arma::mat U(2,3);
    U  << 1.6 << 3.2 << 1.1 << arma::endr
       << 2.9 << 1.0 << 3.1 << arma::endr;
    
    arma::mat mu(2,3);
    mu << 1.0 << 3.0 << 1.0 << arma::endr
       << 2.0 << 1.0 << 3.0 << arma::endr;
    
    arma::cout << "\nU: \n" << U << arma::endl;
    arma::cout << "mu: \n" << mu << arma::endl;
    //
    int nbX = U.n_rows;
    int nbY = U.n_cols;
    
    arma::vec n = arma::sum(mu,1);
    //
    logit logits;
    logits.nbX = nbX;
    logits.nbY = nbY;
    logits.nbParams = 1;
    logits.sigma = 1.0;
    logits.outsideOption = true;
    
    logits.U = U;
    logits.mu = mu;
    //
    double G_val = logits.G(n);
    
    arma::cout << "G -> mu: \n" << logits.mu_sol << arma::endl;
    arma::cout << "G(U): \n" << G_val << arma::endl;
    // 
    arma::mat U_star;
    double Gstar_val = logits.Gstar(U_star, n);
    
    arma::cout << "\n\\nabla G*(\\nabla G(U)): \n" << U_star << arma::endl;
    arma::cout << "G*(mu): \n" << Gstar_val << arma::endl;
    //
    arma::mat H;
    logits.D2G(H, n, true);
    
    arma::cout << "\nD2G: \n" << H << arma::endl;
    //
    arma::mat Hstar;
    logits.D2Gstar(Hstar, n, true);
    
    arma::cout << "D2G*: \n" << Hstar << arma::endl;
    /*
     * Gbar test
     */
    arma::mat mubar(2,3);
    mubar.fill(2);

    arma::mat Ubar_temp, mubar_temp;

    double valGbar = logits.Gbar(U_star,mubar,n,Ubar_temp,mubar_temp);
    
    std::cout << "Gbar val: \n" << valGbar << std::endl;
    /*
     * Simulated logit
     */
    int n_draws = 1000;
    empirical emp_obj;
    
    logits.simul(emp_obj, n_draws, (int) 1777);
    emp_obj.U = U;
    emp_obj.mu = mu;

    emp_obj.G(n);

    arma::cout << "\nG emp -> mu: \n" << emp_obj.mu_sol << arma::endl;
    
    arma::mat Ustar_emp;
    emp_obj.Gstar(n,Ustar_emp);
    
    arma::cout << "\\nabla G*(\\nabla G(U emp)): \n" << Ustar_emp << arma::endl;
    //
    //emp_obj.Gbarx(U.row(0).t(),(mubar.row(0).t())/n(0),Ubar_temp,mubar_temp,0);
    double valGbarEmp = emp_obj.Gbar(Ustar_emp,mubar,n,Ubar_temp,mubar_temp);
    std::cout << "Gbar emp val: \n" << valGbarEmp << std::endl;
    //
    return 0;
}
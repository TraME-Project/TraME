/*
 * cd ~/Desktop/"Google Drive"/GitHub/TraME/src/tests
 * clang++ -O2 -Wall -std=c++11 -I/opt/local/include -I/Library/gurobi650/mac64/include arums_logit_test.cpp -o arums_logit.test -L/Library/gurobi650/mac64/lib -lgurobi_c++ -lgurobi65 -framework Accelerate
 */

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
    //
    double G_val;
    logits.G(G_val,n);
    
    arma::cout << "G -> mu: \n" << logits.mu << arma::endl;
    //
    double Gstar_val;    
    arma::mat U_star;
    logits.Gstar(Gstar_val, U_star, n);
    
    arma::cout << "\\nabla G*(\\nabla G(U)): \n" << U_star << arma::endl;
    //
    arma::mat H;
    logits.D2G(H, n, (int) 1);
    
    arma::cout << "D2G: \n" << H << arma::endl;
    //
    arma::mat Hstar;
    logits.D2Gstar(Hstar, n, (int) 1);
    
    arma::cout << "D2G*: \n" << Hstar << arma::endl;
    //
    int n_draws = 1000;
    empirical emp_obj;
    
    logits.simul(emp_obj,n_draws, (int) 777);
    emp_obj.U = U;

    emp_obj.G(n);

    arma::cout << "G emp -> mu: \n" << emp_obj.mu << arma::endl;
    
    arma::mat Ustar_emp;
    emp_obj.Gstar(n,Ustar_emp);
    
    arma::cout << "\\nabla G*(\\nabla G(U emp)): \n" << Ustar_emp << arma::endl;
    /*
     * Gbar test
     */
    arma::mat mubar(2,3);
    mubar.fill(2);

    arma::mat Ubar_temp, mubar_temp;
    //logits.Gbarx(U.row(0).t(),(mubar.row(0).t())/n(0),Ubar_temp,mubar_temp);
    double valGbar = logits.Gbar(U,mubar,n,Ubar_temp,mubar_temp);
    std::cout << "Gbar val: \n" << valGbar << std::endl;
    //
    //emp_obj.Gbarx(U.row(0).t(),(mubar.row(0).t())/n(0),Ubar_temp,mubar_temp,0);
    double valGbarEmp = emp_obj.Gbar(U,mubar,n,Ubar_temp,mubar_temp);
    std::cout << "Gbar emp val: \n" << valGbarEmp << std::endl;
    //
    return 0;
}
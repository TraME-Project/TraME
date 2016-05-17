/*
 * cd ~/Desktop/SCM/GitHub/TraME/src/tests
 * clang++ -O2 -Wall -std=c++11 -I/opt/local/include -I/Library/gurobi650/mac64/include arums_none_test.cpp -o arums_none.test -L/Library/gurobi650/mac64/lib -lgurobi_c++ -lgurobi65 -framework Accelerate
 */
 
#include "armadillo"

#include "../headers/trame_structs.hpp"

#include "../headers/arums_empirical.hpp"
#include "../headers/arums_none.hpp"

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
    none none_obj;
    none_obj.nbX = nbX;
    none_obj.nbY = nbY;
    none_obj.nbParams = 0;
    
    none_obj.U = U;
    none_obj.mu = mu;
    //
    double valx_Gx1 = none_obj.Gx(U.row(0).t());
    std::cout << "G val: \n" << valx_Gx1 << std::endl;
    //
    arma::mat mubar(2,3);
    mubar.fill(2);

    arma::mat Ubar_temp, mubar_temp;
    double valGbar = none_obj.Gbarx(U.row(0).t(),mubar.row(0).t(),Ubar_temp,mubar_temp);
    std::cout << "Gbar val: \n" << valGbar << std::endl;
    //
    return 0;
}
/*
 * RSC class test
 *
 * Keith O'Hara
 * 05/17/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/tests
 * g++-mp-5 -Wall -O2 -std=c++11 -fopenmp -I/opt/local/include -I/Library/gurobi650/mac64/include arums_rsc_test.cpp -c -o arums_rsc_test.o
 * gfortran-mp-5 -O2 ../prob/prob.f90  -c -o ../prob/prob.o
 * gfortran-mp-5 -O2 ../math/quadpack_double.f90  -c -o ../math/quadpack_double.o
 * gfortran-mp-5 -O2 ../prob/aux.f90  -c -o ../prob/aux.o
 * g++-mp-5 -o arums_rsc_test.test ../prob/prob.o ../math/quadpack_double.o ../prob/aux.o arums_rsc_test.o -L/Library/gurobi650/mac64/lib -lgurobi_c++ -lgurobi65 -lgfortran -fopenmp -framework Accelerate
 */

#include "armadillo"

#include "../headers/trame_structs.hpp"
#include "../headers/trame_aux.hpp"

#include "../headers/arums_empirical.hpp"
#include "../headers/arums_rsc.hpp"

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
    
    arma::vec n = arma::sum(mu,1) + 1.0;
    
    arma::mat z_temp(1,nbY+1);
    z_temp(0,0) = 0.1; z_temp(0,1) = 0.2; z_temp(0,2) = 0.3; z_temp(0,3) = 0.0;
    
    arma::mat zeta = arma::ones(nbX,1) * z_temp;
    
    arma::cout << "zeta: \n" << zeta << arma::endl;
    //
    RSC rsc_obj;
    rsc_obj.U = U;
    rsc_obj.mu = mu;
    
    rsc_obj.build_beta(zeta,2.0,2.0);
    //
    arma::cout << rsc_obj.aux_Influence_lhs << arma::endl;
    arma::cout << rsc_obj.aux_Influence_rhs << arma::endl;
    arma::cout << rsc_obj.aux_DinvPsigma << arma::endl;
    arma::cout << rsc_obj.aux_Psigma << arma::endl;
    //
    std::cout << rsc_obj.cdf(0.5) << std::endl;
    std::cout << rsc_obj.pdf(0.5) << std::endl;
    std::cout << rsc_obj.quantile(0.5) << std::endl;
    std::cout << rsc_obj.pot(0.5) << std::endl;
    //
    arma::vec mu_x;
    
    rsc_obj.Gx(mu_x,(int) 0);
    
    double G_val;
    G_val = rsc_obj.G(n);
    
    arma::cout << rsc_obj.mu_sol << arma::endl;
    std::cout << G_val << std::endl;
    //
    arma::vec U_x;
    
    rsc_obj.Gstarx(U_x,n(0,0),(int) 0);
    
    double Gstar_val;
    Gstar_val = rsc_obj.Gstar(n);
    
    arma::cout << rsc_obj.U_sol << arma::endl;
    std::cout << Gstar_val << std::endl;
    //
    arma::mat hess;
    
    rsc_obj.D2Gstar(hess,n,true);
    
    arma::cout << hess << arma::endl;
    //
    arma::mat nablaGstar;
    arma::mat dtheta = arma::eye(rsc_obj.nbParams,rsc_obj.nbParams);
    
    rsc_obj.dtheta_NablaGstar(nablaGstar,n,&dtheta,true);
    
    arma::cout << nablaGstar << arma::endl;
    //
    return 0;
}

/*
 * sparse matrix test
 *
 * Keith O'Hara
 * 05/17/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/tests
 * clang++ -O2 -Wall -std=c++11 -I/opt/local/include sp_mat_test.cpp -o sp_mat.test -framework Accelerate
 */

#include "armadillo"

#include "../headers/sparse_kron.hpp"

int main()
{
    arma::sp_mat A(5,6);
    arma::sp_mat B(6,5);

    A(0,0) = 1;
    A(1,0) = 2;

    B(0,0) = 3;
    B(0,1) = 4;

    arma::sp_mat C = 2*B;
    arma::sp_mat D = A*C;

    //
    arma::sp_mat kron_AB = kron_sp(A,B);
    
    arma::cout << kron_AB << arma::endl;
    //
    arma::mat F = arma::zeros(5,6);
    F(0,0) = 1;
    F(1,0) = 2;
    
    arma::sp_mat kron_FB = kron_sp(F,B);
    
    arma::cout << kron_FB << arma::endl;
    //
    arma::mat G = arma::zeros(6,5);
    G(0,0) = 3;
    G(0,1) = 4;
    
    arma::sp_mat kron_AG = kron_sp(A,G);
    
    arma::cout << kron_AG << arma::endl;
    //

    // batch insertion of two values at (5, 6) and (9, 9)
    arma::umat locations;
    locations << 5 << 9 << arma::endr
              << 6 << 9 << arma::endr;

    arma::vec values;
    values << 1.5 << 3.2 << arma::endr;

    arma::sp_mat X(locations, values);
    
    //
    int nbY = 3;
    int nbDraws = 10000;
    
    arma::sp_mat A1(nbY,nbDraws);
    arma::sp_mat A21 = arma::speye(nbY,nbY);
    
    //arma::mat A22 = arma::ones(nbDraws,1);
    arma::mat A22 = arma::ones(nbY,1);
    arma::mat A2 = A21 * A22;
    
    arma::cout << A2 << arma::endl;
    
    return 0;
}
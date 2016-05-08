/*
 * test which_max function
 *
 * cd ~/Desktop/"Google Drive"/GitHub/TraME/src/tests
 * clang++ -O2 -Wall -I/opt/local/include which_max.cpp -o which_max.test -framework Accelerate
 */

#include "armadillo"

#include "../headers/which_max.hpp"

int main()
{
    arma::mat A(3,3);
    A  << 3.0 << 4.0 << 2.0 << arma::endr
       << 2.0 << 1.0 << 2.0 << arma::endr
       << 1.0 << 3.0 << 4.0 << arma::endr;
    
    arma::uvec max_vec_1 = which_max(&A, 0);
    arma::uvec max_vec_2 = which_max(&A, 1);
    
    arma::cout << max_vec_1 << arma::endl;
    arma::cout << max_vec_2 << arma::endl;
    
    return 0;
}
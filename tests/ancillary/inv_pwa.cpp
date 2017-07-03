/*
 * inv_pwa test
 *
 * Keith O'Hara
 * 08/01/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/aux
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include inv_pwa.cpp -c -o inv_pwa_test.o
 * g++-mp-7 -O2 -Wall inv_pwa.cpp -o inv_pwa.test -L/opt/local/lib -ltrame -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    arma::vec a(1);
    a(0) = 2;
    arma::mat B = arma::ones(1,1) + 1;
    arma::mat C = arma::ones(1,1);
    double k = 1;
    // this should return a value of 1
    arma::vec res = trame::inv_pwa(a,B,C,k);
    //
    return 0;
}

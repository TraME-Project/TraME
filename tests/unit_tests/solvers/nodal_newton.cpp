/*
 * nodal_newton test
 *
 * Keith O'Hara
 * 02/07/2017
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/equilibrium
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include nodal_newton_test.cpp -o nodal_newton.test -L/opt/local/lib -ltrame -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    //
    // inputs:

    int nbX = 10;
    int nbY = 8;
    double sigma = 1;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat phi = 1.0 + arma::randu(nbX,nbY);
    
    //
    // results

    printf("\n*===================   Start of nodal_newton Test   ===================*\n");
    printf("\n");

    //
    // build
    
    trame::mfe<trame::mmfs::geo> mfe_obj_TU(sigma,false);
    mfe_obj_TU.build(n,m,phi);

    //

    double tol = 1E-06;
    int max_iter = 5000;

    arma::mat mu_TU;
    trame::nodal_newton(mfe_obj_TU,mu_TU);

    trame::nodal_newton(mfe_obj_TU,mu_TU,tol);
    trame::nodal_newton(mfe_obj_TU,mu_TU,max_iter);
    trame::nodal_newton(mfe_obj_TU,mu_TU,tol,max_iter);

    arma::mat U_out, V_out;
    trame::nodal_newton(mfe_obj_TU,mu_TU,U_out,V_out);

    double val_out;
    arma::vec mu_x0_out, mu_0y_out;
    trame::nodal_newton(mfe_obj_TU,mu_TU,mu_x0_out,mu_0y_out,U_out,V_out,val_out,&tol,&max_iter);
    //
    printf("\n*===================    End of nodal_newton Test    ===================*\n");
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
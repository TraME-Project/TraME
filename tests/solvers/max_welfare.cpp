/*
 * max_welfare test
 *
 * Keith O'Hara
 * 24/10/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/equilibrium
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include max_welfare_test.cpp -o max_welfare.test -L/opt/local/lib -ltrame -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs:
    int nbX = 5;
    int nbY = 3;
    double sigma = 1;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat phi = arma::randu(nbX,nbY);
    //
    // results
    printf("\n*===================   Start of max_welfare Test   ===================*\n");
    printf("\n");
    //
    // TU
    trame::arums::logit logit_1(nbX,nbY), logit_2(nbY,nbX);
    
    trame::mfe<trame::mmfs::geo> mfe_obj_TU(sigma,false);
    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu> dse_obj_TU;

    mfe_obj_TU.build(n,m,phi);
    dse_obj_TU.build(n,m,phi,logit_1,logit_2,false);
    //
    //
    //double tol = 1E-06;
    //int max_iter = 5000;

    arma::mat mu_TU_1, mu_TU_2;
    trame::ipfp(mfe_obj_TU,mu_TU_1);

    //trame::max_welfare(dse_obj_TU,mu_TU_2);
    dse_obj_TU.solve(mu_TU_2, (char*) "maxWelfare");

    arma::cout << "Solution of TU-logit problem using ipfp:\n" << mu_TU_1 << arma::endl;
    arma::cout << "Solution of TU-logit problem using max_welfare:\n" << mu_TU_2 << arma::endl;
    //
    // try RSC problem
    /*arma::mat zetaG = arma::ones(nbX,1) * arma::randu(1,nbY+1);
    arma::mat zetaH = arma::ones(nbY,1) * arma::randu(1,nbX+1);

    trame::rsc rsc_1, rsc_2;
    rsc_1.build_beta(zetaG,2.0,2.0);
    rsc_2.build_beta(zetaH,2.0,2.0);
    
    trame::dse<trame::rsc> dse_obj_TU_2;

    dse_obj_TU_2.build_TU(n,m,phi,rsc_1,rsc_2,false);

    arma::mat mu_rsc;
    trame::max_welfare(dse_obj_TU_2,mu_rsc);*/
    //
    printf("\n*===================    End of max_welfare Test    ===================*\n");
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
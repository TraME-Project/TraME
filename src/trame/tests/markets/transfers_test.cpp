/*
 * transfers class test
 *
 * Keith O'Hara
 * 08/01/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/markets
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include -I./../../headers -I/usr/local/include transfers_test.cpp -c -o transfers_test.o
 * g++-mp-5 -O2 -Wall -o transfers.test transfers_test.o -L/opt/local/lib -ltrame -framework Accelerate
 */

#include "trame.hpp"

#include "../../headers/markets/transfers.hpp"
#include "../../headers/markets/MMF.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs:
    int nbX = 18;
    int nbY = 5;
    double sigma = 1;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha  = arma::randu(nbX,nbY);
    arma::mat gamma  = arma::randu(nbX,nbY);
    arma::mat lambda = 1 + arma::randu(nbX,nbY);
    arma::mat zeta   = 1 + arma::randu(nbX,nbY);

    arma::mat phi = alpha + gamma;
    // 
    // build_market_TU_logit
    transfers transf_obj_TU;
    transf_obj_TU.build_TU(phi);

    logit logitM;
    logitM.build(nbX,nbY,sigma,true);

    logit logitW;
    logitM.build(nbX,nbY,sigma,true);

    arma::mat K = arma::exp(phi/(sigma*2));
    MMF mmf_obj_TU;
    mmf_obj_TU.build_TU(n,m,K,false);

    arma::vec x_inv_TU = mmf_obj_TU.marg_x_inv(NULL,m);
    arma::vec y_inv_TU = mmf_obj_TU.marg_y_inv(NULL,x_inv_TU);

    arma::cout << x_inv_TU << arma::endl;
    arma::cout << y_inv_TU << arma::endl;
    //
    // build_market_NTU_logit
    transfers transf_obj_NTU;
    transf_obj_NTU.build_NTU(alpha,gamma);

    arma::mat A = arma::exp(alpha/sigma);
    arma::mat B = arma::exp(gamma/sigma);
    MMF mmf_obj_NTU;
    mmf_obj_NTU.build_NTU(n,m,A,B,false);

    arma::vec x_inv_NTU = mmf_obj_NTU.marg_x_inv(NULL,m);
    arma::vec y_inv_NTU = mmf_obj_NTU.marg_y_inv(NULL,x_inv_NTU);

    arma::cout << x_inv_NTU << arma::endl;
    arma::cout << y_inv_NTU << arma::endl;
    //
    // build_market_LTU_logit
    arma::mat lambda_LTU = lambda/(lambda+zeta);
    arma::mat phi_LTU = (lambda%alpha + zeta%gamma) / (lambda+zeta);

    transfers transf_obj_LTU;
    transf_obj_LTU.build_LTU(lambda_LTU,phi_LTU);

    K = arma::exp(phi_LTU/sigma);
    MMF mmf_obj_LTU;
    mmf_obj_LTU.build_LTU(n,m,lambda_LTU,K,false);

    arma::vec x_inv_LTU = mmf_obj_LTU.marg_x_inv(NULL,m);
    arma::vec y_inv_LTU = mmf_obj_LTU.marg_y_inv(NULL,x_inv_LTU);

    arma::cout << x_inv_LTU << arma::endl;
    arma::cout << y_inv_LTU << arma::endl;
    //
    // results
    printf("\n*===================   Start of MMF Test   ===================*\n");
    printf("\n");
    
    //
    printf("\n*===================    End of MMF Test    ===================*\n");
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
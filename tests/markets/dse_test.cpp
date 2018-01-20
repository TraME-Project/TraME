/*
 * DSE class test
 *
 * Keith O'Hara
 * 08/01/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/trame/tests/markets
 *
 * g++-mp-5 -O2 -Wall -std=c++11 -I/opt/local/include/trame dse_test.cpp -o dse.test -L/opt/local/lib -ltrame -framework Accelerate
 * g++-mp-5 -O2 -Wall -std=c++11 -I./../../include dse_test.cpp -o dse.test -L./../../ -ltrame -framework Accelerate
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    //
    // inputs:
    int nbX = 8;
    int nbY = 5;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha  = arma::randu(nbX,nbY);
    arma::mat gamma  = arma::randu(nbX,nbY);

    arma::mat phi = alpha + gamma;
    //
    // results
    printf("\n*===================   Start of DSE Test   ===================*\n");
    printf("\n");
    //
    // TU
    trame::dse<trame::arums::none,trame::arums::none,trame::transfers::tu> dse_obj_TU;
    dse_obj_TU.build(n,m,phi,false);
    //
    double val;
    arma::vec mux0, mu0y, u, v;
    arma::mat mu_TU, residuals;
    trame::oap_lp(dse_obj_TU, mu_TU, mux0, mu0y, u, v, true, val, residuals);

    std::cout << "Solution of TU-none problem using oap_lp:\n" << std::endl;
    arma::cout << "u:\n" << u << "\n v:\n" << v << "\n resid:\n" << residuals << arma::endl;
    //
    // NTU
    arma::mat alpha_d(2,3);
    alpha_d << 1.6 << 3.2 << 1.1 << arma::endr 
            << 2.9 << 1.0 << 3.1 << arma::endr;

    arma::mat gamma_d(2,3);
    gamma_d << 1.6 << 3.2 << 1.1 << arma::endr 
            << 2.9 << 1.0 << 3.1 << arma::endr;

    arma::mat mu_hat(2,3);
    mu_hat << 1.0 << 3.0 << 1.0 << arma::endr 
           << 2.0 << 1.0 << 3.0 << arma::endr;

    arma::vec n_d = 1.2 * arma::sum(mu_hat,1);
    arma::vec m_d = 1.3 * arma::trans(arma::sum(mu_hat,0));

    int nbX_d = n_d.n_elem;
    int nbY_d = m_d.n_elem;
    //
    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::ntu> dse_obj_NTU;

    trame::arums::logit logit_d_1, logit_d_2;
    logit_d_1.build(nbX_d,nbY_d,1.0,true);
    logit_d_2.build(nbY_d,nbX_d,1.0,true);

    dse_obj_NTU.build(n_d,m_d,alpha_d,gamma_d,logit_d_1,logit_d_2,false);
    //
    arma::mat mu_NTU, U, V;
    // trame::darum(dse_obj_NTU, true, nullptr, mu_NTU, mux0, mu0y, U, V);
    // trame::darum(dse_obj_NTU, mu_NTU);
    dse_obj_NTU.solve(mu_NTU);
    // dse_obj_NTU.solve(mu_NTU,(char*) "darum");

    std::cout << "Solution of NTU-logit problem using DA-RUM:\n" << mu_NTU << std::endl;

    // //
    trame::dse<trame::arums::logit,trame::arums::logit,trame::transfers::tu> dse_obj_TU_2;

    dse_obj_TU_2.build(n_d,m_d,alpha_d + gamma_d,logit_d_1,logit_d_2,false);

    // dse_obj_TU_2.solve(mu_TU);
    dse_obj_TU_2.solve(mu_TU,(char*) "m");
    //
    printf("\n*===================    End of DSE Test    ===================*\n");
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
/*################################################################################
  ##
  ##   Copyright (C) 2015-2018 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
  ##
  ##   This file is part of TraME.
  ##
  ##   TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * transfers::etu class unit test
 *
 * Keith O'Hara
 * 02/07/2017
 *
 * This version:
 * 08/17/2017
 */

#include "trame.hpp"

int main()
{
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();

    //
    // setup:

    int nbX = 10;
    int nbY = 8;
    
    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha = arma::randu(nbX,nbY);
    arma::mat gamma = arma::randu(nbX,nbY);
    arma::mat tau   = arma::randu(nbX,nbY);

    arma::mat U = arma::randu(nbX,nbY);
    arma::mat V = arma::randu(nbX,nbY);
    //
    // results
    printf("\n*===================   Start of transfers::etu Test   ===================*\n");
    printf("\n");
    //

    trame::transfers::etu transfers_obj;
    
    transfers_obj.build(alpha,gamma,tau,false);

    arma::uvec xs(1); xs(0) = 1;
    arma::uvec ys(1); ys(0) = 2;

    // Psi

    transfers_obj.Psi(U,V);

    transfers_obj.Psi(U,V,nullptr,nullptr);
    transfers_obj.Psi(U(xs,ys),V(xs,ys),&xs,&ys);

    double U_xy = U(1,2);
    double V_xy = V(1,2);
    transfers_obj.Psi(U_xy,V(xs,ys),&xs,&ys);
    transfers_obj.Psi(U(xs,ys),V_xy,&xs,&ys);
    transfers_obj.Psi(U_xy,V_xy,xs(0),ys(0));

    // du

    transfers_obj.du_Psi(U,V);

    transfers_obj.du_Psi(U,V,nullptr,nullptr);
    transfers_obj.du_Psi(U(xs,ys),V(xs,ys),&xs,&ys);

    // dparams_Psi

    transfers_obj.dparams_Psi(U,V);
    arma::mat dpars = arma::vectorise(arma::join_cols(arma::join_cols(alpha,gamma),tau));
    transfers_obj.dparams_Psi(U,V,&dpars);

    // Ucal and Vcal

    transfers_obj.Ucal(U);
    transfers_obj.Ucal(U(xs,ys),&xs,&ys);
    transfers_obj.Ucal(U_xy,1,2);

    transfers_obj.Vcal(U);
    transfers_obj.Vcal(U(xs,ys),&xs,&ys);
    transfers_obj.Vcal(U_xy,1,2);

    // UW and VW

    transfers_obj.UW(U);
    transfers_obj.UW(U(xs,ys),&xs,&ys);
    transfers_obj.UW(U_xy,1,2);

    transfers_obj.VW(U);
    transfers_obj.VW(U(xs,ys),&xs,&ys);
    transfers_obj.VW(U_xy,1,2);

    // dw

    transfers_obj.dw_UW(U);
    transfers_obj.dw_UW(U(xs,ys),&xs,&ys);

    transfers_obj.dw_VW(U);
    transfers_obj.dw_VW(U(xs,ys),&xs,&ys);

    // WU and WV

    transfers_obj.WU(U);
    transfers_obj.WU(U(xs,ys),&xs,&ys);

    transfers_obj.WV(U);
    transfers_obj.WV(U(xs,ys),&xs,&ys);

    // generate MMF object

    trame::mmfs::ces mmf_obj = transfers_obj.gen_mmf();
    transfers_obj.gen_mmf(mmf_obj);

    transfers_obj.trans();

    //
    printf("\n*===================    End of transfers::etu Test    ===================*\n");
    printf("\n");
    //
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
        
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);
        
    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";
    //
    return 0;
}

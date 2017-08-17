/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
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
 * mmfs::min class unit test
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
    // setip:

    int nbX = 10;
    int nbY = 8;

    arma::vec n = arma::ones(nbX,1);
    arma::vec m = arma::ones(nbY,1);

    arma::mat alpha = arma::randu(nbX,nbY);
    arma::mat gamma = arma::randu(nbX,nbY);
    //
    // results
    printf("\n*===================   Start of mmfs::min Test   ===================*\n");
    printf("\n");
    //

    trame::mmfs::min mmf_obj;

    mmf_obj.build(alpha,gamma,false);

    arma::uvec xs(1); xs(0) = 1;
    arma::uvec ys(1); ys(0) = 2;

    mmf_obj.M(n,m);

    mmf_obj.M(n,m,nullptr,nullptr);
    mmf_obj.M(n(xs),m(ys),&xs,&ys);

    double n_x = n(1);
    double m_y = m(2);
    mmf_obj.M(n_x,m,&xs,nullptr);
    mmf_obj.M(n,m_y,nullptr,&ys);

    mmf_obj.dmu_x0(n,m);
    mmf_obj.dmu_0y(n,m);

    mmf_obj.dparams_M(n,m);
    arma::mat delta_params_M = alpha;
    mmf_obj.dparams_M(n,m,&delta_params_M);

    mmf_obj.Mx0(n);
    mmf_obj.M0y(m);

    mmf_obj.trans();

    //
    printf("\n*===================    End of mmfs::min Test    ===================*\n");
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

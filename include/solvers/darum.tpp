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
 * Deferred Acceptance (DA) algorithm for general ARUMs with NTU
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/04/2018
 */

// internal darum

template<typename Tg, typename Th, typename Tt, typename std::enable_if<!std::is_same<Tt,transfers::ntu>::value>::type*>
bool
darum_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, 
          const double err_tol, const uint_t max_iter)
{
    printf("darum only works for NTU transfers.\n");
    return false;
}

template<typename Tg, typename Th, typename Tt, typename std::enable_if<std::is_same<Tt,transfers::ntu>::value>::type*>
bool
darum_int(const dse<Tg,Th,Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, 
          const double err_tol, const uint_t max_iter)
{
    bool success = false;

    //

    const uint_t nbX = market.nbX;
    const uint_t nbY = market.nbY;

    arma::mat mu_NR = arma::max(market.n * arma::ones(1,nbY), arma::ones(nbX,1) * market.m.t());
    
    //

    uint_t iter = 0;
    double err = 2*err_tol;

    arma::mat U_P, U_D, mu_P, mu_D;

    while (err > err_tol && iter < max_iter)
    {
        iter++;

        //
        
        market.arums_G.Gbar(market.transfers_obj.alpha,mu_NR,market.n,U_P,mu_P);
        market.arums_H.Gbar(market.transfers_obj.gamma.t(),mu_P.t(),market.m,U_D,mu_D);

        arma::mat mu_diff = mu_P - mu_D.t();
        mu_NR -= mu_diff;

        //

        err = elem_max(arma::abs(mu_diff));
    }

    if (err <= err_tol && iter < max_iter) {
        success = true;
    }

    //
    // return equilibrium objects

    if (mu_out) {
        *mu_out = mu_D.t();
    }

    if (mu_x0_out) {
        *mu_x0_out = market.n - arma::sum(mu_D.t(),1);
    }
    if (mu_0y_out) {
        *mu_0y_out = market.m - arma::trans(arma::sum(mu_D.t(),0));
    }

    if (U_out) {
        *U_out = U_P;
    }
    if (V_out) {
        *V_out = U_D.t();
    }

    //

    return success;
}

// wrappers

template<typename Tg, typename Th, typename Tt>
bool
darum(const dse<Tg,Th,Tt>& market, arma::mat& mu_out)
{
    return darum_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
darum(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double err_tol_inp, const uint_t max_iter_inp)
{
    return darum_int(market,&mu_out,nullptr,nullptr,nullptr,nullptr,err_tol_inp,max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
darum(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    return darum_int(market,&mu_out,nullptr,nullptr,&U_out,&V_out);
}

template<typename Tg, typename Th, typename Tt>
bool
darum(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, 
      const double err_tol_inp, const uint_t max_iter_inp)
{
    return darum_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,err_tol_inp,max_iter_inp);
}

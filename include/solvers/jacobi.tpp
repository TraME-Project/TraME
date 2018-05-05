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
 * Jacobi solver
 *
 * Keith O'Hara
 * 08/25/2016
 *
 * This version:
 * 02/04/2018
 */

template<typename Tg, typename Th, typename Tt>
bool
jacobi_int(const dse<Tg,Th,Tt>& market, const arma::mat* w_low_inp, const arma::mat* w_up_inp, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, 
           const double err_tol, const uint_t max_iter)
{
    bool success = false;

    //

    if (market.need_norm) {
        printf("Jacobi does not yet allow for the case without unmatched agents.\n");
        return false;
    }

    //

    const uint_t nbX = market.nbX;
    const uint_t nbY = market.nbY;

    arma::mat mu_G, mu_H;

    //
    // w_up setup

    arma::mat w;

    if (!w_up_inp) {
        w = w_upper_bound(market);
    } else {
        w = *w_up_inp;

        arma::mat temp_UW = market.transfers_obj.UW(w);
        arma::mat temp_VW = market.transfers_obj.VW(w);

        market.arums_G.G(market.n,temp_UW,mu_G);
        market.arums_H.G(market.m,temp_VW.t(),mu_H);

        if (elem_min(mu_G - mu_H.t()) < 0) {
            printf("jacobi: w_up provided not an actual upper bound.\n");
            return false;
        }
    }

    //
    // w_low setup

    arma::mat w_low;
    if (!w_low_inp) {
        dse<Th,Tg,Tt> market_trans = market.trans();

        w_low = -arma::trans(w_upper_bound(market_trans));
    } else {
        w_low = *w_low_inp;

        arma::mat temp_UW = market.transfers_obj.UW(w_low);
        arma::mat temp_VW = market.transfers_obj.VW(w_low);

        market.arums_G.G(market.n,temp_UW,mu_G);
        market.arums_H.G(market.m,temp_VW.t(),mu_H);

        if (elem_max(mu_G - mu_H.t()) > 0) {
            printf("jacobi: w_low provided not an actual lower bound.\n");
            return false;
        }
    }

    //

    arma::mat U = market.transfers_obj.UW(w);
    arma::mat V = market.transfers_obj.VW(w);

    // market.arums_G.G(market.n,U,mu_G);
    // market.arums_H.G(market.m,V.t(),mu_H);

    // arma::mat Z = mu_G - mu_H.t();

    //

    uint_t iter = 0;
    double err = 2*err_tol;

    trame_jacobi_zeroin_data<Tg,Th,Tt> root_data;

    root_data.U = U;
    root_data.V = V;

    root_data.market_obj = market;

    const arma::mat norm_mat = market.n * market.m.t();

    while (err > err_tol && iter < max_iter) {
        iter++;

        for (uint_t x=0; x < nbX; x++) {
            root_data.x_ind = x;

            for (uint_t y=0; y < nbY; y++) {
                root_data.y_ind = y;

                w(x,y) = zeroin(w_low(x,y), w(x,y), jacobi_zeroin_fn<Tg,Th,Tt>, &root_data, nullptr, nullptr);

                U(x,y) = market.transfers_obj.UW(w(x,y),x,y);
                V(x,y) = market.transfers_obj.VW(w(x,y),x,y);

                root_data.U = U;
                root_data.V = V;
            }
        }

        //

        market.arums_G.G(market.n,U,mu_G);
        market.arums_H.G(market.m,V.t(),mu_H);
        
        err = elem_max( arma::abs(elem_div(mu_G - mu_H.t(),norm_mat)) );
    }

    if (err <= err_tol && iter < max_iter) {
        success = true;
    }

    //

    market.arums_G.G(market.n,U,mu_G);

    if (mu_out) {
        *mu_out = mu_G;
    }

    if (mu_x0_out) {
        *mu_x0_out = market.n - arma::sum(mu_G,1);
    }
    if (mu_0y_out) {
        *mu_0y_out = market.m - arma::trans(arma::sum(mu_G,0));
    }

    if (U_out) {
        *U_out = U;
    }
    if (V_out) {
        *V_out = V;
    }
    //
    return success;
}

//
// wrappers

template<typename Tg, typename Th, typename Tt>
bool
jacobi(const dse<Tg,Th,Tt>& market, arma::mat& mu_out)
{
    return jacobi_int(market,nullptr,nullptr,&mu_out,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
jacobi(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, const double err_tol_inp, const uint_t max_iter_inp)
{
    return jacobi_int(market,nullptr,nullptr,&mu_out,nullptr,nullptr,nullptr,nullptr,err_tol_inp,max_iter_inp);
}

template<typename Tg, typename Th, typename Tt>
bool
jacobi(const dse<Tg,Th,Tt>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    return jacobi_int(market,nullptr,nullptr,&mu_out,nullptr,nullptr,&U_out,&V_out);
}

template<typename Tg, typename Th, typename Tt>
bool
jacobi(const dse<Tg,Th,Tt>& market, const arma::mat& w_low_inp, const arma::mat& w_up_inp, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out,
       const double err_tol_inp, const uint_t max_iter_inp)
{
    return jacobi_int(market,&w_low_inp,&w_up_inp,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,err_tol_inp,max_iter_inp);
}

//
// root-finding function

template<typename Tg, typename Th, typename Tt>
double
jacobi_zeroin_fn(double z, void* opt_data)
{
    trame_jacobi_zeroin_data<Tg,Th,Tt> *d = reinterpret_cast<trame_jacobi_zeroin_data<Tg,Th,Tt>*>(opt_data);

    //

    arma::mat U = d->U;
    arma::mat V = d->V;
    
    U(d->x_ind,d->y_ind) = d->market_obj.transfers_obj.UW(z,d->x_ind,d->y_ind);
    V(d->x_ind,d->y_ind) = d->market_obj.transfers_obj.VW(z,d->x_ind,d->y_ind);

    arma::mat mu_G, mu_H;
    d->market_obj.arums_G.G(d->market_obj.n,U,mu_G);
    d->market_obj.arums_H.G(d->market_obj.m,V.t(),mu_H);

    // arma::mat Z = mu_G - mu_H.t();
    // return Z(d->x_ind,d->y_ind);
    
    return mu_G(d->x_ind,d->y_ind) - mu_H(d->y_ind,d->x_ind);
}

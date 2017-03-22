/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
  ##      Simon Weber
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
 * 02/15/2017
 */

template<typename Ta, typename Tm>
bool jacobi_int(const dse<Ta,Tm>& market, const arma::mat* w_low_inp, const arma::mat* w_up_inp, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, const double* tol_inp, const int* max_iter_inp)
{
    bool success = false;
    //
    if (market.need_norm) {
        printf("Jacobi does not yet allow for the case without unmatched agents.\n");
        return false;
    }
    //
    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;

    Ta arums_G = market.arums_G;
    Ta arums_H = market.arums_H;

    Tm trans_obj = market.trans_obj;

    double tol = (tol_inp) ? *tol_inp : 1E-04;
    int max_iter = (max_iter_inp) ? *max_iter_inp : 10000;

    arma::mat mu_G, mu_H;
    //
    // w_up setup
    arma::mat w;
    if (!w_up_inp) {
        w = w_upper_bound(market);
    } else {
        w = *w_up_inp;

        arma::mat temp_UW = trans_obj.UW(w);
        arma::mat temp_VW = trans_obj.VW(w);

        arums_G.G(n,temp_UW,mu_G);
        arums_H.G(m,temp_VW.t(),mu_H);

        arma::mat Z = mu_G - mu_H.t();

        if (elem_min(Z) < 0) {
            printf("jacobi: w_up provided not an actual upper bound.\n");
            return false;
        }
    }
    //
    // w_low setup
    arma::mat w_low;
    if (!w_low_inp) {
        dse<Ta,Tm> market_trans = market;
        market_trans.trans();

        w_low = -arma::trans(w_upper_bound(market_trans));
    } else {
        w_low = *w_low_inp;

        arma::mat temp_UW = trans_obj.UW(w_low);
        arma::mat temp_VW = trans_obj.VW(w_low);

        arums_G.G(n,temp_UW,mu_G);
        arums_H.G(m,temp_VW.t(),mu_H);

        arma::mat Z = mu_G - mu_H.t();

        if (elem_max(Z) > 0) {
            printf("jacobi: w_low provided not an actual upper bound.\n");
            return false;
        }
    }
    //
    arma::mat U = trans_obj.UW(w);
    arma::mat V = trans_obj.VW(w);

    arums_G.G(n,U,mu_G);
    arums_H.G(m,V.t(),mu_H);

    arma::mat Z = mu_G - mu_H.t();
    //
    int iter = 0;
    double err = 2*tol;

    int x = 0, y = 0;

    trame_jacobi_zeroin_data<Ta,Tm> root_data;

    root_data.x_ind = x;
    root_data.y_ind = y;

    root_data.n = n;
    root_data.m = m;

    root_data.U = U;
    root_data.V = V;

    root_data.arums_G = arums_G;
    root_data.arums_H = arums_H;
    root_data.trans_obj = trans_obj;

    arma::mat norm_mat = n * m.t();

    while (err > tol && iter < max_iter) {
        iter++;
        
        for (x=0; x < nbX; x++) {
            root_data.x_ind = x;

            for (y=0; y < nbY; y++) {
                root_data.y_ind = y;
        
                w(x,y) = zeroin(w_low(x,y), w(x,y), jacobi_zeroin_fn<Ta,Tm>, &root_data, NULL, NULL);
                U(x,y) = trans_obj.UW(w(x,y),x,y);
                V(x,y) = trans_obj.VW(w(x,y),x,y);

                root_data.U = U;
                root_data.V = V;
            }
        }
        //
        arums_G.G(n,U,mu_G);
        arums_H.G(m,V.t(),mu_H);

        Z = mu_G - mu_H.t();
        //
        err = elem_max(arma::abs(elem_div(Z,norm_mat)));
    }

    if (err <= tol && iter < max_iter) {
        success = true;
    }
    //
    arums_G.G(n,U,mu_G);

    if (mu_out) {
        *mu_out = mu_G;
    }

    if (mu_x0_out) {
        *mu_x0_out = n - arma::sum(mu_G,1);
    }
    if (mu_0y_out) {
        *mu_0y_out = m - arma::trans(arma::sum(mu_G,0));
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

// wrappers

template<typename Ta, typename Tm>
bool jacobi(const dse<Ta,Tm>& market, arma::mat& mu_out)
{
    bool res = jacobi_int(market,NULL,NULL,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL);
    
    return res;
}

template<typename Ta, typename Tm>
bool jacobi(const dse<Ta,Tm>& market, arma::mat& mu_out, const double& tol_inp)
{
    bool res = jacobi_int(market,NULL,NULL,&mu_out,NULL,NULL,NULL,NULL,&tol_inp,NULL);
    
    return res;
}

template<typename Ta, typename Tm>
bool jacobi(const dse<Ta,Tm>& market, arma::mat& mu_out, const int& max_iter_inp)
{
    bool res = jacobi_int(market,NULL,NULL,&mu_out,NULL,NULL,NULL,NULL,NULL,&max_iter_inp);
    
    return res;
}

template<typename Ta, typename Tm>
bool jacobi(const dse<Ta,Tm>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp)
{
    bool res = jacobi_int(market,NULL,NULL,&mu_out,NULL,NULL,NULL,NULL,&tol_inp,&max_iter_inp);
    
    return res;
}

template<typename Ta, typename Tm>
bool jacobi(const dse<Ta,Tm>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    bool res = jacobi_int(market,NULL,NULL,&mu_out,NULL,NULL,&U_out,&V_out,NULL,NULL);
    
    return res;
}

template<typename Ta, typename Tm>
bool jacobi(const dse<Ta,Tm>& market, const arma::mat& w_low_inp, const arma::mat& w_up_inp, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, const double* tol_inp, const int* max_iter_inp)
{
    bool res = jacobi_int(market,&w_low_inp,&w_up_inp,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,tol_inp,max_iter_inp);
    
    return res;
}

// root-finding function

template<typename Ta, typename Tm>
double jacobi_zeroin_fn(double z, void* opt_data)
{
    trame_jacobi_zeroin_data<Ta,Tm> *d = reinterpret_cast<trame_jacobi_zeroin_data<Ta,Tm>*>(opt_data);
    //
    double ret = 1.0;

    int x_ind = d->x_ind;
    int y_ind = d->y_ind;

    arma::mat U = d->U;
    arma::mat V = d->V;
    //
    U(x_ind,y_ind) = d->trans_obj.UW(z,x_ind,y_ind);
    V(x_ind,y_ind) = d->trans_obj.VW(z,x_ind,y_ind);
    
    arma::mat mu_G, mu_H;
    d->arums_G.G(d->n,U,mu_G);
    d->arums_H.G(d->m,V.t(),mu_H);

    arma::mat Z = mu_G - mu_H.t();
    ret = Z(x_ind,y_ind);
    //
    return ret;
}

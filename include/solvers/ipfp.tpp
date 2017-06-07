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
 * Iterated Proportional Fitting Procedure (IPFP) for MFE-type markets
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/15/2017
 */

// internal ipfp

template<typename Tt>
bool
ipfp_int(const mfe<Tt>& market, arma::mat* mu_out, arma::vec* mu_x0_out, arma::vec* mu_0y_out, arma::mat* U_out, arma::mat* V_out, arma::vec* u_out, arma::vec* v_out, const double* tol_inp, const int* max_iter_inp, const arma::vec* by_start)
{
    bool success = false;
    //
    //bool noSingles = market.need_norm;

    int nbX = market.nbX;
    int nbY = market.nbY;

    double tol = (tol_inp) ? *tol_inp : 1E-12;
    int max_iter = (max_iter_inp) ? *max_iter_inp : 10000;
    //
    // begin loop
    arma::vec ax(nbX);
    arma::vec by = (by_start) ? *by_start : market.m;

    arma::vec val_old(nbX+nbY);
    arma::vec val_new(nbX+nbY);
    arma::vec val_err(nbX+nbY);

    int iter = 0;
    double err = 2*tol;
    
    while (err > tol && iter < max_iter) {
        iter++;

        val_old = arma::join_cols(ax,by);

        // Solve for 'ax' and then 'by'
        ax = market.marg_x_inv(by);
        by = market.marg_y_inv(ax);

        /* Keith: need to add this later
        if (noSingles) {
            
        }
        */

        val_new = arma::join_cols(ax,by);
        val_err = arma::abs(val_new - val_old);

        err = arma::as_scalar(arma::max(val_err));
    }

    if (err <= tol && iter < max_iter) {
        success = true;
    }
    //
    // Construct the equilibrium outcome based on 'ax' and 'by' obtained above
    if (mu_out || mu_x0_out || mu_0y_out || U_out || V_out || u_out || v_out) {

        arma::mat mu = market.mmfs_obj.M(ax,by);

        if (mu_out) {
            *mu_out = mu;
        }
        //
        if (mu_x0_out || mu_0y_out || U_out || V_out || u_out || v_out) {
            arma::vec mu_x0 = market.mmfs_obj.Mx0(ax);
            arma::vec mu_0y = market.mmfs_obj.M0y(by);

            if (mu_x0_out) {
                *mu_x0_out = mu_x0;
            }
            if (mu_0y_out) {
                *mu_0y_out = mu_0y;
            }
            
            if (U_out) {
                //*U_out = arma::log(mu / arma::repmat(mu_x0,1,mu.n_cols));
                *U_out = arma::log(elem_div(mu,mu_x0));
            }
            if (V_out) {
                //*V_out = arma::trans(arma::log(mu.t() / arma::repmat(mu_0y,1,mu.n_rows)));
                *V_out = arma::trans(arma::log(elem_div(mu.t(),mu_0y)));
            }

            if (u_out) { 
                *u_out = - arma::log(mu_x0);
            }
            if (u_out) {
                *v_out = - arma::log(mu_0y);
            }
        }
    }
    //
    return success;
}

//
// wrappers

template<typename Tt>
bool
ipfp(const mfe<Tt>& market, arma::mat& mu_out)
{
    bool res = ipfp_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
    
    return res;
}

template<typename Tt>
bool
ipfp(const mfe<Tt>& market, arma::mat& mu_out, const double& tol_inp)
{
    bool res = ipfp_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&tol_inp,NULL,NULL);
    
    return res;
}

template<typename Tt>
bool
ipfp(const mfe<Tt>& market, arma::mat& mu_out, const int& max_iter_inp)
{
    bool res = ipfp_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL,&max_iter_inp,NULL);
    
    return res;
}

template<typename Tt>
bool
ipfp(const mfe<Tt>& market, arma::mat& mu_out, const arma::vec& by_start)
{
    bool res = ipfp_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,&by_start);
    
    return res;
}

template<typename Tt>
bool
ipfp(const mfe<Tt>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp)
{
    bool res = ipfp_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&tol_inp,&max_iter_inp,NULL);
    
    return res;
}

template<typename Tt>
bool
ipfp(const mfe<Tt>& market, arma::mat& mu_out, const double& tol_inp, const arma::vec& by_start)
{
    bool res = ipfp_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&tol_inp,NULL,&by_start);
    
    return res;
}

template<typename Tt>
bool
ipfp(const mfe<Tt>& market, arma::mat& mu_out, const int& max_iter_inp, const arma::vec& by_start)
{
    bool res = ipfp_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,NULL,&max_iter_inp,&by_start);
    
    return res;
}

template<typename Tt>
bool
ipfp(const mfe<Tt>& market, arma::mat& mu_out, const double& tol_inp, const int& max_iter_inp, const arma::vec& by_start)
{
    bool res = ipfp_int(market,&mu_out,NULL,NULL,NULL,NULL,NULL,NULL,&tol_inp,&max_iter_inp,&by_start);
    
    return res;
}

template<typename Tt>
bool
ipfp(const mfe<Tt>& market, arma::mat& mu_out, arma::mat& U_out, arma::mat& V_out)
{
    bool res = ipfp_int(market,&mu_out,NULL,NULL,&U_out,&V_out,NULL,NULL,NULL,NULL,NULL);
    
    return res;
}

template<typename Tt>
bool
ipfp(const mfe<Tt>& market, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& U_out, arma::mat& V_out, arma::vec& u_out, arma::vec& v_out, const double* tol_inp, const int* max_iter_inp, const arma::vec* by_start)
{
    bool res = ipfp_int(market,&mu_out,&mu_x0_out,&mu_0y_out,&U_out,&V_out,&u_out,&v_out,tol_inp,max_iter_inp,by_start);
    
    return res;
}

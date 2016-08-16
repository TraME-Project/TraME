/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 Alfred Galichon and the TraME Team
  ##
  ##   This file is part of the R package TraME.
  ##
  ##   The R package TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

// ipfp for logit

bool ipfp (MFE market, bool xFirst, double* tol_inp, arma::vec* by_start, arma::mat& mu, arma::vec& mux0, arma::vec& mu0y, arma::mat& U, arma::mat& V, arma::vec& u, arma::vec& v)
{
    MMF mmf_obj = marker.mmf_obj;

    bool noSingles = need_norm;

    arma::vec n = mmf_obj.n;
    arma::vec m = mmf_obj.m;

    int nbX = mmf_obj.nbX;
    int nbY = mmf_obj.nbY;
    //
    // begin loop
    arma::vec by;
    if (by_start) {
        by = *by_start; 
    } else {
        by = m;
    }

    arma::vec ax(nbX);
    arma::vec val(nbX+nbY);

    double tol;
    if (tol_inp) {
        tol = *tol_inp;
    } else {
        tol = 1E-12;
    }

    double err = 2*tol;
    int iter = 0;
    //
    while (err > tol) {
        iter++;
        val = arma::join_cols(ax,bx);

        // Solve for ax and then by
        ax = mmf_obj.marg_x_inv(NULL,by);
        by = mmf_obj.marg_y_inv(NULL,ax);

        /* Keith: need to add this later
        if (noSingles) {
            
        }
        */

        err = arma::as_scalar(arma::max(arma::abs(arma::join_cols(ax,bs) - val)));
    }
    //
    // Construct the equilibrium outcome based on ax and by obtained from above
    mu = mmf_obj.M(ax,by,NULL,NULL);
    mux0 = mmf_obj.Mx0(ax);
    mu0y = mmf_obj.M0y(by);

    U = arma::log(mu / mux0);
    V = arma::trans(arma::log(mu.t() / mu0y));

    u = - arma::log(mux0);
    v = - arma::log(mu0y);
    //
}

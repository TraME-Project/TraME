/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 the TraME Team:
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
 * general model class
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * This version:
 * 11/26/2016
 */

template<typename Ta>
void model<Ta>::build(const arma::mat& X_inp, const arma::mat& Y_inp)
{
    this->build(X_inp,Y_inp,NULL,NULL,NULL);
}

template<typename Ta>
void model<Ta>::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    need_norm = false;

    nbX = X_inp.n_rows;
    nbY = Y_inp.n_rows;

    dX = X_inp.n_cols;
    dY = Y_inp.n_cols;

    nbParams = dX*dY;
    //
    if (n_inp) {
        n = *n_inp;
    } else {
        n.ones(nbX,1);
    }
    if (m_inp) {
        m = *m_inp;
    } else {
        m.ones(nbY,1);
    }

    phi_xyk = arma::kron(Y_inp,X_inp);
    //
}

template<typename Ta>
void dse<Ta>::build_market_TU(const arma::mat& theta)
{
    market_obj.build_TU(n,m,phi_xyk,need_norm);
}

// general arums input
template<typename Ta>
void model<Ta>::build_market_TU(const arma::mat& theta, const Ta& arums_G_inp, const Ta& arums_H_inp)
{
    market_obj.build_TU(n,m,phi_xyk,arums_G_inp,arums_H_inp,need_norm);
}

// empirical version
template<typename Ta> template<typename T>
void model<Ta>::build_TU(const arma::mat& theta, T arums_G_inp, T arums_H_inp, int nbDraws, int seed)
{
    market_obj.build_TU(n,m,phi_xyk,arums_G_inp,arums_H_inp,nbDraws,seed,need_norm);
}

template<typename Ta>
bool model<Ta>::mme(arma::mat& mu_hat)
{
    printf("general mme not defined\n");
    return false;
}

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
 * general model class
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * This version:
 * 03/22/2017
 */

// first method to build
inline 
void model_base::build(const arma::cube& phi_xyk_inp)
{
    this->build_int(phi_xyk_inp,NULL,NULL);
}

inline 
void model_base::build(const arma::cube& phi_xyk_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    this->build_int(phi_xyk_inp,&n_inp,&m_inp);
}

inline 
void model_base::build_int(const arma::cube& phi_xyk_inp, const arma::vec* n_inp, const arma::vec* m_inp)
{
    need_norm = false;

    nbX = phi_xyk_inp.n_rows;
    nbY = phi_xyk_inp.n_cols;
    nbParams = phi_xyk_inp.n_slices;
    //
    n = (n_inp) ? *n_inp : arma::ones(nbX,1);
    m = (m_inp) ? *m_inp : arma::ones(nbY,1);

    phi_xyk = phi_xyk_inp;
    //
}

// second method to build
inline 
void model_base::build(const arma::mat& X_inp, const arma::mat& Y_inp)
{
    this->build_int(X_inp,Y_inp,NULL,NULL);
}

inline 
void model_base::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    this->build_int(X_inp,Y_inp,&n_inp,&m_inp);
}

inline 
void model_base::build_int(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec* n_inp, const arma::vec* m_inp)
{
    need_norm = false;

    nbX = X_inp.n_rows;
    nbY = Y_inp.n_rows;

    dX = X_inp.n_cols;
    dY = Y_inp.n_cols;

    nbParams = dX*dY;
    //
    n = (n_inp) ? *n_inp : arma::ones(nbX,1);
    m = (m_inp) ? *m_inp : arma::ones(nbY,1);
    //
    arma::mat phi_xy_temp = arma::kron(Y_inp,X_inp);
    arma::cube phi_xyk_temp(phi_xy_temp.memptr(),nbX,nbY,nbParams,false); // share memory

    phi_xyk = phi_xyk_temp;
}

// gradients

inline 
void model_base::dparam(const arma::mat* dparams_inp, arma::mat& dparamsPsi_out)
{
    this->dparam(dparams_inp,dparamsPsi_out,NULL,NULL);
}

inline 
void model_base::dparam(const arma::mat* dparams_inp, arma::mat& dparamsPsi_out, arma::mat* dparamsG_out, arma::mat* dparamsH_out)
{
    arma::mat dparams_mat = (dparams_inp) ? *dparams_inp : arma::eye(nbParams,nbParams);

    dparamsPsi_out = Phi_xy() * dparams_mat;
    //
    if (dparamsG_out) {
        *dparamsG_out = arma::zeros(0,dparams_mat.n_cols);
    }
    if (dparamsH_out) {
        *dparamsH_out = arma::zeros(0,dparams_mat.n_cols);
    }
}

// Keith: should probably switch this to be a member variable
inline
arma::mat model_base::Phi_xy()
{
    // mirror R's approach to creating a matrix from an array; take each slice and vectorise that matrix
    arma::mat phi_xy_mat(nbX*nbY,nbParams);
    for(int k = 0; k < nbParams; k++) {
        phi_xy_mat.col(k) = arma::vectorise(phi_xyk.slice(k));
    }

    return phi_xy_mat;
}

inline
arma::mat model_base::Phi_xy_theta(const arma::mat& theta)
{
    arma::mat ret = arma::reshape(Phi_xy() * theta,nbX,nbY);
    return ret;
}

inline
void model_base::init_param(arma::mat& params)
{
    params.zeros(nbParams,1);
}

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
 * Derived classes to provide wrappers to the TraME library
 *
 * Keith O'Hara
 * 11/10/2016
 */

template<typename Ta>
void dse_empirical_R::build_ETU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, Ta arums_G_inp, Ta arums_H_inp, bool need_norm_inp)
{
    try {
        trame::empirical empirical_1 = arums_G_inp.simul();
        trame::empirical empirical_2 = arums_H_inp.simul();

        this->build_LTU(n_inp,m_inp,lambda_inp,phi_inp,empirical_1,empirical_2,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

template<typename Ta>
void dse_empirical_R::build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, Ta arums_G_inp, Ta arums_H_inp, bool need_norm_inp)
{
    try {
        trame::empirical empirical_1 = arums_G_inp.simul();
        trame::empirical empirical_2 = arums_H_inp.simul();

        this->build_NTU(n_inp,m_inp,alpha_inp,gamma_inp,empirical_1,empirical_2,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

template<typename Ta>
void dse_empirical_R::build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, Ta arums_G_inp, Ta arums_H_inp, bool need_norm_inp)
{
    try {
        trame::empirical empirical_1 = static_cast<trame::empirical&>(arums_G_inp);
        trame::empirical empirical_2 = static_cast<trame::empirical&>(arums_H_inp);

        this->build_TU(n_inp,m_inp,phi_inp,empirical_1,empirical_2,need_norm_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "trame: C++ exception (unknown reason)" );
    }
}

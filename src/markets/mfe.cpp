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
 * Matching Function Equilibrium (MFE) market
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 06/06/2017
 */

#include "trame.hpp"

namespace trame
{

//
// builds

// template<>
// void
// mfe<mmfs::ces>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp)
// {
//     nbX = n_inp.n_elem;
//     nbY = m_inp.n_elem;
    
//     n = n_inp;
//     m = m_inp;

//     mmfs_obj.build(alpha_inp,gamma_inp,tau_inp,need_norm);
// }

// template<>
// void
// mfe<mmfs::cd>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp)
// {
//     nbX = n_inp.n_elem;
//     nbY = m_inp.n_elem;
    
//     n = n_inp;
//     m = m_inp;
    
//     mmfs_obj.build(lambda_inp,phi_inp/sigma,need_norm);
// }

// template<>
// void
// mfe<mmfs::min>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp)
// {
//     nbX = n_inp.n_elem;
//     nbY = m_inp.n_elem;
    
//     n = n_inp;
//     m = m_inp;

//     mmfs_obj.build(alpha_inp/sigma,gamma_inp/sigma,need_norm);
// }

// template<>
// void
// mfe<mmfs::geo>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp)
// {
//     nbX = n_inp.n_elem;
//     nbY = m_inp.n_elem;
    
//     n = n_inp;
//     m = m_inp;

//     mmfs_obj.build(phi_inp/sigma,need_norm);
// }

}

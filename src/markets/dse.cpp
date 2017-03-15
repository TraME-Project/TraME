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
 * Demand-Supply Equilibrium (DSE) market
 * template specialization
 *
 * Keith O'Hara
 * 08/17/2016
 *
 * This version:
 * 03/14/2017
 */

#include "trame.hpp"

namespace trame
{

// we specialize because cupids_lp is only define for empirical classes
template<>
bool dse<empirical,tu>::solve(arma::mat& mu_sol, const char* solver)
{
    bool res = false;
    const char sig = (solver != NULL) ? solver[0] : char(0);
    
    if (solver) { // not NULL
        if (sig=='c') {
            res = cupids_lp(*this,mu_sol);
        }
        if (sig=='e') {
            res = eap_nash(*this,mu_sol);
        }
        if (sig=='j') {
            res = jacobi(*this,mu_sol);
        }
        if (sig=='m') {
            res = max_welfare(*this,mu_sol);
        }
        if (sig=='o') {
            res = oap_lp(*this,mu_sol);
        }
    } /*else {
        if (NTU) {
            res = darum(*this,mu_sol);
        } else if (TU) {
            res = max_welfare(*this,mu_sol);
        } else {
            res = jacobi(*this,mu_sol);
        }
    }*/
    //
    return res;
}

}

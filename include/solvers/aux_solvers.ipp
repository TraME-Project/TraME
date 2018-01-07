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
 * auxiliary functions for equilibrium solvers
 *
 * Keith O'Hara
 * 08/23/2016
 *
 * This version:
 * 07/26/2016
 */

/*
 * build disaggregate epsilon function; used in cupids_lp
 * takes U_xy and arums as input; returns U_iy as output.
 */

inline
int 
build_disaggregate_epsilon(const arma::vec& n, const trame::arums::empirical& arums_emp_inp, arma::mat& epsilon_iy, arma::mat& epsilon0_i, arma::mat& I_ix)
{
    const int nbX = arums_emp_inp.nbX;
    const int nbY = arums_emp_inp.nbY;

    const int n_draws = arums_emp_inp.aux_n_draws;
    const int nbI = nbX * n_draws;

    arma::vec I_01(nbX);
    I_ix.zeros(nbI,nbX);

    arma::mat epsilons = arma::zeros(nbI,nbY+1);

    //

    for (int x=0; x < nbX; x++) {
        arma::mat epsilon;

        if (arums_emp_inp.x_homogeneous) {
            epsilon = arums_emp_inp.atoms.slice(0);
        } else {
            epsilon = arums_emp_inp.atoms.slice(x);
        }

        epsilons.rows(x*n_draws,(x+1)*n_draws-1) = std::move(epsilon);
        
        //

        I_01.zeros();
        I_01(x) = 1;
        
        I_ix.rows(x*n_draws,(x+1)*n_draws-1) = arma::repmat(I_01.t(),n_draws,1); // Keith: check use of byrow here
    }

    //

    epsilon_iy = std::move(epsilons.cols(0,nbY-1));
    epsilon0_i = std::move(epsilons.col(nbY));
    
    //

    return n_draws;
}

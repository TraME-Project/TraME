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
 * probit additive random utility model (ARUM) class
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 06/10/2017
 */

#include "trame.hpp"

//
// build functions

trame::arums::probit::probit(int nbX_inp, int nbY_inp)
{   
    this->build(nbX_inp, nbY_inp);
}

trame::arums::probit::probit(int nbX_inp, int nbY_inp, bool outside_option_inp)
{   
    this->build_int(nbX_inp, nbY_inp, NULL, outside_option_inp);
}

trame::arums::probit::probit(int nbX_inp, int nbY_inp, double rho_inp, bool outside_option_inp)
{   
    this->build_int(nbX_inp, nbY_inp, &rho_inp, outside_option_inp);
}

void
trame::arums::probit::build(int nbX_inp, int nbY_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
}

void
trame::arums::probit::build(int nbX_inp, int nbY_inp, bool outside_option_inp)
{   
    this->build_int(nbX_inp, nbY_inp, NULL, outside_option_inp);
}

void
trame::arums::probit::build(int nbX_inp, int nbY_inp, double rho_inp, bool outside_option_inp)
{   
    this->build_int(nbX_inp, nbY_inp, &rho_inp, outside_option_inp);
}

void
trame::arums::probit::build_int(int nbX_inp, int nbY_inp, double* rho_inp, bool outside_option_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
    outside_option = outside_option_inp;

    if (rho_inp) {
        rho = *rho_inp;
    }
    //
    if (outside_option_inp) {
        aux_nb_options = nbY + 1;
    } else {
        aux_nb_options = nbY;
    }
    dim_params = (nbX_inp * aux_nb_options * (aux_nb_options-1))/2;
}

void
trame::arums::probit::unifCorrelCovMatrices()
{
    this->unifCorrelCovMatrices(rho);
}

void
trame::arums::probit::unifCorrelCovMatrices(double rho_inp)
{
    arma::mat Sig = rho_inp * arma::ones(aux_nb_options,aux_nb_options) + (1-rho_inp) * arma::eye(aux_nb_options,aux_nb_options);
    
    if (outside_option) {
        Sig.col(aux_nb_options-1).fill(0);
        Sig.row(aux_nb_options-1).fill(0);
        Sig(aux_nb_options-1,aux_nb_options-1) = 1;
    }
    //
    Covar.set_size(aux_nb_options,aux_nb_options,nbX); // note: this is different to the R code

    for (int i=0; i<nbX; i++) {
        Covar.slice(i) = Sig;
    }
    //
}

//
// simulation

trame::arums::empirical
trame::arums::probit::simul()
const
{
    empirical emp_obj;
    
    this->simul_int(emp_obj,NULL,NULL);
    //
    return emp_obj;
}

trame::arums::empirical
trame::arums::probit::simul(int nbDraws, int seed)
const
{
    empirical emp_obj;
    
    this->simul_int(emp_obj,&nbDraws,&seed);
    //
    return emp_obj;
}

void
trame::arums::probit::simul(empirical& obj_out)
const
{
    this->simul_int(obj_out,NULL,NULL);
}

void
trame::arums::probit::simul(empirical& obj_out, int nbDraws, int seed)
const
{
    this->simul_int(obj_out,&nbDraws,&seed);
}

void
trame::arums::probit::simul_int(empirical& obj_out, int* nbDraws, int* seed_val)
const
{
    int n_draws = 0;
    
    if (nbDraws) {
        n_draws = *nbDraws;
    } else {
#ifdef TRAME_DEFAULT_SIM_DRAWS
        n_draws = TRAME_DEFAULT_SIM_DRAWS;
#else
        n_draws = 1000;
#endif
    }
    //
    if (seed_val) {
        arma::arma_rng::set_seed(*seed_val);
    }
    //
    arma::vec V;
    arma::mat Q, SqrtCovar;
    arma::cube atoms(n_draws,aux_nb_options,nbX);
    
    for (int j=0; j<nbX; j++) {
        eig_sym(V, Q, Covar.slice(j));
        SqrtCovar = Q * arma::diagmat(1.0/arma::sqrt(V)) * Q.t();
        //
        atoms.slice(j) = arma::randn(n_draws,aux_nb_options) * SqrtCovar;
    }
    //
    obj_out.nbX = nbX;
    obj_out.nbY = nbY;
    obj_out.dim_params = atoms.n_elem;
    obj_out.atoms = atoms;
    obj_out.aux_nbDraws = n_draws;
    obj_out.x_homogeneous = false;
    obj_out.outside_option = outside_option;

    if (outside_option) {
        obj_out.nb_options = nbY + 1;
    } else {
        obj_out.nb_options = nbY;
    }
    //
    if (seed_val) {
        arma::arma_rng::set_seed_random(); // need to reset the seed
    }
}

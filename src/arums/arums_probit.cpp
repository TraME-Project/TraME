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
 * 07/25/2017
 */

#include "ancillary/ancillary.hpp"
#include "arums/arums.hpp"

//
// build functions

trame::arums::probit::probit(const int nbX_inp, const int nbY_inp)
{   
    this->build(nbX_inp, nbY_inp);
}

trame::arums::probit::probit(const int nbX_inp, const int nbY_inp, const bool outside_option_inp)
{   
    this->build_int(nbX_inp, nbY_inp, nullptr, outside_option_inp);
}

trame::arums::probit::probit(const int nbX_inp, const int nbY_inp, const double rho_inp, const bool outside_option_inp)
{   
    this->build_int(nbX_inp, nbY_inp, &rho_inp, outside_option_inp);
}

void
trame::arums::probit::build(const int nbX_inp, const int nbY_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
}

void
trame::arums::probit::build(const int nbX_inp, const int nbY_inp, const bool outside_option_inp)
{   
    this->build_int(nbX_inp, nbY_inp, nullptr, outside_option_inp);
}

void
trame::arums::probit::build(const int nbX_inp, const int nbY_inp, const double rho_inp, const bool outside_option_inp)
{   
    this->build_int(nbX_inp, nbY_inp, &rho_inp, outside_option_inp);
}

void
trame::arums::probit::build_int(const int nbX_inp, const int nbY_inp, const double* rho_inp, const bool outside_option_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
    outside_option = outside_option_inp;

    rho = (rho_inp) ? *rho_inp : 0.0;
    //
    aux_nb_options = (outside_option_inp) ? nbY + 1 : nbY;

    dim_params = (nbX_inp * aux_nb_options * (aux_nb_options-1))/2; // Keith: check this
}

void
trame::arums::probit::unifCorrelCovMatrices()
{
    this->unifCorrelCovMatrices(rho);
}

void
trame::arums::probit::unifCorrelCovMatrices(const double rho_inp)
{
    arma::mat Sig = rho_inp * arma::ones(aux_nb_options,aux_nb_options) + (1.0 - rho_inp) * arma::eye(aux_nb_options,aux_nb_options);
    
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
    
    this->simul_int(emp_obj,nullptr,nullptr);
    //
    return emp_obj;
}

trame::arums::empirical
trame::arums::probit::simul(const int n_draws, const int seed)
const
{
    empirical emp_obj;
    
    this->simul_int(emp_obj,&n_draws,&seed);
    //
    return emp_obj;
}

void
trame::arums::probit::simul(empirical& obj_out)
const
{
    this->simul_int(obj_out,nullptr,nullptr);
}

void
trame::arums::probit::simul(empirical& obj_out, const int n_draws, const int seed)
const
{
    this->simul_int(obj_out,&n_draws,&seed);
}

void
trame::arums::probit::simul_int(empirical& obj_out, const int* n_draws_inp, const int* seed_val)
const
{
    int n_draws = 0;
    
    if (n_draws_inp) {
        n_draws = *n_draws_inp;
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
    obj_out.build(nbX,nbY,atoms,false,outside_option);
    //
    if (seed_val) {
        arma::arma_rng::set_seed_random(); // need to reset the seed
    }
}

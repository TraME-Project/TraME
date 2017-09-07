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
 * No hetergeneity (none) additive random utility model (ARUM) class
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

trame::arums::none::none(const int nbX_inp, const int nbY_inp)
{   
    this->build(nbX_inp,nbY_inp);
}

void
trame::arums::none::build(const int nbX_inp, const int nbY_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
    dim_params = 0;
}

//
// indirect utility

double
trame::arums::none::G(const arma::vec& n)
{   
    return this->G(n,U,mu_sol);
}

double
trame::arums::none::G(const arma::vec& n, const arma::mat& U_inp, arma::mat& mu_out)
const
{   
    double val = 0.0;
    
    mu_out.set_size(nbX,nbY);
    arma::mat mu_x_temp;
    //
    for (int i=0; i<nbX; i++) {
        double val_x = Gx(U_inp.row(i).t(),mu_x_temp);
        //
        val += n(i)*val_x;
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double
trame::arums::none::Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out)
const
{
    const arma::uvec temp_vec = which_max(U_x_inp, 0);
    const int y = temp_vec(0);
    //
    mu_x_out.zeros(nbY,1);
    
    if (y < nbY) {
        mu_x_out(y) = 1;
    }
    //
    return std::max(elem_max(U_x_inp), 0.0);
}

// just to conform with other arums classes
double 
trame::arums::none::Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, const int x)
const
{
    double val_x = 0.0;

    if (U_x_inp.n_rows > 1 && U_x_inp.n_cols > 1) {
        val_x = this->Gx(U_x_inp.row(x).t(),mu_x_out);
    } else if (U_x_inp.n_rows == 1 && U_x_inp.n_cols > 1) { 
        val_x = this->Gx(U_x_inp.t(),mu_x_out);
    } else {
        val_x = this->Gx(U_x_inp,mu_x_out);
    }
    //
    return val_x;
}

//
// Fenchel transform of G

double
trame::arums::none::Gstar(const arma::vec& n)
{   
    printf("Gstar not yet defined for the no arums case.\n");

    return 0.0;
}

double
trame::arums::none::Gstar(const arma::vec& n, const arma::mat& mu_inp, arma::mat& U_out)
const
{   
    printf("Gstar not yet defined for the no arums case.\n");

    return 0.0;
}

double
trame::arums::none::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, const int x)
const
{   
    printf("Gstarx not yet defined for the no arums case.\n");

    return 0.0;
}

//
// Gbar is used by DARUM

double
trame::arums::none::Gbar(const arma::mat& Ubar, const arma::mat& mubar, const arma::vec& n, arma::mat& U_out, arma::mat& mu_out)
const
{   
    double val = 0.0;
    
    U_out.set_size(nbX,nbY);
    mu_out.set_size(nbX,nbY);
    //
    arma::mat U_x_temp, mu_x_temp;
    
    for (int i=0; i<nbX; i++) {
        double val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),U_x_temp,mu_x_temp);
        //
        val += n(i)*val_temp;
        U_out.row(i) = arma::trans(U_x_temp);
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double
trame::arums::none::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out)
const
{
    int count_int = 0;
    const int nb_y0 = Ubar_x.n_elem;

    mu_x_out.set_size(nb_y0,1);
    //
    arma::uvec srt_ind = arma::sort_index(Ubar_x,"descend");
    //
    mu_x_out.set_size(nb_y0,1);
    double cumul = arma::as_scalar(mubar_x(srt_ind(count_int)));
    
    while ((count_int < nb_y0-1) && (cumul < 1.0) && (Ubar_x(srt_ind(count_int)) > 0)) {
        mu_x_out(srt_ind(count_int)) = mubar_x(srt_ind(count_int));

        count_int++;

        cumul += mubar_x(srt_ind(count_int)); // Keith: is this in the correct place?
    }
    //
    if (Ubar_x(srt_ind(count_int)) > 0) {
        mu_x_out(srt_ind(count_int)) = mubar_x(srt_ind(count_int)) + 1 - cumul;
    }

    U_x_out = arma::zeros(nb_y0,1);
    //
    return arma::accu(mu_x_out % Ubar_x);
}

// just to conform with other arums classes
double
trame::arums::none::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, const int x)
const
{
    return this->Gbarx(Ubar_x, mubar_x, U_x_out, mu_x_out);
}

//
// D2G

arma::mat
trame::arums::none::D2G(const arma::vec& n, const bool x_first)
const
{   
    arma::mat ret;

    printf("D2G not yet defined for the no arums case.\n");

    return ret;
}

void
trame::arums::none::D2G(arma::mat &H, const arma::vec& n, const bool x_first)
const
{   
    printf("D2G not yet defined for the no arums case.\n");
}

arma::mat
trame::arums::none::D2G(const arma::vec& n, const arma::mat& U_inp, const bool x_first)
const
{   
    arma::mat ret;

    printf("D2G not yet defined for the no arums case.\n");

    return ret;
}

void
trame::arums::none::D2G(arma::mat &H, const arma::vec& n, const arma::mat& U_inp, const bool x_first)
const
{   
    printf("D2G not yet defined for the no arums case.\n");
}

//
// D2Gstar

arma::mat
trame::arums::none::D2Gstar(const arma::vec& n, const bool x_first)
const
{   
    arma::mat ret;

    printf("D2Gstar not yet defined for the no arums case.\n");

    return ret;
}

void
trame::arums::none::D2Gstar(arma::mat &H, const arma::vec& n, const bool x_first)
const
{   
    printf("D2Gstar not yet defined for the no arums case.\n");
}

arma::mat
trame::arums::none::D2Gstar(const arma::vec& n, const arma::mat& mu_inp, const bool x_first)
const
{   
    arma::mat ret;

    printf("D2Gstar not yet defined for the no arums case.\n");

    return ret;
}

void
trame::arums::none::D2Gstar(arma::mat &H, const arma::vec& n, const arma::mat& mu_inp, const bool x_first)
const
{   
    printf("D2Gstar not yet defined for the no arums case.\n");
}

//
// dparams gradient

arma::mat
trame::arums::none::dparams_NablaGstar(const arma::vec& n, const arma::mat* dparams_inp, const bool x_first)
const
{
    arma::mat ret;
    this->dparams_NablaGstar(ret,n,mu_sol,dparams_inp,x_first);
    //
    return ret;
}

void
trame::arums::none::dparams_NablaGstar(arma::mat &ret, const arma::vec& n, const arma::mat* dparams_inp, const bool x_first)
const
{
    this->dparams_NablaGstar(ret,n,mu_sol,dparams_inp,x_first);
}

arma::mat
trame::arums::none::dparams_NablaGstar(const arma::vec& n, const arma::mat& mu_inp, const arma::mat* dparams_inp, const bool x_first)
const
{
    arma::mat ret;
    this->dparams_NablaGstar(ret,n,mu_inp,dparams_inp,x_first);
    //
    return ret;
}

void
trame::arums::none::dparams_NablaGstar(arma::mat &ret, const arma::vec& n, const arma::mat& mu_inp, const arma::mat* dparams_inp, const bool x_first)
const
{
    ret = arma::zeros(nbX*nbY,1);
}

//
// simulation

trame::arums::empirical
trame::arums::none::simul()
const
{
    empirical emp_obj;
    
    this->simul_int(emp_obj,nullptr,nullptr);
    //
    return emp_obj;
}

trame::arums::empirical
trame::arums::none::simul(const int n_draws, const int seed)
const
{
    empirical emp_obj;
    
    this->simul_int(emp_obj,&n_draws,&seed);
    //
    return emp_obj;
}

void
trame::arums::none::simul(empirical& obj_out)
const
{
    this->simul_int(obj_out,nullptr,nullptr);
}

void
trame::arums::none::simul(empirical& obj_out, const int n_draws, const int seed)
const
{
    this->simul_int(obj_out,&n_draws,&seed);
}

void
trame::arums::none::simul_int(empirical& obj_out, const int* n_draws_inp, const int* seed_val)
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
    obj_out.build(nbX,nbY,arma::zeros(n_draws,nbY+1,nbX),false,false); // Keith: is the 'outside_option = false' ?
    //
    if (seed_val) {
        arma::arma_rng::set_seed_random(); // need to reset the seed
    }
}

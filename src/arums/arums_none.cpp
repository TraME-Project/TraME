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
 * none additive random utility model (ARUM) class
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 02/11/2017
 */

#include "trame.hpp"

trame::arums::none::none(int nbX_inp, int nbY_inp)
{   
    this->build(nbX_inp,nbY_inp);
}

void trame::arums::none::build(int nbX_inp, int nbY_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
    nbParams = 0;
}

double trame::arums::none::G(const arma::vec& n)
{   
    double val = this->G(n,U,mu_sol);
    //
    return val;
}

double trame::arums::none::G(const arma::vec& n, const arma::mat& U_inp, arma::mat& mu_out)
{   
    double val=0.0, val_x;
    
    mu_out.set_size(nbX,nbY);
    arma::mat mu_x_temp;
    //
    for (int i=0; i<nbX; i++) {
        val_x = Gx(U_inp.row(i).t(),mu_x_temp);
        //
        val += n(i)*val_x;
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double trame::arums::none::Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out)
{
    arma::uvec temp_vec = which_max(U_x_inp, 0);
    int y = temp_vec(0);
    //
    mu_x_out.zeros(nbY,1);
    
    if (y < nbY) {
        mu_x_out(y) = 1;
    }
    //
    double val_x = std::max(elem_max(U_x_inp), (double) 0.0);
    //
    return val_x;
}

// just to conform with other arums classes
double trame::arums::none::Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, int x)
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

double trame::arums::none::Gstar(const arma::vec& n)
{   
    printf("Gstar not yet defined for no arums case.\n");

    return 0.0;
}

double trame::arums::none::Gstar(const arma::vec& n, const arma::mat& mu_inp, arma::mat& U_out)
{   
    printf("Gstar not yet defined for no arums case.\n");

    return 0.0;
}

double trame::arums::none::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x)
{   
    printf("Gstarx not yet defined for no arums case.\n");

    return 0.0;
}

double trame::arums::none::Gbar(const arma::mat& Ubar, const arma::mat& mubar, const arma::vec& n, arma::mat& U_out, arma::mat& mu_out)
{   
    double val=0.0, val_temp;
    
    U_out.set_size(nbX,nbY);
    mu_out.set_size(nbX,nbY);
    arma::mat U_x_temp, mu_x_temp;
    //
    for (int i=0; i<nbX; i++) {
        val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),U_x_temp,mu_x_temp);
        //
        val += n(i)*val_temp;
        U_out.row(i) = arma::trans(U_x_temp);
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double trame::arums::none::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out)
{
    int count_int=0;
    int nbY0 = Ubar_x.n_elem;
    //
    //arma::mat srt = arma::sort(Ubar_x,"descend");
    arma::uvec srt_ind = arma::sort_index(Ubar_x,"descend");
    //
    mu_x_out.set_size(nbY0,1);
    double cumul = arma::as_scalar(mubar_x(srt_ind(count_int)));
    //
    while ((count_int < nbY0-1) && (cumul < 1.0) && (Ubar_x(srt_ind(count_int)) > 0)) {
        mu_x_out(srt_ind(count_int)) = mubar_x(srt_ind(count_int));
        count_int++;
        cumul += mubar_x(srt_ind(count_int)); // Keith: is this in the correct place?
    }
    //
    if (Ubar_x(srt_ind(count_int)) > 0) {
        mu_x_out(srt_ind(count_int)) = mubar_x(srt_ind(count_int)) + 1 - cumul;
    }
    //
    U_x_out = arma::zeros(nbY0,1);
    double valx = arma::accu(mu_x_out % Ubar_x);
    //
    return valx;
}

// just to conform with other arums classes
double trame::arums::none::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, int x)
{
    double val_x = this->Gbarx(Ubar_x, mubar_x, U_x_out, mu_x_out);
    //
    return val_x;
}

arma::vec trame::arums::none::dtheta_NablaGstar()
{
    arma::vec ret = arma::zeros(nbX*nbY,1);
    return ret;
}

trame::arums::empirical trame::arums::none::simul()
{
    empirical emp_obj;
    
    this->simul(emp_obj,NULL,NULL);
    //
    return emp_obj;
}

trame::arums::empirical trame::arums::none::simul(int* nbDraws, int* seed)
{
    empirical emp_obj;
    
    this->simul(emp_obj,nbDraws,seed);
    //
    return emp_obj;
}

void trame::arums::none::simul(empirical& obj_out)
{
    this->simul(obj_out,NULL,NULL);
}

void trame::arums::none::simul(empirical& obj_out, int* nbDraws, int* seed_val)
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
    arma::cube atoms(n_draws,nbY+1,nbX);
    atoms.fill(0);
    //
    obj_out.nbX = nbX;
    obj_out.nbY = nbY;
    obj_out.nbParams = atoms.n_elem;
    obj_out.atoms = atoms;
    obj_out.aux_nbDraws = n_draws;
    obj_out.xHomogenous = false;
    //
    if (seed_val) {
        arma::arma_rng::set_seed_random(); // need to reset the seed
    }
}

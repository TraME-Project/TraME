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
 * none class
 *
 * Keith O'Hara
 * 08/08/2016
 */

#include "trame.hpp"

trame::none::none(int nbX_inp, int nbY_inp)
{   
    this->build(nbX_inp,nbY_inp);
}

void trame::none::build(int nbX_inp, int nbY_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
    nbParams = 0;
}

double trame::none::G(arma::vec n)
{   
    double val = this->G(n,U,mu_sol);
    //
    return val;
}

double trame::none::G(arma::vec n, const arma::mat& U_inp, arma::mat& mu_out)
{   
    int i;
    double val=0.0, val_x;
    
    mu_out.set_size(nbX,nbY);
    arma::mat mu_x_temp;
    //
    for (i=0; i<nbX; i++) {
        val_x = Gx(U_inp.row(i).t(),mu_x_temp);
        //
        val += n(i)*val_x;
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double trame::none::Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out)
{
    arma::uvec temp_vec = which_max(U_x_inp, (int) 0);
    int y = temp_vec(0);
    //
    mu_x_out.set_size(nbY,1);
    mu_x_out.zeros();
    
    if (y < nbY) {
        mu_x_out(y) = 1;
    }
    //
    double val_x = std::max(elem_max(U_x_inp), (double) 0.0);
    //
    return val_x;
}

// just to conform with other arums classes
double trame::none::Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, int x)
{
    double val_x = 0.0;

    if (U_x_inp.n_rows > 1 && U_x_inp.n_cols > 1) {
        val_x = this->Gx(U_x_inp.row(x).t(),mu_x_out);
    } if (U_x_inp.n_rows == 1 && U_x_inp.n_cols > 1) { 
        val_x = this->Gx(U_x_inp.t(),mu_x_out);
    } else {
        val_x = this->Gx(U_x_inp,mu_x_out);
    }
    //
    return val_x;
}

double trame::none::Gstar(arma::vec n, const arma::mat& mu_inp, arma::mat& U_out)
{   
    printf("Gstar not yet defined for no arums case.\n");

    return 0.0;
}

double trame::none::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x)
{   
    printf("Gstarx not yet defined for no arums case.\n");

    return 0.0;
}

double trame::none::Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_out, arma::mat& mu_out)
{   
    int i;
    double val=0.0, val_temp;
    
    U_out.set_size(nbX,nbY);
    mu_out.set_size(nbX,nbY);
    arma::mat U_x_temp, mu_x_temp;
    //
    for (i=0; i<nbX; i++) {
        val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),U_x_temp,mu_x_temp);
        //
        val += n(i)*val_temp;
        U_out.row(i) = arma::trans(U_x_temp);
        mu_out.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double trame::none::Gbarx(arma::mat Ubarx, arma::mat mubarx, arma::mat& U_x_out, arma::mat& mu_x_out)
{
    int count_int=0;
    int nbY0 = Ubarx.n_elem;
    //
    //arma::mat srt = arma::sort(Ubarx,"descend");
    arma::uvec srt_ind = arma::sort_index(Ubarx,"descend");
    //
    mu_x_out.set_size(nbY0,1);
    double cumul = arma::as_scalar(mubarx(srt_ind(count_int)));
    //
    while ((count_int < nbY0-1) & (cumul < 1.0) & (Ubarx(srt_ind(count_int)) > 0)) {
        mu_x_out(srt_ind(count_int)) = mubarx(srt_ind(count_int));
        count_int++;
        cumul += mubarx(srt_ind(count_int)); // Keith: is this in the correct place?
    }
    //
    if (Ubarx(srt_ind(count_int)) > 0) {
        mu_x_out(srt_ind(count_int)) = mubarx(srt_ind(count_int)) + 1 - cumul;
    }
    //
    U_x_out = arma::zeros(nbY0,1);
    //
    double valx = arma::accu(mu_x_out % Ubarx);
    //
    return valx;
}

arma::vec trame::none::dtheta_NablaGstar()
{
    arma::vec ret = arma::zeros(nbX*nbY,1);
    return ret;
}

void trame::none::simul(empirical &ret, int nbDraws, int seed_val)
{
    arma::arma_rng::set_seed(seed_val);
    //
    arma::cube atoms(nbDraws,nbY+1,nbX);
    atoms.fill(0);
    //
    ret.nbX = nbX;
    ret.nbY = nbY;
    ret.nbParams = atoms.n_elem;
    ret.atoms = atoms;
    ret.aux_nbDraws = nbDraws;
    ret.xHomogenous = false;
    //
    arma::arma_rng::set_seed_random(); // need to reset the seed
}

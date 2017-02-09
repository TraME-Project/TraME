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
 * logit class
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 02/09/2017
 */

#include "trame.hpp"

trame::logit::logit(int nbX_inp, int nbY_inp)
{   
    this->build(nbX_inp,nbY_inp, (double) 1.0, (bool) true);
}

trame::logit::logit(int nbX_inp, int nbY_inp, double sigma_inp, bool outsideOption_inp)
{   
    this->build(nbX_inp,nbY_inp,sigma_inp,outsideOption_inp);
}

void trame::logit::build(int nbX_inp, int nbY_inp)
{   
    this->build(nbX_inp,nbY_inp, (double) 1.0, (bool) true);
}

void trame::logit::build(int nbX_inp, int nbY_inp, double sigma_inp, bool outsideOption_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
    nbParams = 1;
    sigma = sigma_inp;
    outsideOption = outsideOption_inp;
}

double trame::logit::G(const arma::vec& n)
{   
    double val = this->G(n,U,mu_sol);
    //
    return val;
}

double trame::logit::G(const arma::vec& n, const arma::mat& U_inp, arma::mat& mu_out)
{   
    arma::mat denom;
    
    arma::mat expU = arma::exp(U_inp / sigma);
    //
    if (outsideOption) {
        denom = 1 + arma::sum(expU,1);
    } else {
        denom = arma::sum(expU,1);
    }
    //
    double val = sigma*arma::accu(n % arma::log(denom));
    mu_out = elem_prod(n/denom, expU);
    //
    return val;
}

double trame::logit::Gstar(const arma::vec& n)
{
    double val = this->Gstar(n,mu_sol,U_sol);
    //
    return val;
}

double trame::logit::Gstar(const arma::vec& n, const arma::mat& mu_inp, arma::mat& U_out)
{
    double val = 0.0;

    arma::mat n_repd = arma::repmat(n,1,nbY);
    //
    if (outsideOption) {
        arma::mat mu_x_0 = n - arma::sum(mu_inp,1);
        
        val   = sigma * ( arma::accu(mu_inp % arma::log(mu_inp/n_repd)) + arma::accu(mu_x_0 % arma::log(mu_x_0/n)) );
        U_out = sigma * arma::log(elem_div(mu_inp, mu_x_0));
    } else {
        val   = sigma * arma::accu(mu_inp % arma::log(mu_inp/n_repd));
        U_out = sigma * arma::log(mu_inp / n_repd);
    }
    //
    return val;
}

double trame::logit::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out)
{
    double val_x = 0.0;
    //
    if (outsideOption) {
        double mu0 = 1 - arma::accu(mu_x_inp);
        
        val_x   = sigma * ( mu0 * std::log(mu0) + arma::accu(mu_x_inp % arma::log(mu_x_inp)) );
        U_x_out = sigma * arma::log(mu_x_inp / mu0);
    } else {
        val_x   = sigma * arma::accu(mu_x_inp % arma::log(mu_x_inp));
        U_x_out = sigma * arma::log(mu_x_inp);
    }
    //
    return val_x;
}

// just to conform with other arums classes
double trame::logit::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x)
{
    double val_x = 0.0;

    if (mu_x_inp.n_rows > 1 && mu_x_inp.n_cols > 1) {
        val_x = this->Gstarx(mu_x_inp.row(x).t(),U_x_out);
    } if (mu_x_inp.n_rows == 1 && mu_x_inp.n_cols > 1) { 
        val_x = this->Gstarx(mu_x_inp.t(),U_x_out);
    } else {
        val_x = this->Gstarx(mu_x_inp,U_x_out);
    }
    //
    return val_x;
}

double trame::logit::Gbar(const arma::mat& Ubar, const arma::mat& mubar, const arma::vec& n, arma::mat& U_out, arma::mat& mu_out)
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

double differMargX(double z, const trame::trame_zeroin_data& opt_data)
{
    arma::mat temp_mat = arma::min(z * opt_data.exp_Ubar_X, opt_data.mubar_X);
    double ret = z + arma::accu(temp_mat) - 1;
    //
    return ret;
}

double trame::logit::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out)
{
    double valx = 0.0;

    if (outsideOption) {
        arma::mat exp_Ubar_X = arma::exp(Ubar_x/sigma);
        
        trame_zeroin_data root_data;
        root_data.exp_Ubar_X = exp_Ubar_X;
        root_data.mubar_X = mubar_x;
        
        double mu_x_0 = zeroin(0.0, 1.0, differMargX, root_data, NULL, NULL);
        //
        mu_x_out = arma::min(mu_x_0 * exp_Ubar_X, mubar_x);
        U_x_out  = sigma * arma::log(mu_x_out/mu_x_0);
        //
        valx = arma::accu(mu_x_out % Ubar_x) - sigma*(mu_x_0*std::log(mu_x_0) + arma::accu(mu_x_out % arma::log(mu_x_out)));
    }
    //
    return valx;
}

// just to conform with other arums classes
double trame::logit::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, int x)
{
    double val_x = this->Gbarx(Ubar_x, mubar_x, U_x_out, mu_x_out);
    //
    return val_x;
}

arma::mat trame::logit::D2G(const arma::vec& n, bool xFirst)
{
    arma::mat ret;

    this->D2G(ret,n,xFirst);
    //
    return ret;
}

void trame::logit::D2G(arma::mat& H, const arma::vec& n, bool xFirst)
{
    // NOTE: the formula is the same regardless of whether outsideOption is true
    int x, y, yprime;
    
    arma::mat mu_xy = mu_sol;
    H.zeros(nbX*nbY,nbX*nbY);
    //
    for (x = 0; x < nbX; x++) {
        for (y = 0; y < nbY; y++) {
            for (yprime = 0; yprime < nbY; yprime++) {
                if (xFirst) {
                    if (y==yprime) {
                        H(x + nbX*y, x + nbX*yprime) = mu_xy(x,y)*(1 - mu_xy(x,y)/n(x)) / sigma;
                    } else {
                        H(x + nbX*y, x + nbX*yprime) = - mu_xy(x,y)*mu_xy(x,yprime) / (n(x)*sigma);
                    }
                } else {
                    if (y==yprime) {
                        H(y + nbY*x, yprime + nbY*x) = mu_xy(x,y)*(1 - mu_xy(x,y)/n(x)) / sigma;
                    } else {
                        H(y + nbY*x, yprime + nbY*x) = - mu_xy(x,y)*mu_xy(x,yprime) / (n(x)*sigma);
                    }
                }
            }
        }
    }
}

void trame::logit::D2G(arma::mat& H, const arma::mat& U_inp, const arma::vec& n, bool xFirst)
{
    // NOTE: the formula is the same regardless of whether outsideOption is true
    int x, y, yprime;

    arma::mat mu_xy;
    this->G(n,U_inp,mu_xy);
    
    H.zeros(nbX*nbY,nbX*nbY);
    //
    for (x = 0; x < nbX; x++) {
        for (y = 0; y < nbY; y++) {
            for (yprime = 0; yprime < nbY; yprime++) {
                if (xFirst) {
                    if (y==yprime) {
                        H(x + nbX*y, x + nbX*yprime) = mu_xy(x,y)*(1 - mu_xy(x,y)/n(x)) / sigma;
                    } else {
                        H(x + nbX*y, x + nbX*yprime) = - mu_xy(x,y)*mu_xy(x,yprime) / (n(x)*sigma);
                    }
                } else {
                    if (y==yprime) {
                        H(y + nbY*x, yprime + nbY*x) = mu_xy(x,y)*(1 - mu_xy(x,y)/n(x)) / sigma;
                    } else {
                        H(y + nbY*x, yprime + nbY*x) = - mu_xy(x,y)*mu_xy(x,yprime) / (n(x)*sigma);
                    }
                }
            }
        }
    }
}

arma::mat trame::logit::D2Gstar(const arma::vec& n, bool xFirst)
{
    arma::mat ret;

    this->D2Gstar(ret,n,xFirst);
    //
    return ret;
}

void trame::logit::D2Gstar(arma::mat &H, const arma::vec& n, bool xFirst)
{
    // NOTE: the formula is the same regardless of whether outsideOption == 1 or 0
    int x, y, yprime;
    
    arma::mat mu_x_0 = n - arma::sum(mu_sol,1);
    arma::mat mu_x_0_recip = arma::ones(mu_x_0.n_rows,mu_x_0.n_cols) / mu_x_0; // reciprocal of mu_x_0
    arma::mat mu_xy_recip = arma::ones(mu_sol.n_rows,mu_sol.n_cols) / mu_sol;
    
    H.zeros(nbX*nbY,nbX*nbY);
    //
    for (x = 0; x < nbX; x++) {
        for (y = 0; y < nbY; y++) {
            for (yprime = 0; yprime < nbY; yprime++) {
                if (xFirst) {
                    if (y==yprime) {
                        H(x + nbX*y, x + nbX*yprime) = mu_x_0_recip(x) + mu_xy_recip(x,y);
                    }else{
                        H(x + nbX*y, x + nbX*yprime) = mu_x_0_recip(x);
                    }
                } else {
                    if (y==yprime) {
                        H(y + nbY*x, yprime + nbY*x) = mu_x_0_recip(x) + mu_xy_recip(x,y);
                    } else {
                        H(y + nbY*x, yprime + nbY*x) = mu_x_0_recip(x);
                    }
                }
            }
        }
    }
    H *= sigma;
}

arma::mat trame::logit::dtheta_NablaGstar(const arma::vec& n, arma::mat* dtheta_inp, bool xFirst)
{
    arma::mat ret;

    this->dtheta_NablaGstar(ret,n,dtheta_inp,xFirst);
    //
    return ret;
}

void trame::logit::dtheta_NablaGstar(arma::mat &ret, const arma::vec& n, arma::mat* dtheta_inp, bool xFirst)
{
    arma::mat logmu_temp, mu_x_0;

    arma::mat dtheta = (dtheta_inp) ? *dtheta_inp : arma::ones(1,1);
    
    if (dtheta.n_elem==0) {
        ret.zeros(nbX*nbY,0);
    } else {
        if (outsideOption) {
            mu_x_0 = arma::repmat(n - arma::sum(mu,1),1,mu.n_cols);
            
            if (xFirst) {
                logmu_temp = arma::vectorise(arma::log(mu/mu_x_0));
            } else {
                logmu_temp = arma::vectorise(arma::trans(arma::log(mu/mu_x_0)));
            }
            //
            ret = arma::vectorise(dtheta) % logmu_temp;
        } else {
            if (xFirst) {
                logmu_temp = arma::vectorise(log(mu));
            } else {
                logmu_temp = arma::vectorise(arma::trans(arma::log(mu)));
            }
            //
            ret = arma::vectorise(dtheta) % logmu_temp;
        }
    }
}

trame::empirical trame::logit::simul()
{
    empirical emp_obj;
    
    this->simul(emp_obj,NULL,NULL);
    //
    return emp_obj;
}

trame::empirical trame::logit::simul(int* nbDraws, int* seed)
{
    empirical emp_obj;
    
    this->simul(emp_obj,nbDraws,seed);
    //
    return emp_obj;
}

void trame::logit::simul(empirical& obj_out)
{
    this->simul(obj_out,NULL,NULL);
}

void trame::logit::simul(empirical& obj_out, int* nbDraws, int* seed_val)
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
    // note: digamma(1) \approx -0.5772156649
    arma::cube atoms;
    if (outsideOption) {
        atoms = -0.5772156649 - sigma * arma::log( - arma::log(arma::randu(n_draws,nbY+1,nbX)) );
    } else {
        atoms = -0.5772156649 - sigma * arma::log( - arma::log(arma::randu(n_draws,nbY,nbX)) );
    }
    //
    obj_out.nbX = nbX;
    obj_out.nbY = nbY;

    obj_out.nbParams = atoms.n_elem;
    obj_out.atoms = atoms;
    obj_out.aux_nbDraws = n_draws;

    obj_out.xHomogenous = false;
    obj_out.outsideOption = outsideOption;

    if (outsideOption) {
        obj_out.nbOptions = nbY + 1;
    } else {
        obj_out.nbOptions = nbY;
    }
    //
    if (seed_val) {
        arma::arma_rng::set_seed_random(); // need to reset the seed
    }
}

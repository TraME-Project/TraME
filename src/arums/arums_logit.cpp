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
 * logit additive random utility model (ARUM) class
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 02/21/2017
 */

#include "trame.hpp"

trame::arums::logit::logit(int nbX_inp, int nbY_inp)
{   
    this->build(nbX_inp,nbY_inp);
}

trame::arums::logit::logit(int nbX_inp, int nbY_inp, double sigma_inp, bool outsideOption_inp)
{   
    this->build(nbX_inp,nbY_inp,sigma_inp,outsideOption_inp);
}

void 
trame::arums::logit::build(int nbX_inp, int nbY_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
    nbParams = 1;
}

void 
trame::arums::logit::build(int nbX_inp, int nbY_inp, double sigma_inp, bool outsideOption_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
    nbParams = 1;
    sigma = sigma_inp;
    outsideOption = outsideOption_inp;
}

double 
trame::arums::logit::G(const arma::vec& n)
{   
    return this->G(n,U,mu_sol);
}

double 
trame::arums::logit::G(const arma::vec& n, const arma::mat& U_inp, arma::mat& mu_out)
{   
    arma::mat expU = arma::exp(U_inp / sigma);
    arma::vec denom = (outsideOption) ? (1.0 + arma::sum(expU,1)) : 0.0 + arma::sum(expU,1); // the '0.0 + ' fixes a compiling bug
    //
    mu_out = elem_prod(n/denom, expU);
    double val = sigma*arma::accu(n % arma::log(denom));
    //
    return val;
}

double 
trame::arums::logit::Gstar(const arma::vec& n)
{
    return this->Gstar(n,mu_sol,U_sol);
}

double 
trame::arums::logit::Gstar(const arma::vec& n, const arma::mat& mu_inp, arma::mat& U_out)
{
    double val = 0.0;

    arma::mat n_repd = arma::repmat(n,1,nbY);
    //
    if (outsideOption) {
        arma::vec mu_x0 = n - arma::sum(mu_inp,1);
        
        val   = sigma * ( arma::accu(mu_inp % arma::log(mu_inp/n_repd)) + arma::accu(mu_x0 % arma::log(mu_x0/n)) );
        U_out = sigma * arma::log(elem_div(mu_inp, mu_x0));
    } else {
        val   = sigma * arma::accu(mu_inp % arma::log(mu_inp/n_repd));
        U_out = sigma * arma::log(mu_inp / n_repd);
    }
    //
    return val;
}

double 
trame::arums::logit::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out)
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
double 
trame::arums::logit::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x)
{
    double val_x = 0.0;

    if (mu_x_inp.n_rows > 1 && mu_x_inp.n_cols > 1) {
        val_x = this->Gstarx(mu_x_inp.row(x).t(),U_x_out);
    } else if (mu_x_inp.n_rows == 1 && mu_x_inp.n_cols > 1) { 
        val_x = this->Gstarx(mu_x_inp.t(),U_x_out);
        // val_x = this->Gstarx(mu_x_inp,U_x_out);
    } else {
        val_x = this->Gstarx(mu_x_inp,U_x_out);
    }
    //
    return val_x;
}

double 
trame::arums::logit::Gbar(const arma::mat& Ubar, const arma::mat& mubar, const arma::vec& n, arma::mat& U_out, arma::mat& mu_out)
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

double 
trame::arums::logit::differMargX(double z, void* opt_data)
{
    trame_logit_zeroin_data *d = reinterpret_cast<trame_logit_zeroin_data*>(opt_data);

    arma::mat temp_mat = arma::min(z * d->exp_Ubar_X, d->mubar_X);
    double ret = z + arma::accu(temp_mat) - 1;
    //
    return ret;
}

double 
trame::arums::logit::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out)
{
    double valx = 0.0;

    if (outsideOption) {
        arma::mat exp_Ubar_X = arma::exp(Ubar_x/sigma);
        
        trame_logit_zeroin_data root_data;
        root_data.exp_Ubar_X = exp_Ubar_X;
        root_data.mubar_X = mubar_x;
        
        double mu_x0 = zeroin(0.0, 1.0, differMargX, &root_data, NULL, NULL);
        //
        mu_x_out = arma::min(mu_x0 * exp_Ubar_X, mubar_x);
        U_x_out  = sigma * arma::log(mu_x_out/mu_x0);
        //
        valx = arma::accu(mu_x_out % Ubar_x) - sigma*(mu_x0*std::log(mu_x0) + arma::accu(mu_x_out % arma::log(mu_x_out)));
    }
    //
    return valx;
}

// just to conform with other arums classes
double 
trame::arums::logit::Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, int x)
{
    return this->Gbarx(Ubar_x, mubar_x, U_x_out, mu_x_out);
}

arma::mat 
trame::arums::logit::D2G(const arma::vec& n, bool x_first)
{
    arma::mat ret;
    this->D2G(ret,n,U,x_first);
    //
    return ret;
}

void 
trame::arums::logit::D2G(arma::mat& H, const arma::vec& n, bool x_first)
{
    // NOTE: the formula is the same regardless of whether outsideOption is true
    
    arma::mat mu_xy = mu_sol;
    H.zeros(nbX*nbY,nbX*nbY);
    //
    for (int x = 0; x < nbX; x++) {
        for (int y = 0; y < nbY; y++) {
            for (int y_prime = 0; y_prime < nbY; y_prime++) {
                if (x_first) {
                    if (y==y_prime) {
                        H(x + nbX*y, x + nbX*y_prime) = mu_xy(x,y)*(1 - mu_xy(x,y)/n(x)) / sigma;
                    } else {
                        H(x + nbX*y, x + nbX*y_prime) = - mu_xy(x,y)*mu_xy(x,y_prime) / (n(x)*sigma);
                    }
                } else {
                    if (y==y_prime) {
                        H(y + nbY*x, y_prime + nbY*x) = mu_xy(x,y)*(1 - mu_xy(x,y)/n(x)) / sigma;
                    } else {
                        H(y + nbY*x, y_prime + nbY*x) = - mu_xy(x,y)*mu_xy(x,y_prime) / (n(x)*sigma);
                    }
                }
            }
        }
    }
}

arma::mat 
trame::arums::logit::D2G(const arma::vec& n, const arma::mat& U_inp, bool x_first)
{
    arma::mat ret;
    this->D2G(ret,n,U_inp,x_first);
    //
    return ret;
}

void 
trame::arums::logit::D2G(arma::mat& H, const arma::vec& n, const arma::mat& U_inp, bool x_first)
{
    // NOTE: the formula is the same regardless of whether outsideOption is true

    arma::mat mu_xy;
    this->G(n,U_inp,mu_xy);
    
    H.zeros(nbX*nbY,nbX*nbY);
    //
    for (int x = 0; x < nbX; x++) {
        for (int y = 0; y < nbY; y++) {
            for (int y_prime = 0; y_prime < nbY; y_prime++) {
                if (x_first) {
                    if (y==y_prime) {
                        H(x + nbX*y, x + nbX*y_prime) = mu_xy(x,y)*(1 - mu_xy(x,y)/n(x)) / sigma;
                    } else {
                        H(x + nbX*y, x + nbX*y_prime) = - mu_xy(x,y)*mu_xy(x,y_prime) / (n(x)*sigma);
                    }
                } else {
                    if (y==y_prime) {
                        H(y + nbY*x, y_prime + nbY*x) = mu_xy(x,y)*(1 - mu_xy(x,y)/n(x)) / sigma;
                    } else {
                        H(y + nbY*x, y_prime + nbY*x) = - mu_xy(x,y)*mu_xy(x,y_prime) / (n(x)*sigma);
                    }
                }
            }
        }
    }
}

arma::mat 
trame::arums::logit::D2Gstar(const arma::vec& n, bool x_first)
{
    arma::mat ret;
    this->D2Gstar(ret,n,mu,x_first);
    //
    return ret;
}

void 
trame::arums::logit::D2Gstar(arma::mat &ret, const arma::vec& n, bool x_first)
{
    this->D2Gstar(ret,n,mu,x_first);
}

arma::mat 
trame::arums::logit::D2Gstar(const arma::vec& n, const arma::mat& mu_inp, bool x_first)
{
    arma::mat ret;
    this->D2Gstar(ret,n,mu_inp,x_first);
    //
    return ret;
}

void 
trame::arums::logit::D2Gstar(arma::mat &H, const arma::vec& n, const arma::mat& mu_inp, bool x_first)
{
    // NOTE: the formula is the same regardless of whether outsideOption == 1 or 0

    arma::vec mu_x0 = n - arma::sum(mu_inp,1);
    arma::vec mu_x0_recip = arma::ones(mu_x0.n_rows,1) / mu_x0; // reciprocal of mu_x0
    arma::mat mu_xy_recip = arma::ones(mu_inp.n_rows,mu_inp.n_cols) / mu_inp;
    
    H.zeros(nbX*nbY,nbX*nbY);
    //
    for (int x = 0; x < nbX; x++) {
        for (int y = 0; y < nbY; y++) {
            for (int y_prime = 0; y_prime < nbY; y_prime++) {
                if (x_first) {
                    if (y==y_prime) {
                        H(x + nbX*y, x + nbX*y_prime) = mu_x0_recip(x) + mu_xy_recip(x,y);
                    }else{
                        H(x + nbX*y, x + nbX*y_prime) = mu_x0_recip(x);
                    }
                } else {
                    if (y==y_prime) {
                        H(y + nbY*x, y_prime + nbY*x) = mu_x0_recip(x) + mu_xy_recip(x,y);
                    } else {
                        H(y + nbY*x, y_prime + nbY*x) = mu_x0_recip(x);
                    }
                }
            }
        }
    }
    H *= sigma;
}

arma::mat 
trame::arums::logit::dparams_NablaGstar(const arma::vec& n, arma::mat* dparams_inp, bool x_first)
{
    arma::mat ret;
    this->dparams_NablaGstar(ret,n,mu_sol,dparams_inp,x_first);
    //
    return ret;
}

void 
trame::arums::logit::dparams_NablaGstar(arma::mat &ret, const arma::vec& n, arma::mat* dparams_inp, bool x_first)
{
    this->dparams_NablaGstar(ret,n,mu_sol,dparams_inp,x_first);
}

arma::mat 
trame::arums::logit::dparams_NablaGstar(const arma::vec& n, const arma::mat& mu_inp, arma::mat* dparams_inp, bool x_first)
{
    arma::mat ret;
    this->dparams_NablaGstar(ret,n,mu_inp,dparams_inp,x_first);
    //
    return ret;
}

void 
trame::arums::logit::dparams_NablaGstar(arma::mat &ret, const arma::vec& n, const arma::mat& mu_inp, arma::mat* dparams_inp, bool x_first)
{
    arma::mat logmu_temp, mu_x0;

    arma::mat dparams = (dparams_inp) ? *dparams_inp : arma::ones(1,1);
    
    if (dparams.n_elem==0) {
        ret.zeros(nbX*nbY,0);
    } else {
        if (outsideOption) {
            mu_x0 = arma::repmat(n - arma::sum(mu_inp,1),1,mu_inp.n_cols);
            
            if (x_first) {
                logmu_temp = arma::vectorise(arma::log(mu_inp/mu_x0));
            } else {
                logmu_temp = arma::vectorise(arma::trans(arma::log(mu_inp/mu_x0)));
            }
            //
            ret = arma::vectorise(dparams) % logmu_temp;
        } else {
            if (x_first) {
                logmu_temp = arma::vectorise(log(mu_inp));
            } else {
                logmu_temp = arma::vectorise(arma::trans(arma::log(mu_inp)));
            }
            //
            ret = arma::vectorise(dparams) % logmu_temp;
        }
    }
}

trame::arums::empirical 
trame::arums::logit::simul()
{
    empirical emp_obj;
    
    this->simul(emp_obj,NULL,NULL);
    //
    return emp_obj;
}

trame::arums::empirical 
trame::arums::logit::simul(int* nbDraws, int* seed)
{
    empirical emp_obj;
    
    this->simul(emp_obj,nbDraws,seed);
    //
    return emp_obj;
}

void 
trame::arums::logit::simul(empirical& obj_out)
{
    this->simul(obj_out,NULL,NULL);
}

void 
trame::arums::logit::simul(empirical& obj_out, int* nbDraws, int* seed_val)
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

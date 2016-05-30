/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 Alfred Galichon and the TraME Team
  ##
  ##   This file is part of the R package TraME.
  ##
  ##   The R package TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

#include "zeroin.hpp"

// logit class
class logit
{
    public:
        // build_logit objects
        int nbX;
        int nbY;
        int nbParams;
        double sigma;
        bool outsideOption;
        
        // equilibrium objects
        arma::mat mu;
        arma::mat U;
        
        arma::mat U_sol;
        arma::mat mu_sol;
        
        // member functions
        void build(int nbX_b, int nbY_b, int nbParams_b, double sigma_b, bool outsideOption_b);
        double G(arma::vec n);
        double G(arma::mat &mu_ret, arma::vec n);
        
        double Gstar(arma::vec n);
        double Gstar(arma::mat &U_ret, arma::vec n);
        double Gstarx(arma::mat &Ux_ret, arma::mat mux);
        
        double Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_inp, arma::mat& mu_inp);
        double Gbarx(arma::mat Ubarx, arma::mat mubarx, arma::mat& Ux_inp, arma::mat& mux_inp);
        
        void D2G(arma::mat &H, arma::vec n, bool xFirst);
        void D2Gstar(arma::mat &H, arma::vec n, bool xFirst);
        void dtheta_NablaGstar(arma::mat &ret, arma::vec n, arma::mat dtheta, bool xFirst);
        
        void simul(empirical &ret, int nbDraws, int seed);
        
    //private:
};

// for convenience:
void logit::build(int nbX_b, int nbY_b, int nbParams_b, double sigma_b, bool outsideOption_b)
{   
    nbX = nbX_b;
    nbY = nbY_b;
    nbParams = nbParams_b;
    sigma = sigma_b;
    outsideOption = outsideOption_b;
}

double logit::G(arma::vec n)
{   
    arma::mat denom;
    
    arma::mat expU = arma::exp(U / sigma);
    //
    if (outsideOption) {
        denom = 1 + arma::sum(expU,1);
    } else {
        denom = arma::sum(expU,1);
    }
    //
    double val = sigma*arma::accu(n % arma::log(denom));
    mu_sol = arma::repmat(n/denom,1,nbY) % expU;
    //
    return val;
}

double logit::G(arma::mat &mu_ret, arma::vec n)
{   
    arma::mat denom;
    
    arma::mat expU = arma::exp(U / sigma);
    //
    if (outsideOption) {
        denom = 1 + arma::sum(expU,1);
    } else {
        denom = arma::sum(expU,1);
    }
    //
    double val = sigma*arma::accu(n % arma::log(denom));
    mu_ret = arma::repmat(n/denom,1,nbY) % expU;
    //
    return val;
}

double logit::Gstar(arma::vec n)
{
    double val=0.0;
    
    arma::mat mu_x_0;
    arma::mat n_repd = arma::repmat(n,1,nbY);
    //
    if(outsideOption){
        mu_x_0 = n - arma::sum(mu_sol,1);
        
        val = sigma * ( arma::accu(mu_sol % arma::log(mu_sol/n_repd)) + arma::accu(mu_x_0 % arma::log(mu_x_0/n)) );
        U_sol   = sigma * arma::log(mu_sol / arma::repmat(mu_x_0,1,nbY));
    }else{
        val = sigma * arma::accu(mu_sol % arma::log(mu_sol/n_repd));
        U_sol   = sigma * arma::log(mu_sol / n_repd);
    }
    //
    return val;
}

double logit::Gstar(arma::mat &U_ret, arma::vec n)
{
    double val=0.0;
    
    arma::mat mu_x_0;
    arma::mat n_repd = arma::repmat(n,1,nbY);
    //
    if(outsideOption){
        mu_x_0 = n - arma::sum(mu_sol,1);
        
        val   = sigma * ( arma::accu(mu_sol % arma::log(mu_sol/n_repd)) + arma::accu(mu_x_0 % arma::log(mu_x_0/n)) );
        U_ret = sigma * arma::log(mu_sol / arma::repmat(mu_x_0,1,nbY));
    }else{
        val   = sigma * arma::accu(mu_sol % arma::log(mu_sol/n_repd));
        U_ret = sigma * arma::log(mu_sol / n_repd);
    }
    //
    return val;
}

double logit::Gstarx(arma::mat &Ux, arma::mat mux)
{
    double val_x=0.0;
    
    arma::mat mu0;
    //
    if(outsideOption){
        mu0 = 1 - arma::accu(mux);
        
        val_x   = sigma * ( arma::accu(mu0 % arma::log(mu0)) + arma::accu(mux % arma::log(mux)) );
        Ux = sigma * arma::log(mux / mu0);
    }else{
        val_x   = sigma * arma::accu(mux % arma::log(mux));
        Ux = sigma * arma::log(mux);
    }
    //
    return val_x;
}

double logit::Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_inp, arma::mat& mu_inp)
{   
    int i;
    double val=0.0, val_temp;
    
    U_inp.set_size(nbX,nbY);
    mu_inp.set_size(nbX,nbY);
    arma::mat Ux_temp, mux_temp;
    //
    for(i=0; i<nbX; i++){
        val_temp = logit::Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),Ux_temp,mux_temp);
        //
        val += n(i)*val_temp;
        U_inp.row(i) = arma::trans(Ux_temp);
        mu_inp.row(i) = arma::trans(n(i)*mux_temp);
    }
    //
    return val;
}

double differMargX(double z, const trame_zeroin_data& opt_data)
{
    arma::mat temp_mat = arma::min(z * opt_data.expUbarX, opt_data.mubarX);
    double ret = z + arma::accu(temp_mat) - 1;
    //
    return ret;
}

double logit::Gbarx(arma::mat Ubarx, arma::mat mubarx, arma::mat& Ux_inp, arma::mat& mux_inp)
{
    double valx=0.0;
    if(outsideOption){
        double tol_zero = 1E-12;
        double max_iter = 10000;
        
        arma::mat expUbarX = arma::exp(Ubarx/sigma);
        
        trame_zeroin_data root_data;
        root_data.expUbarX = expUbarX;
        root_data.mubarX = mubarx;
        
        double mu_x_0 = zeroin(0.0, 1.0, differMargX, root_data, tol_zero, max_iter);
        //
        mux_inp = arma::min(mu_x_0 * expUbarX, mubarx);
        Ux_inp  = sigma * arma::log(mux_inp/mu_x_0);
        //
        valx = arma::accu(mux_inp % Ubarx) - sigma*(mu_x_0*std::log(mu_x_0) + arma::accu(mux_inp % arma::log(mux_inp)));
    }
    return valx;
}

void logit::D2G(arma::mat &H, arma::vec n, bool xFirst)
{
    // NOTE: the formula is the same regardless of whether outsideOption is true
    int x, y, yprime;
    
    arma::mat muxy = mu_sol;
    H.zeros(nbX*nbY,nbX*nbY);
    //
    for(x = 0; x < nbX; x++){
        for(y = 0; y < nbY; y++){
            for(yprime = 0; yprime < nbY; yprime++){
                if(xFirst){
                    if(y==yprime){
                        H(x + nbX*y, x + nbX*yprime) = muxy(x,y)*(1 - muxy(x,y)/n(x)) / sigma;
                    }else{
                        H(x + nbX*y, x + nbX*yprime) = - muxy(x,y)*muxy(x,yprime) / (n(x)*sigma);
                    }
                }else{
                    if(y==yprime){
                        H(y + nbY*x, yprime + nbY*x) = muxy(x,y)*(1 - muxy(x,y)/n(x)) / sigma;
                    }else{
                        H(y + nbY*x, yprime + nbY*x) = - muxy(x,y)*muxy(x,yprime) / (n(x)*sigma);
                    }
                }
            }
        }
    }
}

void logit::D2Gstar(arma::mat &H, arma::vec n, bool xFirst)
{
    // NOTE: the formula is the same regardless of whether outsideOption == 1 or 0
    int x, y, yprime;
    
    arma::mat mu_x_0 = n - arma::sum(mu_sol,1);
    arma::mat mu_x_0_recip = arma::ones(mu_x_0.n_rows,mu_x_0.n_cols) / mu_x_0; // reciprocal of mu_x_0
    arma::mat muxy_recip = arma::ones(mu_sol.n_rows,mu_sol.n_cols) / mu_sol;
    
    H.zeros(nbX*nbY,nbX*nbY);
    //
    for(x = 0; x < nbX; x++){
        for(y = 0; y < nbY; y++){
            for(yprime = 0; yprime < nbY; yprime++){
                if(xFirst){
                    if(y==yprime){
                        H(x + nbX*y, x + nbX*yprime) = mu_x_0_recip(x) + muxy_recip(x,y);
                    }else{
                        H(x + nbX*y, x + nbX*yprime) = mu_x_0_recip(x);
                    }
                }else{
                    if(y==yprime){
                        H(y + nbY*x, yprime + nbY*x) = mu_x_0_recip(x) + muxy_recip(x,y);
                    }else{
                        H(y + nbY*x, yprime + nbY*x) = mu_x_0_recip(x);
                    }
                }
            }
        }
    }
    H *= sigma;
}

void logit::dtheta_NablaGstar(arma::mat &ret, arma::vec n, arma::mat dtheta, bool xFirst)
{
    arma::mat logmu_temp, mu_x_0;
    
    if(dtheta.n_elem==0){
        ret.zeros(nbX*nbY,0);
    }else{
        if(outsideOption){
            mu_x_0 = arma::repmat(n - arma::sum(mu,1),1,mu.n_cols);
            
            if(xFirst){
                logmu_temp = arma::vectorise(arma::log(mu/mu_x_0));
            }else{
                logmu_temp = arma::vectorise(arma::trans(arma::log(mu/mu_x_0)));
            }
            //
            ret = arma::vectorise(dtheta) % logmu_temp;
        }else{
            if(xFirst){
                logmu_temp = arma::vectorise(log(mu));
            }else{
                logmu_temp = arma::vectorise(arma::trans(arma::log(mu)));
            }
            //
            ret = arma::vectorise(dtheta) % logmu_temp;
        }
    }
}

void logit::simul(empirical &ret, int nbDraws, int seed_val)
{
    arma::arma_rng::set_seed(seed_val);
    
    // note: digamma(1) \approx -0.5772156649
    arma::cube atoms;
    if(outsideOption){
        atoms = -0.5772156649 - sigma * arma::log( - arma::log(arma::randu(nbDraws,nbY+1,nbX)) );
    }else{
        atoms = -0.5772156649 - sigma * arma::log( - arma::log(arma::randu(nbDraws,nbY,nbX)) );
    }
    //
    ret.nbX = nbX;
    ret.nbY = nbY;
    ret.nbParams = atoms.n_elem;
    ret.atoms = atoms;
    ret.aux_nbDraws = nbDraws;
    ret.xHomogenous = false;
    ret.outsideOption = outsideOption;
    //
    arma::arma_rng::set_seed_random(); // need to reset the seed
}

// this generated weird assembly code...
/*double differMargX_old(double z, void *opt_data)
{
    
    trame_opt_data *opt_data_temp = reinterpret_cast<trame_opt_data*>(opt_data);
    arma::mat* expUbarX = opt_data_temp->expUbarX;
    arma::mat* mubarX = opt_data_temp->mubarX;
    //
    arma::mat temp_mat = arma::min(z * (*expUbarX), (*mubarX));
    double ret = z + arma::accu(temp_mat) - 1;
    //
    return ret;
}*/

// this generated weird assembly code...
/*double logit::Gbarx_old(arma::mat Ubarx, arma::mat mubarx, arma::mat& Ux_inp, arma::mat& mux_inp)
{
    double valx=0.0;
    if(outsideOption){
        double tol_zero = 1E-12;
        double max_iter = 10000;
        
        arma::mat expUbarX = arma::exp(Ubarx/sigma);
        
        trame_opt_data root_data;
        root_data.expUbarX = &expUbarX;
        root_data.mubarX = &mubarx;
        
        double mu_x_0 = zeroin_old(0.0, 1.0, differMargX, &root_data, tol_zero, max_iter);
        //
        mux_inp = arma::min(mu_x_0 * expUbarX, mubarx);
        Ux_inp  = sigma * arma::log(mux/mu_x_0);
        //
        valx = arma::accu(mux % Ubarx) - sigma*(mu_x_0*std::log(mu_x_0) + arma::accu(mux % arma::log(mux)));
    }
    return valx;
}*/

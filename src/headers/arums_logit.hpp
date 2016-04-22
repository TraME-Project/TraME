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
        int outsideOption; // 1 = yes, 0 = no
        
        // equilibrium objects
        arma::mat mu;
        arma::mat mux;
        arma::mat U;
        arma::mat Ux;
        
        // member functions
        void build(int nbX_b, int nbY_b, int nbParams_b, double sigma_b, int outsideOption_b);
        void G(arma::mat &val, arma::vec n);
        void Gstar(arma::mat &val, arma::vec n);
        void Gstar(arma::mat &val, arma::mat &U_ret, arma::vec n);
        void Gstarx(arma::mat &valx, arma::mat mux);
        void Gstarx(arma::mat &valx, arma::mat &Ux_ret, arma::mat mux);
        void D2G(arma::mat &H, arma::vec n, int xFirst);
        void D2Gstar(arma::mat &H, arma::vec n, int xFirst);
        void dtheta_NablaGstar(arma::mat &ret, arma::mat mu, arma::vec n, arma::mat dtheta, int xFirst);
        void Gbarx(double &valx, arma::mat Ubarx, arma::mat mubarx);
        //void simul(arma::mat &ret, int nbDraws, int seed);
        
        
    //private:
};

// for convenience:
void logit::build(int nbX_b, int nbY_b, int nbParams_b, double sigma_b, int outsideOption_b)
{   
    nbX = nbX_b;
    nbY = nbY_b;
    nbParams = nbParams_b;
    sigma = sigma_b;
    outsideOption = outsideOption_b;
}

void logit::G(arma::mat &val, arma::vec n)
{   
    arma::mat denom;
    
    arma::mat expU = arma::exp(U / sigma);
    //
    if(outsideOption==1){
        denom = 1 + arma::sum(expU,1);
    }else{
        denom = arma::sum(expU,1);
    }
    //
    val = sigma*arma::sum(n * arma::log(denom));
    mu  = (double(n)/denom) % expU;
}

void logit::Gstar(arma::mat &val, arma::vec n)
{
    arma::mat mux0;
    //
    if(outsideOption==1){
        mux0 = n - arma::sum(mu,1);
        
        val = sigma * ( arma::accu(mu % arma::log(mu/n)) + arma::accu(mux0 % arma::log(mux0/n)) );
        U   = sigma * arma::log(mu / mux0);
    }else{
        val = sigma * arma::accu(mu % arma::log(mu/n));
        U   = sigma * arma::log(mu / n);
    }
}

void logit::Gstar(arma::mat &val, arma::mat &U_ret, arma::vec n)
{
    arma::mat mux0;
    //
    if(outsideOption==1){
        mux0 = n - arma::sum(mu,1);
        
        val   = sigma * ( arma::accu(mu % arma::log(mu/n)) + arma::accu(mux0 % arma::log(mux0/n)) );
        U_ret = sigma * arma::log(mu / mux0);
    }else{
        val   = sigma * arma::accu(mu % arma::log(mu/n));
        U_ret = sigma * arma::log(mu / n);
    }
}

void logit::Gstarx(arma::mat &valx, arma::mat mux)
{
    arma::mat mu0;
    //
    if(outsideOption==1){
        mu0 = 1 - arma::accu(mux);
        
        valx = sigma * ( arma::accu(mu0 % arma::log(mu0)) + arma::accu(mux % arma::log(mux)) );
        Ux   = sigma * arma::log(mux / mu0);
    }else{
        valx = sigma * arma::accu(mux % arma::log(mux));
        Ux   = sigma * arma::log(mux);
    }
}

void logit::Gstarx(arma::mat &valx, arma::mat &Ux_ret, arma::mat mux)
{
    arma::mat mu0;
    //
    if(outsideOption==1){
        mu0 = 1 - arma::accu(mux);
        
        valx   = sigma * ( arma::accu(mu0 % arma::log(mu0)) + arma::accu(mux % arma::log(mux)) );
        Ux_ret = sigma * arma::log(mux / mu0);
    }else{
        valx   = sigma * arma::accu(mux % arma::log(mux));
        Ux_ret = sigma * arma::log(mux);
    }
}

void logit::D2G(arma::mat &H, arma::vec n, int xFirst)
{
    // NOTE: the formula is the same regardless of whether outsideOption == 1 or 0
    int x, y, yprime;
    
    arma::mat muxy = mu;
    H.zeros(nbX*nbY,nbX*nbY);
    //
    for(x = 0; x < nbX; x++){
        for(y = 0; y < nbY; y++){
            for(yprime = 0; yprime < nbY; yprime++){
                if(xFirst==1){
                    if(y==yprime){
                        H(x + nbX*y, x + nbX*yprime) = muxy(x,y)*(1 - muxy(x,y)/n(x)) / sigma;
                    }else{
                        H(x + nbX*y, x + nbX*yprime) = - muxy(x,y)*muxy(x,yprime) / (n(x)*sigma);
                    }
                }else{
                    if(y==yprime){
                        H(y + nbY*x, yprime + nbY*x) = muxy(x,y)*(1 - muxy(x,y)/n(x)) / sigma;
                    }else{
                        H(y + nbY*x, yprime + nbY*x) = - muxy(x,y)*muxy(x,yprime) / (n(x)*arums$sigma);
                    }
                }
            }
        }
    }
}

void logit::D2Gstar(arma::mat &H, arma::vec n, int xFirst)
{
    // NOTE: the formula is the same regardless of whether outsideOption == 1 or 0
    int x, y, yprime;
    
    arma::mat mux0 = n - arma::sum(mu,1);
    arma::mat mux0_recip = arma::ones(mux0.n_rows,mux0.n_cols) / mux0; // reciprocal of mux0
    arma::mat muxy_recip = arma::ones(mu.n_rows,mu.n_cols) / mu;
    
    H.zeros(nbX*nbY,nbX*nbY);
    //
    for(x = 0; x < nbX; x++){
        for(y = 0; y < nbY; y++){
            for(yprime = 0; yprime < nbY; yprime++){
                if(xFirst==1){
                    if(y==yprime){
                        H(x + nbX*y, x + nbX*yprime) = mux0_recip(x) + muxy_recip(x,y);
                    }else{
                        H(x + nbX*y, x + nbX*yprime) = mux0_recip(x);
                    }
                }else{
                    if(y==yprime){
                        H(y + nbY*x, yprime + nbY*x) = mux0_recip(x) + muxy_recip(x,y);
                    }else{
                        H(y + nbY*x, yprime + nbY*x) = mux0_recip(x);
                    }
                }
            }
        }
    }
    H *= sigma;
}

void logit::dtheta_NablaGstar(arma::mat &ret, arma::mat mu, arma::vec n, arma::mat dtheta, int xFirst)
{
    arma::mat logmu_temp, mux0;
    
    if(dtheta.n_elem==0){
        ret.zeros(nbX*nbY,0); 
    }else{
        if(outsideOption==1){
            mux0 = arma::repmat(n - arma::sum(mu,1),1,mu.n_cols);
            
            if(xFirst==1){
                logmu_temp = arma::vectorise(arma::log(mu/mux0));
            }else{
                logmu_temp = arma::vectorise(arma::trans(arma::log(mu/mux0)));
            }
            //
            ret = arma::vectorise(dtheta) % logmu_temp;
        }else{
            if(xFirst==1){
                logmu_temp = arma::vectorise(log(mu));
            }else{
                logmu_temp = arma::vectorise(arma::trans(arma::log(mu)));
            }
            //
            ret = arma::vectorise(dtheta) % logmu_temp;
        }
    }
}

double differMargX(double z, trame_opt_data &opt_data)
{
    arma::mat expUbarX = opt_data.expUbarX;
    arma::mat mubarX = opt_data.mubarX;
    //
    arma::mat temp_mat = arma::min(z * expUbarX, mubarX);
    double ret = z + arma::accu(temp_mat) - 1;
    //
    return ret;
}

void logit::Gbarx(double &valx, arma::mat Ubarx, arma::mat mubarx)
{
    if(outsideOption==1){
        double tol_zero = 1E-08;
        double max_iter = 1000;
        
        arma::mat expUbarX = arma::exp(Ubarx/sigma);
        
        trame_opt_data root_data;
        root_data.expUbarX = expUbarX;
        root_data.mubarX = mubarx;
        
        double mux0 = zeroin(0.0, 1.0, differMargX, &root_data, tol_zero, max_iter);
        //
        mux = arma::min(mux0 * expUbarX, mubarx);
        Ux  = sigma * arma::log(mux/mux0);
        //
        valx = arma::accu(mux % Ubarx) - sigma*(mux0*std::log(mux0) + arma::accu(mux % arma::log(mux));
    }
}
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

//#include "trame_aux.hpp"

// RSC class
class RSC
{
    public:
        // build_logit objects
        int nbX;
        int nbY;
        int nbParams;
        bool outsideOption;
        
        double* dist_pars;
        
        arma::mat zeta;
        arma::mat aux_ord;
        
        arma::cube aux_Influence_lhs;
        arma::cube aux_Influence_rhs;
        
        arma::cube aux_DinvPsigma; 
        arma::cube aux_Psigma;
        
        double   (*aux_cdf_eps)(double x, double* dist_pars);
        double (*aux_quant_eps)(double x, double* dist_pars);
       
        double (*aux_pdf_eps)(double x, double* dist_pars);
        double (*aux_pot_eps)(double x, double* dist_pars);
        
        // input objects
        arma::mat mu;
        arma::mat U;
        
        // equilibrium objects
        arma::mat mu_sol;
        
        // member functions
        void build(arma::mat zeta_b, bool outsideOption_b, 
                   double   (*aux_cdf_eps_b)(double x, double* dist_pars),
                   double (*aux_quant_eps_b)(double x, double* dist_pars),
                   double   (*aux_pdf_eps_b)(double x, double* dist_pars),
                   double   (*aux_pot_eps_b)(double x, double* dist_pars));

        void build_beta(arma::mat zeta_b, double alpha, double beta);
        
        double G(arma::vec n);
        double Gx(arma::vec& mu_x, int x);
        
        double cdf (double x);
        double pdf (double x);
        double quantile (double x);
        double pot (double x);
        
    //private:
};

void RSC::build(arma::mat zeta_b, bool outsideOption_b, 
                double   (*aux_cdf_eps_b)(double x, double* dist_pars),
                double (*aux_quant_eps_b)(double x, double* dist_pars),
                double   (*aux_pdf_eps_b)(double x, double* dist_pars),
                double   (*aux_pot_eps_b)(double x, double* dist_pars))
{
    if (!outsideOption_b) {
        return;
    }
    //
    int i,j;
    //
    aux_cdf_eps   = aux_cdf_eps_b;
    aux_quant_eps = aux_quant_eps_b;
    aux_pdf_eps   = aux_pdf_eps_b;
    aux_pot_eps   = aux_pot_eps_b;
    ///
    nbX = zeta_b.n_rows;
    nbY = zeta_b.n_cols - 1;
    
    zeta = zeta_b;
    //
    aux_ord = arma::zeros(nbX,nbY+1);
    
    arma::mat D =  arma::eye(nbY+1,nbY+1) + arma::join_cols(arma::zeros(1,nbY+1), arma::join_rows(arma::eye(nbY,nbY),arma::zeros(nbY,1)));
    arma::mat D_inv = arma::inv(D);
    
    arma::mat neg_ones(nbY,1);
    neg_ones.fill(-1);
    arma::mat N_temp = arma::join_rows(arma::eye(nbY,nbY),neg_ones);
    //
    arma::uvec ordx_temp;
    arma::mat Psigmax(nbY+1,nbY+1);
    
    aux_Influence_lhs.set_size(nbY,nbY+1,nbX);
    aux_Influence_rhs.set_size(nbY+1,nbY+1,nbX);
    aux_Psigma.set_size(nbY+1,nbY+1,nbX);
    aux_DinvPsigma.set_size(nbY+1,nbY+1,nbX);
    //
    for (i=0; i<nbX; i++) {
        ordx_temp = arma::sort_index(zeta_b.row(i));
        aux_ord.row(i) = arma::conv_to< arma::rowvec >::from(ordx_temp);
        
        Psigmax.zeros();
        for (j=0; j<nbY+1; j++) {
            Psigmax(j,ordx_temp(j)) = 1;
        }
        
        aux_Influence_lhs.slice(i) = N_temp * (Psigmax.t() * D_inv);
        aux_Influence_rhs.slice(i) = D * Psigmax;
        
        aux_Psigma.slice(i) = arma::eye(nbY+1,nbY+1) * Psigmax;
        aux_DinvPsigma.slice(i) = D_inv * Psigmax;
    }
    //
    outsideOption = true;
}

// epsilon is a beta(alpha,beta) distribution
void RSC::build_beta(arma::mat zeta_b, double alpha, double beta)
{
    dist_pars = new double[2];
    dist_pars[0] = alpha;
    dist_pars[1] = beta;
    
    RSC::build(zeta_b,true,pbeta,qbeta,dbeta,iqbeta);
}

double RSC::G(arma::vec n)
{   
    int i;
    double val=0.0, val_x_temp;
    
    mu_sol.set_size(nbX,nbY);
    arma::vec mu_x;
    //
    for(i=0; i<nbX; i++){
        val_x_temp = Gx(mu_x,i);
        //
        val += n(i)*val_x_temp;
        mu_sol.row(i) = arma::trans(n(i)*mu_x);
    }
    //
    return val;
}

double RSC::Gx(arma::vec& mu_x, int x)
{
    int nbAlt = nbY + 1;
    int i,j,y,z;
    
    double val_x=0, E_eps_temp=0, E_eps_temp_next=0, cumul_temp=0; 
    double run_max=0, run_min=0, run_temp=0;
    double mu_x_tilde_y, e_y;
    
    arma::vec U_x = U.row(x).t();
    
    arma::vec mu_x_tilde = arma::zeros(nbAlt,1);
    arma::vec U_x_tilde = arma::join_cols(arma::vectorise(U_x),arma::zeros(1,1));
    //
    for (i=0; i<nbAlt; i++) {
        y = aux_ord(x,i);
        run_max = quantile(0);
        //
        j = 0;
        
        while (j < i) {
            z = aux_ord(x,j);
            
            if (zeta(x,z) != zeta(x,y)) {
                run_temp = (U_x_tilde(y) - U_x_tilde(z)) / (zeta(x,z) - zeta(x,y)); 
                run_max = std::max(run_max,run_temp);
            } else {
                run_max = INFINITY;
            }
            
            j++;
        }
        //
        run_min = quantile(1);
        //
        j = nbAlt-1;
        
        while (j > i) {
            z = aux_ord(x,j);
            
            if (zeta(x,z) != zeta(x,y)) {
                run_temp = (U_x_tilde(y) - U_x_tilde(z)) / (zeta(x,z) - zeta(x,y)); 
                run_min = std::min(run_min,run_temp);
            }
            
            j--;
        }
        //
        if (run_min > run_max) {
            mu_x_tilde_y = std::max(cdf(run_min) - cdf(run_max),0.0);
        } else {
            mu_x_tilde_y = 0;
        }
        
        mu_x_tilde(y) = mu_x_tilde_y;
        //
        if (mu_x_tilde_y > 0) {
            cumul_temp += mu_x_tilde_y;
            e_y = quantile(cumul_temp);
            
            E_eps_temp_next = e_y * cumul_temp - pot(cumul_temp);
            
            val_x += mu_x_tilde_y*U_x_tilde(y) + zeta(x,y)*(E_eps_temp_next - E_eps_temp);
        }
    }
    //
    mu_x = mu_x_tilde.rows(0,nbAlt-2);
    //
    return val_x;
}


double RSC::cdf(double x)
{
    double res;
    
    res = (*aux_cdf_eps)(x, dist_pars);
    
    return res;
}

double RSC::pdf(double x)
{
    double res;
    
    res = (*aux_pdf_eps)(x, dist_pars);
    
    return res;
}

double RSC::quantile(double x)
{
    double res;
    
    res = (*aux_quant_eps)(x, dist_pars);
    
    return res;
}

double RSC::pot(double x)
{
    double res;
    
    res = (*aux_pot_eps)(x, dist_pars);
    
    return res;
}
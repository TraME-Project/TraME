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

#include "trame_aux.hpp"

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
        
        // equilibrium objects
        arma::mat mu;
        arma::mat mux;
        arma::mat U;
        arma::mat Ux;
        
        // member functions
        void build(arma::mat sigma_b, bool outsideOption_b);
        
        
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
    arma::uvec ordx_temp;
    //
    aux_cdf_eps   = aux_cdf_eps_b;
    aux_quant_eps = aux_quant_eps_b;
    aux_pdf_eps   = aux_pdf_eps_b;
    aux_pot_eps   = aux_pot_eps_b;
    //
    nbX = zeta.n_rows;
    nbX = zeta.n_cols - 1;
    //
    aux_ord = arma::zeros(nbX,nbY+1)
    
    arma::mat D =  arma::eye(nbY+1,nbY+1) + arma::join_cols(arma::zeros(1,nbY+1), arma::join_rows(arma::eye(nbY,nbY),arma::zeros(1,nbY)));
    arma::mat D_inv = arma::inv(D);
    
    arma::mat neg_ones(nbY,1);
    neg_ones.fill(-1);
    arma::mat N_temp = arma::join_rows(arma::eye(nbY,nbY),neg_ones);
    //
    aux_Influence_lhs.set_size(nbY,nbY+1,nbX);
    aux_Influence_rhs.set_size(nbY+1,nbY+1,nbX);
    aux_Psigma.set_size(nbY+1,nbY+1,nbX);
    aux_DinvPsigma.set_size(nbY+1,nbY+1,nbX);
    //
    for (i=0; i<nbX; i++) {
        ordx_temp = arma::sort_index(zeta.row(i));
        aux_ord.row(i) = ordx_temp;
        
        Psigmax.zeros()
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

void RSC::build_beta(arma::mat zeta, double alpha, double beta)
{
    double dist_pars[2];
    dist_pars[0] = alpha;
    dist_pars[1] = beta;
    
    RSC::build(zeta,true,pbeta,qbeta,dbeta,iqbeta);
}
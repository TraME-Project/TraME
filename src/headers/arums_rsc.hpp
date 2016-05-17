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
class RSC
{
    public:
        // build_logit objects
        int nbX;
        int nbY;
        int nbParams;
        double zeta;
        bool outsideOption;
        
        double aux_Influence_lhs;
        double aux_Influence_rhs;
        
        double aux_DinvPsigma; 
        double aux_Psigma;
        
        double aux_cdf_eps;
        double aux_quant_eps;
       
        double aux_pdf_eps;
        double aux_pot_eps;
        
        // equilibrium objects
        arma::mat mu;
        arma::mat mux;
        arma::mat U;
        arma::mat Ux;
        
        // member functions
        //void build(int nbX_b, int nbY_b, int nbParams_b, double sigma_b, bool outsideOption_b);
        
        
    //private:
};
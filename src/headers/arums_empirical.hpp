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

//#include "zeroin.hpp"

// logit class
class empirical
{
    public:
        // build_logit objects
        int nbX;
        int nbY;
        
        int nbParams;
        int aux_nbDraws;
        
        bool xHomogenous;
        bool outsideOption;
        
        arma::cube atoms;
        
        // equilibrium objects
        arma::mat mux;
        
        // member functions
        void build(int nbX_b, int nbY_b, arma::cube atoms_b, bool xHomogenous_b, bool outsideOption_b);
        double Gx(arma::mat Ux, int x);
        
    //private:
};

void empirical::build(int nbX_b, int nbY_b, arma::cube atoms_b, bool xHomogenous_b, bool outsideOption_b)
{   
    nbX = nbX_b;
    nbY = nbY_b;
    
    atoms = atoms_b;
    
    nbParams = atoms_b.n_elem;
    aux_nbDraws = atoms.n_rows;
    
    xHomogenous = xHomogenous_b;
    outsideOption = outsideOption_b;
}

double empirical::Gx(arma::mat Ux, int x)
{   
    arma::mat Uxs, Utilde;
    
    if(outsideOption){
        Uxs = arma::join_cols(arma::vectorize(Ux));
    }else{
        Uxs = Ux;
    }
    
    if(xHomogenous){
        Utilde = arma::ones(aux_nbDraws,1) * Uxs.t() + atoms.slice(0);
    }else{
        Utilde = arma::ones(aux_nbDraws,1) * Uxs.t() + atoms.slice(x-1);
    }
    //
    int tt;
    arma::vec argmaxs = arma::max(Utilde,1);
    arma::uvec argmax_inds = which_max(&Utilde, 1);
    
    double thesum = 0.0;
    for(tt=0; tt < aux_nbDraws; tt++){
        thesum += argmaxs(tt,0);
    }
    double valx = thesum/(double)(aux_nbDraws);
    //
    mux.set_size(nbY,1);
    arma::uvec temp_find;
    for(tt=0; tt < nbY; tt++){
        temp_find = arma::find(argmax_inds==tt)
        mux(tt,0) = (double)(temp_find.n_elem)/(double)(aux_nbDraws);
    }
    //
    return valx;
}
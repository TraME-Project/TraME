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

#include "which_max.hpp"
#include "generic_lp.hpp"

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
        arma::mat U;
        arma::mat mu;
        arma::mat mux;
        
        // member functions
        void build(int nbX_b, int nbY_b, arma::cube atoms_b, bool xHomogenous_b, bool outsideOption_b);
        
        double G(arma::vec n);
        
        double Gstar(arma::vec n, arma::mat& U_inp);
        
        double Gx(arma::mat Ux, int x);
        double Gx(arma::mat Ux, arma::mat& mux_inp, int x);
        
        double Gstarx(arma::mat mux, arma::mat& Ux_inp, int x);
        
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

double empirical::G(arma::vec n)
{   
    int i;
    double val=0.0, valx_temp;
    
    mu.set_size(nbX,nbY);
    arma::mat mux_temp;
    //
    for(i=0; i<nbX; i++){
        valx_temp = empirical::Gx(U.row(i).t(),mux_temp,i);
        //
        val += n(i)*valx_temp;
        mu.row(i) = arma::trans(n(i)*mux_temp);
    }
    //
    return val;
}

double empirical::Gstar(arma::vec n, arma::mat& U_inp)
{   
    int i;
    double val=0.0, val_temp;
    
    U_inp.set_size(nbX,nbY);
    arma::mat Ux_temp;
    //
    for(i=0; i<nbX; i++){
        val_temp = empirical::Gstarx((mu.row(i).t())/n(i),Ux_temp,i);
        //
        val += n(i)*val_temp;
        U_inp.row(i) = arma::trans(Ux_temp);
    }
    //
    return val;
}

double empirical::Gx(arma::mat Ux, int x)
{   
    arma::mat Uxs, Utilde;
    
    if(outsideOption){
        Uxs = arma::join_cols(arma::vectorise(Ux),arma::zeros(1,1));
    }else{
        Uxs = Ux;
    }
    
    if(xHomogenous){
        Utilde = arma::ones(aux_nbDraws,1) * Uxs.t() + atoms.slice(0);
    }else{
        Utilde = arma::ones(aux_nbDraws,1) * Uxs.t() + atoms.slice(x);
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
        temp_find = arma::find(argmax_inds==tt);
        mux(tt,0) = (double)(temp_find.n_elem)/(double)(aux_nbDraws);
    }
    //
    return valx;
}

double empirical::Gx(arma::mat Ux, arma::mat& mux_inp, int x)
{   
    arma::mat Uxs, Utilde;
    
    if(outsideOption){
        Uxs = arma::join_cols(arma::vectorise(Ux),arma::zeros(1,1));
    }else{
        Uxs = Ux;
    }
    
    if(xHomogenous){
        Utilde = arma::ones(aux_nbDraws,1) * Uxs.t() + atoms.slice(0);
    }else{
        Utilde = arma::ones(aux_nbDraws,1) * Uxs.t() + atoms.slice(x);
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
    mux_inp.set_size(nbY,1);
    arma::uvec temp_find;
    for(tt=0; tt < nbY; tt++){
        temp_find = arma::find(argmax_inds==tt);
        mux_inp(tt,0) = (double)(temp_find.n_elem)/(double)(aux_nbDraws);
    }
    //
    return valx;
}

double empirical::Gstarx(arma::mat mux, arma::mat& Ux_inp, int x)
{   
    int nbOptions, jj;
    double valx=0.0;
    arma::mat Phi, Ux_temp;
    
    if(outsideOption){
        nbOptions = nbY + 1;
    }else{
        nbOptions = nbY;
    }
    
    if(xHomogenous){
        Phi = atoms.slice(0);
    }else{
        Phi = atoms.slice(x);
    }
    //
    arma::vec p = arma::ones(aux_nbDraws,1)/aux_nbDraws;
    arma::mat q;
    
    if(outsideOption){
        arma::mat temp_q(1,1); 
        temp_q(0,0) = 1 - arma::accu(mux);
        q = arma::join_cols(arma::vectorise(mux),temp_q);
    }else{
        q = arma::vectorise(mux);
    }
    //
    arma::vec obj_grbi = arma::vectorise(Phi);
    
    arma::mat A1 = arma::kron(arma::ones(1,nbOptions),arma::eye(aux_nbDraws,aux_nbDraws));
    arma::mat A2 = arma::kron(arma::eye(nbOptions,nbOptions),arma::ones(1,aux_nbDraws));
    
    arma::mat A_grbi = arma::join_cols(A1,A2);
    
    arma::vec rhs_grbi = arma::join_cols(p,q);
    
    char* sense_grbi = new char[A_grbi.n_rows];
    for(jj=0;jj<A_grbi.n_rows;jj++){
        sense_grbi[jj] = '=';
    }
    
    bool LP_optimal;
    int modelSense = 1; // maximize
    double objval;
    
    arma::mat sol_mat(obj_grbi.n_elem,2);
    arma::mat dual_mat(A_grbi.n_rows,2);
    
    try {
        LP_optimal = generic_LP(&obj_grbi, &A_grbi, modelSense, &rhs_grbi, sense_grbi, NULL, NULL, NULL, NULL, objval, sol_mat, dual_mat);
        //
        arma::mat u = dual_mat.col(0).rows(0,aux_nbDraws-1);
    
        if(outsideOption){
            Ux_temp = dual_mat.col(0).rows(aux_nbDraws,aux_nbDraws+nbY);
            Ux_inp = - Ux_temp.rows(0,nbY-1) + arma::as_scalar(Ux_temp.row(nbY));
        }else{
            Ux_temp = dual_mat.col(0).rows(aux_nbDraws,aux_nbDraws+nbY-1);
            Ux_inp = -Ux_temp - arma::accu(p % u);
        }
        //
        valx = -objval;
    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //
    delete[] sense_grbi;
    //
    return valx;
}
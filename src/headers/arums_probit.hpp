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

class probit
{
    public:
        // build objects
        int nbX;
        int nbY;
        int nbParams;
        int aux_nbOptions;
        
        bool outsideOption;
        
        double rho;
        
        arma::cube Covar;
        
        // member functions
        void build(int nbX_inp, int nbY_inp, bool outsideOption_inp);
        void simul(empirical &ret, int nbDraws, int seed);
        void unifCorrelCovMatrices();
        arma::cube unifCorrelCovMatrices(double rho);
};

// for convenience:
void probit::build(int nbX_inp, int nbY_inp, bool outsideOption_inp)
{   
    nbX = nbX_inp;
    nbY = nbY_inp;
    outsideOption = outsideOption_inp;
    //
    if (outsideOption_inp) {
        aux_nbOptions = nbY + 1;
    } else {
        aux_nbOptions = nbY;
    }
    nbParams = (nbX_inp * aux_nbOptions * (aux_nbOptions-1))/2;
}

void probit::simul(empirical &ret, int nbDraws, int seed_val)
{
    arma::arma_rng::set_seed(seed_val);
    //
    int j;
    arma::vec V;
    arma::mat Q, SqrtCovar;
    arma::cube atoms(nbDraws,aux_nbOptions,nbX);
    
    for (j=0; j<nbX; j++) {
        eig_sym(V, Q, Covar.slice(j));
        SqrtCovar = Q * arma::diagmat(1.0/arma::sqrt(V)) * Q.t();
        //
        atoms.slice(j) = arma::randn(nbDraws,aux_nbOptions) * SqrtCovar;
    }
    //
    ret.nbX = nbX;
    ret.nbY = nbY;
    ret.nbParams = atoms.n_elem;
    ret.atoms = atoms;
    ret.aux_nbDraws = nbDraws;
    ret.xHomogenous = false;
    ret.outsideOption = outsideOption;
    if (outsideOption) {
        ret.nbOptions = nbY + 1;
    } else {
        ret.nbOptions = nbY;
    }
    //
    arma::arma_rng::set_seed_random(); // need to reset the seed
}

void probit::unifCorrelCovMatrices()
{
    int i;
    //
    arma::mat Sig = rho * arma::ones(aux_nbOptions,aux_nbOptions) + (1-rho) * arma::eye(aux_nbOptions,aux_nbOptions);
    //
    if (outsideOption) {
        Sig.col(aux_nbOptions-1).fill(0);
        Sig.row(aux_nbOptions-1).fill(0);
        Sig(aux_nbOptions-1,aux_nbOptions-1) = 1;
    }
    //
    Covar.set_size(aux_nbOptions,aux_nbOptions,nbX); // note: this is different to the R code
    for (i=0; i<nbX; i++) {
        Covar.slice(i) = Sig;
    }
}

arma::cube probit::unifCorrelCovMatrices(double rho_inp)
{
    int i;
    //
    arma::mat Sig = rho_inp * arma::ones(aux_nbOptions,aux_nbOptions) + (1-rho_inp) * arma::eye(aux_nbOptions,aux_nbOptions);
    //
    if (outsideOption) {
        Sig.col(aux_nbOptions-1).fill(0);
        Sig.row(aux_nbOptions-1).fill(0);
        Sig(aux_nbOptions-1,aux_nbOptions-1) = 1;
    }
    //
    arma::cube Covar_ret(aux_nbOptions,aux_nbOptions,nbX); // note: this is different to the R code
    for (i=0; i<nbX; i++) {
        Covar_ret.slice(i) = Sig;
    }
    //
    return Covar_ret;
}
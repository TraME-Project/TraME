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

class none
{
    public:
        // build_logit objects
        int nbX;
        int nbY;
        int nbParams;
        
        // equilibrium objects
        arma::mat mu;
        arma::mat mux;
        arma::mat U;
        arma::mat Ux;
        
        // member functions
        void build(int nbX_b, int nbY_b);
        double Gx(arma::mat Ux);
        arma::vec none::dtheta_NablaGstar(arma::mat mu, arma::vec n);
        
        void simul(empirical &ret, int nbDraws, int seed);
};

// for convenience:
void none::build(int nbX_b, int nbY_b)
{   
    nbX = nbX_b;
    nbY = nbY_b;
    nbParams = 0;
}

double none::Gx(arma::mat Ux_inp)
{
    arma::uvec y = which_max(&Ux, (int) 0);
    //
    mux.set_size(nbY,1);
    if(y < nbY){
        mux(y) = 1;
    }
    //
    double valx = std::max(arma::as_scalar(arma::max(arma::vectorise(Ux))), (double) 0.0);
    //
    return valx;
}

arma::vec none::dtheta_NablaGstar(arma::mat mu, arma::vec n)
{
    arma::vec ret = arma::zeros(nbX*nbY,1);
    return ret;
}

void none::simul(empirical &ret, int nbDraws, int seed_val)
{
    arma::arma_rng::set_seed(seed_val);
    //
    arma::cube atoms(nbDraws,nbY+1,nbX);
    atoms.fill(0);
    //
    ret.nbX = nbX;
    ret.nbY = nbY;
    ret.nbParams = atoms.n_elem;
    ret.atoms = atoms;
    ret.aux_nbDraws = nbDraws;
    ret.xHomogenous = false;
    //
    arma::arma_rng::set_seed_random(); // need to reset the seed
}

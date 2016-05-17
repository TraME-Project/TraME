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
        arma::vec dtheta_NablaGstar();
        double Gbarx(arma::mat Ubarx, arma::mat mubarx, arma::mat& Ux_inp, arma::mat& mux_inp);
        
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
    arma::uvec temp_vec = which_max(&Ux_inp, (int) 0);
    int y = temp_vec(0);
    //
    mux.zeros(nbY,1);
    if(y < nbY){
        mux(y) = 1;
    }
    //
    double valx = std::max(arma::as_scalar(arma::max(arma::vectorise(Ux_inp))), (double) 0.0);
    //
    return valx;
}

arma::vec none::dtheta_NablaGstar()
{
    arma::vec ret = arma::zeros(nbX*nbY,1);
    return ret;
}

double none::Gbarx(arma::mat Ubarx, arma::mat mubarx, arma::mat& Ux_inp, arma::mat& mux_inp)
{
    int count_int=0;
    int nbY0 = Ubarx.n_elem;
    //
    //arma::mat srt = arma::sort(Ubarx,"descend");
    arma::uvec srt_ind = arma::sort_index(Ubarx,"descend");
    //
    mux_inp.set_size(nbY0,1);
    double cumul = arma::as_scalar(mubarx(srt_ind(count_int)));
    //
    while((count_int < nbY0-1) & (cumul < 1.0) & (Ubarx(srt_ind(count_int)) > 0)){
        mux_inp(srt_ind(count_int)) = mubarx(srt_ind(count_int));
        count_int++;
        cumul += mubarx(srt_ind(count_int)); // Keith: is this in the correct place?
    }
    //
    if(Ubarx(srt_ind(count_int)) > 0){
        mux_inp(srt_ind(count_int)) = mubarx(srt_ind(count_int)) + 1 - cumul;
    }
    //
    Ux = arma::zeros(nbY0,1);
    //
    double valx = arma::accu(mux % Ubarx);
    //
    return valx;
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

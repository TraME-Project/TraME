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
        
        // input objects
        arma::mat mu;
        arma::mat U;
        
        // equilibrium objects
        arma::mat mu_sol;
        arma::mat U_sol;
        
        // member functions
        void build(int nbX_b, int nbY_b);
        double G(arma::vec n);
        double Gx(arma::mat Ux, arma::mat& mu_x_inp);

        double Gstar(arma::mat& U_inp, arma::vec n);
        
        double Gbar(arma::mat Ubarx, arma::mat mubarx, arma::vec n, arma::mat& Ux_inp, arma::mat& mu_x_inp);
        double Gbarx(arma::mat Ubarx, arma::mat mubarx, arma::mat& Ux_inp, arma::mat& mu_x_inp);
        
        arma::vec dtheta_NablaGstar();
        
        void simul(empirical &ret, int nbDraws, int seed);
};

// for convenience:
void none::build(int nbX_b, int nbY_b)
{   
    nbX = nbX_b;
    nbY = nbY_b;
    nbParams = 0;
}

double none::G(arma::vec n)
{   
    int i;
    double val=0.0, val_x;
    
    mu_sol.set_size(nbX,nbY);
    arma::mat mu_x_temp;
    //
    for(i=0; i<nbX; i++){
        val_x = none::Gx(U.row(i).t(),mu_x_temp);
        //
        val += n(i)*val_x;
        mu_sol.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double none::Gx(arma::mat Ux, arma::mat& mu_x_inp)
{
    arma::uvec temp_vec = which_max(&Ux, (int) 0);
    int y = temp_vec(0);
    //
    mu_x_inp.set_size(nbY,1);
    mu_x_inp.zeros();
    
    if (y < nbY) {
        mu_x_inp(y) = 1;
    }
    //
    double val_x = std::max(arma::as_scalar(arma::max(arma::vectorise(Ux))), (double) 0.0);
    //
    return val_x;
}

double none::Gstar(arma::mat& U_inp, arma::vec n)
{   
    printf("Gstar not yet defined for no arums case.\n");

    return 0.0;
}

double none::Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_inp, arma::mat& mu_inp)
{   
    int i;
    double val=0.0, val_temp;
    
    U_inp.set_size(nbX,nbY);
    mu_inp.set_size(nbX,nbY);
    arma::mat Ux_temp, mu_x_temp;
    //
    for(i=0; i<nbX; i++){
        val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),Ux_temp,mu_x_temp);
        //
        val += n(i)*val_temp;
        U_inp.row(i) = arma::trans(Ux_temp);
        mu_inp.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double none::Gbarx(arma::mat Ubarx, arma::mat mubarx, arma::mat& Ux_inp, arma::mat& mu_x_inp)
{
    int count_int=0;
    int nbY0 = Ubarx.n_elem;
    //
    //arma::mat srt = arma::sort(Ubarx,"descend");
    arma::uvec srt_ind = arma::sort_index(Ubarx,"descend");
    //
    mu_x_inp.set_size(nbY0,1);
    double cumul = arma::as_scalar(mubarx(srt_ind(count_int)));
    //
    while((count_int < nbY0-1) & (cumul < 1.0) & (Ubarx(srt_ind(count_int)) > 0)){
        mu_x_inp(srt_ind(count_int)) = mubarx(srt_ind(count_int));
        count_int++;
        cumul += mubarx(srt_ind(count_int)); // Keith: is this in the correct place?
    }
    //
    if(Ubarx(srt_ind(count_int)) > 0){
        mu_x_inp(srt_ind(count_int)) = mubarx(srt_ind(count_int)) + 1 - cumul;
    }
    //
    Ux_inp = arma::zeros(nbY0,1);
    //
    double valx = arma::accu(mu_x_inp % Ubarx);
    //
    return valx;
}

arma::vec none::dtheta_NablaGstar()
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

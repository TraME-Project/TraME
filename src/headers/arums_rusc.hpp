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

// RUSC class
class RUSC
{
    public:
        // build_logit objects
        int nbX;
        int nbY;
        int nbParams;
        bool outsideOption;
        
        arma::mat zeta;
        arma::mat aux_ord;
        
        arma::cube aux_A;
        arma::mat  aux_b;
        arma::vec  aux_c; 
        
        // input objects
        arma::mat mu;
        arma::mat U;
        
        // equilibrium objects
        arma::mat mu_sol;
        arma::mat U_sol;
        
        // member functions
        void build(arma::mat zeta_b, bool outsideOption_b);
        
        double G(arma::vec n);
        double Gx(arma::vec& mu_x, int x);
        
        double Gstar(arma::vec n);
        double Gstarx(arma::vec& U_x, double n_x, int x);
        
        double Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_inp, arma::mat& mu_inp);
        double Gbarx(arma::vec Ubarx, arma::vec mubarx, arma::mat& Ux_inp, arma::mat& mux_inp, int x);
        
        void simul(empirical &ret, int nbDraws, int seed);
        
    //private:
};

void RUSC::build(arma::mat zeta_b, bool outsideOption_b)
{
    if (!outsideOption_b) {
        printf("outsideOption=F not implemented yet on RUSC arums\n");
        return;
    }
    //
    int i,j;
    //
    nbX = zeta_b.n_rows;
    nbY = zeta_b.n_cols - 1;
    nbParams = zeta_b.n_elem;
    
    zeta = zeta_b;
    //
    aux_ord = arma::zeros(nbX,nbY+1);

    aux_A.set_size(nbY,nbY+1,nbX);
    aux_A.zeros();

    aux_b = arma::zeros(nbX,nbY);
    aux_c = arma::zeros(nbY,1);
    //
    double z_0;
    arma::mat max_z, max_z0_mat;
    arma::vec z_x_temp, z_x, max_z0;
    arma::uvec ordx_temp;
    
    for (i=0; i<nbX; i++) {
        z_x_temp = zeta.row(i).t();

        z_x = z_x_temp.shed_row(nbY);
        z_0 = z_x(nbY);

        max_z  = arma::max(z_x * arma::ones(1,nbY), arma::ones(nbY,1) * z_x.t());
        max_z0 = arma::max(z,z_0);
        max_z0_mat = max_z0 * arma::ones(1,nbY);

        A_x = max_z0_mat + max_z0_mat.t() - max_z - z_0; 
        //
        aux_A.slice(i) = A_x;
        aux_b.row(i) = z_0 - max_z0.t();
        aux_c(i) = -z_0/2;
        
        ordx_temp = arma::sort_index(zeta_b.row(i));
        aux_ord.row(i) = arma::conv_to< arma::rowvec >::from(ordx_temp);
    }
    //
    outsideOption = true;
}

double RUSC::G(arma::vec n)
{   
    int i;
    double val=0.0, val_x_temp;
    
    mu_sol.set_size(nbX,nbY);
    arma::vec mu_x;
    //
    for (i=0; i<nbX; i++) {
        val_x_temp = Gx(mu_x,i);
        //
        val += n(i)*val_x_temp;
        mu_sol.row(i) = arma::trans(n(i)*mu_x);
    }
    //
    return val;
}

double RUSC::Gx(arma::vec& mu_x, int x)
{
    int nbAlt = nbY + 1;
    int i,j,y,z;
    
    double val_x; 
    double run_max=0, run_min=0, run_temp=0;
    
    arma::vec U_x = U.row(x).t();
    
    arma::vec mu_x_tilde = arma::zeros(nbAlt,1);
    arma::vec U_x_tilde = arma::join_cols(arma::vectorise(U_x),arma::zeros(1,1));
    //
    for (i=0; i<nbAlt; i++) {
        y = aux_ord(x,i);
        run_max = 0.0;
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
        run_min = 1;
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
        mu_x_tilde(y) = std::max(run_min - run_max,0.0);
    }
    //
    mu_x = mu_x_tilde.rows(0,nbAlt-2);
    //
    val_x = arma::accu(mu_x % (U_x - aux_b.row(x))) - mu_x.t() * aux_A.slice(x) * mu_x/2 - aux_c(x);
    //
    return val_x;
}

double RUSC::Gstar(arma::vec n)
{   
    int i;
    double val=0.0, val_x_temp;
    
    U_sol.set_size(nbX,nbY);
    arma::vec U_x_temp;
    //
    for (i=0; i<nbX; i++) {
        //val_x_temp = Gstarx((mu.row(i).t())/n(i),U_x_temp,i);
        val_x_temp = Gstarx(U_x_temp,n(i),i);
        //
        val += n(i)*val_x_temp;
        U_sol.row(i) = arma::trans(U_x_temp);
    }
    //
    return val;
}

double RUSC::Gstarx(arma::vec& U_x, double n_x, int x)
{
    double val_x = 0;
    arma::vec mu_x = (mu_sol.row(x).t())/n_x; // we divide by n(x)
    
    arma::vec A_mu = arma::trans(aux_A.slice(x) * mu_x); 

    U_x = A_mu + aux_b.row(x).t();
    val_x = arma::accu(mu_x % A_mu)/2 + arma::accu(mu_x % aux_b.row(x).t()) + aux_c(x);
    
    return val_x;
}

double RUSC::Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_inp, arma::mat& mu_inp)
{   
    int i;
    double val=0.0, val_temp;
    
    if (!TRAME_PRESOLVED_GBAR){
        presolve_LP_Gbar();
    }
    
    U_inp.set_size(nbX,nbY);
    mu_inp.set_size(nbX,nbY);
    arma::mat Ux_temp, mux_temp;
    //
    for(i=0; i<nbX; i++){
        val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),Ux_temp,mux_temp,i);
        //
        val += n(i)*val_temp;
        U_inp.row(i) = arma::trans(Ux_temp);
        mu_inp.row(i) = arma::trans(n(i)*mux_temp);
    }
    //
    return val;
}

double RUSC::Gbarx(arma::mat U_bar_x, arma::mat mu_bar_x, arma::mat& U_x_inp, arma::mat& mu_x_inp)
{
    int nbAlt = nbY + 1;

    arma::vec obj_grbi = arma::join_cols(aux_b.row(x).t(),arma::zeros(1,1));
    arma::mat A_grbi = arma::ones(1,nbAlt);
    arma::vec rhs_grbi = arma::ones(1,1);

    arma::mat Q_grbi = arma::zeros(nbAlt,nbAlt);
    Q_grbi.submat(0,0,nbAlt-2,nbAlt-2) = aux_A.slice(x) / 2;

    arma::vec lb_grbi = arma::zeros(1,nbAlt);
    arma::vec ub_grbi = arma::join_cols(mu_bar_x,arma::ones(1,1));
    
    char* sense_grbi = new char[1];
    sense_grbi[0] = '=';
    
    bool LP_optimal;
    int modelSense = 0; // minimize
    double objval;
    
    arma::mat sol_mat(obj_grbi.n_elem,2);
    arma::mat dual_mat(1,2);
    
    try {
        LP_optimal = generic_LP((int) A_grbi.n_rows, (int) A_grbi.n_cols, obj_grbi.memptr(), A_grbi.memptr(), modelSense, rhs_grbi.memptr(), sense_grbi, Q_grbi.memptr(), lb_grbi.memptr(), ub_grbi.memptr(), NULL, objval, sol_mat, dual_mat);
        //
        if (LP_optimal) {
            mu_x_inp = sol_mat(arma::span(0,nbAlt-2),0);

            A_mu = arma::trans(aux_A.slice(x) * mu_x_inp);

            U_x_inp = A_mu + aux_b.row(x).t();
            //
            valx = -objval - aux_c(x);
        } else {
            std::cout << "Non-optimal value found during optimization" << std::endl;
        }
#if !defined(TRAME_USE_GUROBI_C)
    } catch(GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << std::endl;
        std::cout << e.getMessage() << std::endl;
#endif
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    
    return val_x;
}

void RSC::simul(empirical &ret, int nbDraws, int seed_val)
{
    int i;
    arma::arma_rng::set_seed(seed_val);
    //
    arma::cube atoms(nbDraws,nbY+1,nbX);
    
    for (i=0; i<nbX; i++) {
        atoms.slice(i) = arma::randu(nbDraws,1) * zeta.row(i);
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

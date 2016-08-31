/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
  ##      Simon Weber
  ##
  ##   This file is part of TraME.
  ##
  ##   TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * RUSC class
 *
 * Keith O'Hara
 * 08/08/2016
 */

#include "trame.hpp"

trame::rusc::rusc(arma::mat zeta_inp, bool outsideOption_inp)
{
    this->build(zeta_inp, outsideOption_inp);
}

void trame::rusc::build(arma::mat zeta_inp, bool outsideOption_inp)
{
    if (!outsideOption_inp) {
        printf("outsideOption=F not implemented yet on RUSC arums\n");
        return;
    }
    //
    int i;
    //
    nbX = zeta_inp.n_rows;
    nbY = zeta_inp.n_cols - 1;
    nbParams = zeta_inp.n_elem;
    
    zeta = zeta_inp;
    //
    aux_ord = arma::zeros(nbX,nbY+1);

    aux_A.set_size(nbY,nbY+1,nbX);
    aux_A.zeros();

    aux_b = arma::zeros(nbX,nbY);
    aux_c = arma::zeros(nbY,1);
    //
    double z_0;
    arma::mat A_x, max_z, max_z0_mat;
    arma::vec z_x, max_z0;
    arma::uvec ordx_temp;
    
    for (i=0; i<nbX; i++) {
        z_x = zeta.row(i).t();

        z_x.shed_rows(nbY,zeta.n_cols-1);
        z_0 = z_x(nbY);

        max_z  = arma::max(z_x * arma::ones(1,nbY), arma::ones(nbY,1) * z_x.t());
        max_z0 = arma::max(z_x,z_0);
        max_z0_mat = max_z0 * arma::ones(1,nbY);

        A_x = max_z0_mat + max_z0_mat.t() - max_z - z_0; 
        //
        aux_A.slice(i) = A_x;
        aux_b.row(i) = z_0 - max_z0.t();
        aux_c(i) = -z_0/2;
        
        ordx_temp = arma::sort_index(zeta_inp.row(i));
        aux_ord.row(i) = arma::conv_to< arma::rowvec >::from(ordx_temp);
    }
    //
    outsideOption = true;
}

double trame::rusc::G(arma::vec n)
{   
    double val = this->G(n,U,mu_sol);
    //
    return val;
}

double trame::rusc::G(arma::vec n, const arma::mat& U_inp, arma::mat& mu_out)
{   
    int i;
    double val=0.0, val_x_temp;
    
    mu_out.set_size(nbX,nbY);
    arma::vec mu_x;
    //
    for (i=0; i<nbX; i++) {
        val_x_temp = Gx(U_inp.row(i).t(), mu_x, i);
        //
        val += n(i)*val_x_temp;
        mu_out.row(i) = arma::trans(n(i)*mu_x);
    }
    //
    return val;
}

double trame::rusc::Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, int x)
{
    int nbAlt = nbY + 1;
    int i,j,y,z;
    
    double val_x; 
    double run_max=0, run_min=0, run_temp=0;
    
    arma::vec mu_x_tilde = arma::zeros(nbAlt,1);
    arma::vec U_x_tilde = arma::join_cols(arma::vectorise(U_x_inp),arma::zeros(1,1));
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
    mu_x_out = mu_x_tilde.rows(0,nbAlt-2);
    //
    val_x = arma::accu(mu_x_out % (U_x_inp - aux_b.row(x))) - arma::as_scalar(mu_x_out.t() * aux_A.slice(x) * mu_x_out/2) - aux_c(x);
    //
    return val_x;
}

double trame::rusc::Gstar(arma::vec n)
{
    double val = this->Gstar(n,mu_sol,U_sol);
    //
    return val;
}

double trame::rusc::Gstar(arma::vec n, const arma::mat& mu_inp, arma::mat& U_out)
{   
    int i;
    double val=0.0, val_x_temp;
    
    U_out.set_size(nbX,nbY);
    arma::vec U_x_temp;
    //
    for (i=0; i<nbX; i++) {
        val_x_temp = Gstarx((mu_inp.row(i).t())/n(i),U_x_temp,i);
        //
        val += n(i)*val_x_temp;
        U_out.row(i) = arma::trans(U_x_temp);
    }
    //
    return val;
}

double trame::rusc::Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x)
{
    double val_x = 0;
    
    arma::vec A_mu = arma::trans(aux_A.slice(x) * mu_x_inp); 

    U_x_out = A_mu + aux_b.row(x).t();
    val_x = arma::accu(mu_x_inp % A_mu)/2 + arma::accu(mu_x_inp % aux_b.row(x).t()) + aux_c(x);
    
    return val_x;
}

double trame::rusc::Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_out, arma::mat& mu_out)
{   
    int i;
    double val=0.0, val_temp;
    
    U_out.set_size(nbX,nbY);
    mu_out.set_size(nbX,nbY);
    arma::mat U_x_temp, mux_temp;
    //
    for (i=0; i<nbX; i++) {
        val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),U_x_temp,mux_temp,i);
        //
        val += n(i)*val_temp;
        U_out.row(i) = arma::trans(U_x_temp);
        mu_out.row(i) = arma::trans(n(i)*mux_temp);
    }
    //
    return val;
}

double trame::rusc::Gbarx(arma::mat Ubar_x, arma::mat mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, int x)
{
    int nbAlt = nbY + 1;
    double val_x = 0.0;

    arma::vec obj_grbi = arma::join_cols(aux_b.row(x).t(),arma::zeros(1,1));
    arma::mat A_grbi = arma::ones(1,nbAlt);
    arma::vec rhs_grbi = arma::ones(1,1);

    arma::mat Q_grbi = arma::zeros(nbAlt,nbAlt);
    Q_grbi.submat(0,0,nbAlt-2,nbAlt-2) = aux_A.slice(x) / 2;

    arma::vec lb_grbi = arma::zeros(1,nbAlt);
    arma::vec ub_grbi = arma::join_cols(mubar_x,arma::ones(1,1));
    
    char* sense_grbi = new char[1];
    sense_grbi[0] = '=';
    
    bool LP_optimal;
    int modelSense = 0; // minimize
    double objval;
    
    arma::mat sol_mat(obj_grbi.n_elem,2);
    arma::mat dual_mat(1,2);
    
    try {
        //LP_optimal = generic_LP((int) A_grbi.n_rows, (int) A_grbi.n_cols, obj_grbi.memptr(), A_grbi.memptr(), modelSense, rhs_grbi.memptr(), sense_grbi, Q_grbi.memptr(), lb_grbi.memptr(), ub_grbi.memptr(), NULL, objval, sol_mat, dual_mat);
        LP_optimal = generic_LP((int) A_grbi.n_rows, (int) A_grbi.n_cols, obj_grbi.memptr(), A_grbi.memptr(), modelSense, rhs_grbi.memptr(), sense_grbi, Q_grbi.memptr(), lb_grbi.memptr(), ub_grbi.memptr(), NULL, objval, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));
        //
        if (LP_optimal) {
            mu_x_out = sol_mat(arma::span(0,nbAlt-2),0);

            arma::mat A_mu = arma::trans(aux_A.slice(x) * mu_x_out);

            U_x_out = A_mu + aux_b.row(x).t();
            //
            val_x = -objval - aux_c(x);
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

trame::empirical trame::rusc::simul()
{
    empirical emp_obj;
    
    this->simul(emp_obj,NULL,NULL);
    //
    return emp_obj;
}

trame::empirical trame::rusc::simul(int* nbDraws, int* seed)
{
    empirical emp_obj;
    
    this->simul(emp_obj,nbDraws,seed);
    //
    return emp_obj;
}

void trame::rusc::simul(empirical& obj_out)
{
    this->simul(obj_out,NULL,NULL);
}

void trame::rusc::simul(empirical& obj_out, int* nbDraws, int* seed_val)
{
    int n_draws = 0;
    if (nbDraws) {
        n_draws = *nbDraws;
    } else {
#ifdef TRAME_DEFAULT_SIM_DRAWS
        n_draws = TRAME_DEFAULT_SIM_DRAWS;
#else
        n_draws = 1000;
#endif
    }
    //
    if (seed_val) {
        arma::arma_rng::set_seed(*seed_val);
    }
    //
    arma::cube atoms(n_draws,nbY+1,nbX);
    
    for (int i=0; i<nbX; i++) {
        atoms.slice(i) = arma::randu(n_draws,1) * zeta.row(i);
    }
    //
    obj_out.nbX = nbX;
    obj_out.nbY = nbY;
    obj_out.nbParams = atoms.n_elem;
    obj_out.atoms = atoms;
    obj_out.aux_nbDraws = n_draws;
    obj_out.xHomogenous = false;
    obj_out.outsideOption = outsideOption;
    if (outsideOption) {
        obj_out.nbOptions = nbY + 1;
    } else {
        obj_out.nbOptions = nbY;
    }
    //
    arma::arma_rng::set_seed_random(); // need to reset the seed
}

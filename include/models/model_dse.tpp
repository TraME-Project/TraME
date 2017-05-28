/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
  ##      Alfred Galichon
  ##      Keith O'Hara
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
 * general model class
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * This version:
 * 05/27/2017
 */

// Note: 'theta' refers to model parameters;'params' refers to the structural parameters

//
// first method to build

template<typename Tg, typename Th, typename Tt>
inline 
void 
model<dse<Tg,Th,Tt>>::build(const arma::cube& phi_xyk_inp)
{
    this->build_int(phi_xyk_inp,NULL,NULL);
}

template<typename Tg, typename Th, typename Tt>
inline 
void 
model<dse<Tg,Th,Tt>>::build(const arma::cube& phi_xyk_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    this->build_int(phi_xyk_inp,&n_inp,&m_inp);
}

template<typename Tg, typename Th, typename Tt>
inline 
void 
model<dse<Tg,Th,Tt>>::build_int(const arma::cube& phi_xyk_inp, const arma::vec* n_inp, const arma::vec* m_inp)
{
    need_norm = false;

    nbX = phi_xyk_inp.n_rows;
    nbY = phi_xyk_inp.n_cols;
    dim_theta = phi_xyk_inp.n_slices;
    //
    n = (n_inp) ? *n_inp : arma::ones(nbX,1);
    m = (m_inp) ? *m_inp : arma::ones(nbY,1);

    // phi_xyk = phi_xyk_inp;
    model_data = cube_to_mat(phi_xyk_inp); // Phi_xy = matrix(phi_xyk)
}

//
// second method to build

template<typename Tg, typename Th, typename Tt>
inline 
void 
model<dse<Tg,Th,Tt>>::build(const arma::mat& X_inp, const arma::mat& Y_inp)
{
    this->build_int(X_inp,Y_inp,NULL,NULL);
}

template<typename Tg, typename Th, typename Tt>
inline
void
model<dse<Tg,Th,Tt>>::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    this->build_int(X_inp,Y_inp,&n_inp,&m_inp);
}

template<typename Tg, typename Th, typename Tt>
inline
void
model<dse<Tg,Th,Tt>>::build_int(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec* n_inp, const arma::vec* m_inp)
{
    need_norm = false;

    nbX = X_inp.n_rows;
    nbY = Y_inp.n_rows;

    dX = X_inp.n_cols;
    dY = Y_inp.n_cols;

    dim_theta = dX*dY;
    //
    n = (n_inp) ? *n_inp : arma::ones(nbX,1);
    m = (m_inp) ? *m_inp : arma::ones(nbY,1);
    //
    model_data = model_build_int(market_obj,X_inp,Y_inp);
    arma::mat phi_xy_temp = arma::kron(Y_inp,X_inp);
    arma::cube phi_xyk_temp(phi_xy_temp.memptr(),nbX,nbY,dim_theta,false); // share memory

    // phi_xyk = phi_xyk_temp;
    model_data = cube_to_mat(phi_xyk_temp);
}

//
// build markets (TU case only right now)

template<typename Tg, typename Th, typename Tt>
void 
model<dse<Tg,Th,Tt>>::model_to_market(const arma::mat& theta)
{
    market_obj.build(n,m,Phi_xy_theta(theta),NULL,need_norm);
}

// template<typename Tm>
// void 
// model<dse<Tg,Th,Tt>>::build_market_TU(const arma::mat& theta)
// {
//     market_obj.build(n,m,Phi_xy_theta(theta),NULL,need_norm);
// }

// general arums input
// template<typename Tm>
// template<typename Ta, typename Tb>
// void 
// model<dse<Tg,Th,Tt>>::build_market_TU(const arma::mat& theta, const Ta& arums_G_inp, const Tb& arums_H_inp)
// {
//     market_obj.build(n,m,Phi_xy_theta(theta),arums_G_inp,arums_H_inp,need_norm);
// }

// empirical version
// template<typename Tm> 
// template<typename Ta, typename Tb>
// void 
// model<dse<Tg,Th,Tt>>::build_market_TU(const arma::mat& theta, const Ta& arums_G_inp, const Tb& arums_H_inp, int nbDraws, int seed)
// {
//     market_obj.build(n,m,Phi_xy_theta(theta),arums_G_inp,arums_H_inp,nbDraws,seed,need_norm);
// }

//
// gradients

template<typename Tg, typename Th, typename Tt>
void 
model<dse<Tg,Th,Tt>>::dtheta(const arma::mat* delta_theta_inp, arma::mat& dtheta_Psi_out)
{
    this->dtheta(delta_theta_inp,dtheta_Psi_out,NULL,NULL);
}

template<typename Tg, typename Th, typename Tt>
void 
model<dse<Tg,Th,Tt>>::dtheta(const arma::mat* delta_theta_inp, arma::mat& dtheta_Psi_out, arma::mat* dtheta_G_out, arma::mat* dtheta_H_out)
{
    arma::mat delta_theta_mat = (delta_theta_inp) ? *delta_theta_inp : arma::eye(dim_theta,dim_theta);

    dtheta_Psi_out = model_data * delta_theta_mat;
    //
    if (dtheta_G_out) {
        *dtheta_G_out = arma::zeros(0,delta_theta_mat.n_cols);
    }
    if (dtheta_H_out) {
        *dtheta_H_out = arma::zeros(0,delta_theta_mat.n_cols);
    }
}

template<typename Tg, typename Th, typename Tt>
void 
model<dse<Tg,Th,Tt>>::dtheta_mu(const arma::mat& theta, const arma::mat* delta_theta, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& dmu_out)
{
    model_to_market(theta); // need to replace this later with general 'parametric_market'
    //
    arma::mat dtheta_Psi;
    dtheta(delta_theta,dtheta_Psi,NULL,NULL);
    //
    model_dmu(market_obj,dtheta_Psi,mu_out,mu_x0_out,mu_0y_out,dmu_out);
}

//
// MLE

template<typename Tg, typename Th, typename Tt>
bool 
model<dse<Tg,Th,Tt>>::mle(const arma::mat& mu_hat, arma::mat& theta_hat, const arma::mat* theta_0_inp)
{
    bool success = false;
    //
    double err_tol = 1E-06;
    int max_iter = 5000;

    arma::vec theta_0;
    (theta_0_inp) ? theta_0 = *theta_0_inp : initial_theta(theta_0);

    model_to_market(theta_0);

    bool by_individual = true;
    double scale = std::max(arma::accu(n),arma::accu(m));

    arma::vec mu_hat_x0 = n - arma::sum(mu_hat,1);
    arma::vec mu_hat_0y = m - arma::trans(arma::sum(mu_hat,0));
    //
    // add optimization data
    trame_model_mle_opt_data<dse<Tg,Th,Tt>> opt_data;
    
    opt_data.model_obj = *this;
    opt_data.by_individual = by_individual;
    opt_data.scale = scale;
    
    opt_data.mu_hat = mu_hat;
    opt_data.mu_hat_x0 = mu_hat_x0;
    opt_data.mu_hat_0y = mu_hat_0y;
    //
    double obj_val = 0;

    success = model_mle_optim(theta_0,log_likelihood,&opt_data,&obj_val,&err_tol,&max_iter);
    //
    theta_hat = theta_0;

    return success;
}

//
// MME

template<typename Tg, typename Th, typename Tt>
bool 
model<dse<Tg,Th,Tt>>::mme(const arma::mat& mu_hat, arma::mat& theta_hat)
{
    bool success = false;
    //
    // double xtol_rel = (xtol_rel_inp) ? *xtol_rel_inp : 1E-04;
    // int max_eval = (max_eval_inp) ? *max_eval_inp : 1E05;
    double err_tol = 1E-06;
    int max_iter = 5000;

    arma::vec theta_0;
    initial_theta(theta_0);

    arma::mat dtheta_Psi;
    dtheta(NULL,dtheta_Psi);

    model_to_market(theta_0);

    arma::mat kron_term = dtheta_Psi;
    arma::mat C_hat = arma::vectorise(arma::trans(arma::vectorise(mu_hat)) * kron_term);
    //
    // add optimization data
    trame_model_mme_opt_data<dse<Tg,Th,Tt>> opt_data;
    
    opt_data.market = market_obj;
    opt_data.dim_theta = dim_theta;
    opt_data.C_hat = C_hat;
    opt_data.kron_term = kron_term;
    //
    arma::vec sol_vec = arma::join_cols(arma::vectorise(kron_term * theta_0)/2.0,theta_0);

    double obj_val = 0;

    success = model_mme_optim(sol_vec,model_mme_opt_objfn,&opt_data,&obj_val,&err_tol,&max_iter);
    //
    arma::mat U = arma::reshape(sol_vec.rows(0,nbX*nbY-1),nbX,nbY);
    theta_hat = sol_vec.rows(nbX*nbY,nbX*nbY+dim_theta-1);

    //double val_ret = obj_val;
    //
    return success;
}

//
// solve wrappers

template<typename Tg, typename Th, typename Tt>
bool 
model<dse<Tg,Th,Tt>>::solve(arma::mat& mu_sol)
{
    bool res = market_obj.solve(mu_sol,NULL);
    //
    return res;
}

template<typename Tg, typename Th, typename Tt>
bool 
model<dse<Tg,Th,Tt>>::solve(arma::mat& mu_sol, const char* solver)
{
    bool res = market_obj.solve(mu_sol,solver);
    //
    return res;
}

template<typename Tg, typename Th, typename Tt>
bool 
model<dse<Tg,Th,Tt>>::solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver)
{
    bool res = market_obj.solve(mu_sol,U,V,solver);
    //
    return res;
}

//
// internal

template<typename Tg, typename Th, typename Tt>
inline
arma::mat 
model<dse<Tg,Th,Tt>>::Phi_xy_theta(const arma::mat& theta)
{
    arma::mat ret = arma::reshape(model_data * theta,nbX,nbY);
    return ret;
}

template<typename Tg, typename Th, typename Tt>
inline
void 
model<dse<Tg,Th,Tt>>::initial_theta(arma::mat& params)
{
    params.zeros(dim_theta,1);
}

//
// optimization-related functions

template<typename Tg, typename Th, typename Tt>
bool 
model<dse<Tg,Th,Tt>>::model_mle_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, double* value_out, double* err_tol_inp, int* max_iter_inp)
{
    bool success = generic_optim(init_out_vals,opt_objfn,opt_data,value_out,err_tol_inp,max_iter_inp);
    //
    return success;
}

template<typename Tg, typename Th, typename Tt>
double 
model<dse<Tg,Th,Tt>>::log_likelihood(const arma::vec& vals_inp, arma::vec* grad_vec, void* opt_data)
{
    trame_model_mle_opt_data<dse<Tg,Th,Tt>> *d = reinterpret_cast<trame_model_mle_opt_data<dse<Tg,Th,Tt>>*>(opt_data);

    bool by_individual = d->by_individual;

    int nbX = d->model_obj.nbX;
    int nbY = d->model_obj.nbY;

    double scale = d->scale;

    arma::mat mu_hat = d->mu_hat;
    arma::vec mu_hat_x0 = d->mu_hat_x0;
    arma::vec mu_hat_0y = d->mu_hat_0y;
    //
    arma::vec mu_x0, mu_0y;
    arma::mat mu, dmu;

    d->model_obj.dtheta_mu(vals_inp,NULL,mu,mu_x0,mu_0y,dmu);
    //
    double ret_val = 0.0;

    if (by_individual) {
        ret_val = - arma::accu(2.0*mu_hat % arma::log(mu)) - arma::accu(mu_hat_x0 % arma::log(mu_x0)) - arma::accu(mu_hat_0y % arma::log(mu_0y));

        if (grad_vec) {
            arma::mat term_1 = arma::trans( elem_sub(2.0*mu_hat/arma::reshape(mu,nbX,nbY), mu_hat_x0/mu_x0) );
            arma::mat term_2 = mu_hat_0y / mu_0y;

            arma::mat term_grad = elem_prod(arma::vectorise(arma::trans(elem_sub(term_1,term_2))),dmu);

            *grad_vec = (- arma::trans(arma::sum(term_grad,0))) / scale;
        }
    } else {
        double N = arma::accu(arma::join_cols(arma::vectorise(mu),arma::join_cols(mu_x0,mu_0y)));

        ret_val = - arma::accu(mu_hat % arma::log(mu / N)) - arma::accu(mu_hat_x0 % arma::log(mu_x0 / N)) - arma::accu(mu_hat_0y % arma::log(mu_0y / N));

        if (grad_vec) {
            arma::mat term_1 = arma::trans( elem_sub(mu_hat/arma::reshape(mu,nbX,nbY), mu_hat_x0/mu_x0) );
            arma::mat term_2 = mu_hat_0y / mu_0y;

            arma::mat term_grad = elem_prod(arma::vectorise(arma::trans(elem_sub(term_1,term_2))),dmu);

            arma::vec term_3 = (arma::accu(arma::join_cols(arma::vectorise(mu_hat),arma::join_cols(mu_hat_x0,mu_hat_0y))) / N)  * arma::trans(arma::sum(dmu,0));

            *grad_vec = (- arma::trans(arma::sum(term_grad,0)) - term_3) / scale;
        }
    }
    //
    return ret_val / scale;
}

template<typename Tg, typename Th, typename Tt>
bool 
model<dse<Tg,Th,Tt>>::model_mme_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, double* value_out, double* err_tol_inp, int* max_iter_inp)
{
    bool success = generic_optim(init_out_vals,opt_objfn,opt_data,value_out,err_tol_inp,max_iter_inp);
    //
    return success;
}

template<typename Tg, typename Th, typename Tt>
double 
model<dse<Tg,Th,Tt>>::model_mme_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void* opt_data)
{
    trame_model_mme_opt_data<dse<Tg,Th,Tt>> *d = reinterpret_cast<trame_model_mme_opt_data<dse<Tg,Th,Tt>>*>(opt_data);
    //
    int nbX = d->market.nbX;
    int nbY = d->market.nbY;
    int dim_theta = d->dim_theta;
    arma::mat C_hat = d->C_hat;
    arma::mat kron_term = d->kron_term;
    //
    arma::mat U = arma::reshape(vals_inp.rows(0,nbX*nbY-1),nbX,nbY);

    arma::mat theta = vals_inp.rows(nbX*nbY,dim_theta + nbX*nbY - 1);
    arma::mat phi_mat = arma::reshape(kron_term * theta,nbX,nbY);
    //
    arma::mat mu_G, mu_H;

    double val_G = d->market.arums_G.G(d->market.n,U,mu_G);
    double val_H = d->market.arums_H.G(d->market.m,arma::trans(phi_mat - U),mu_H);
    //
    double ret = val_G + val_H - arma::accu(theta%C_hat);

    if (grad) {
        arma::vec grad_U = arma::vectorise(mu_G - mu_H.t());
        arma::vec grad_theta = arma::vectorise( arma::trans(arma::vectorise(mu_H.t())) * kron_term ) - C_hat;
        
        *grad = arma::join_cols(grad_U,grad_theta);
    }
    //
    return ret;
}

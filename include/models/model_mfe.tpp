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
 * model<mfe> class
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * This version:
 * 06/06/2017
 */

// Note: 'theta' refers to model parameters;'params' refers to the structural parameters

//
// build

template<typename Tt>
inline
void
model<mfe<Tt>>::build(const arma::mat& X_inp, const arma::mat& Y_inp)
{
    this->build_int(X_inp,Y_inp,NULL,NULL,NULL);
}

template<typename Tt>
inline
void
model<mfe<Tt>>::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    this->build_int(X_inp,Y_inp,&n_inp,&m_inp,NULL);
}

template<typename Tt>
inline
void
model<mfe<Tt>>::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp, const double& sigma_inp)
{
    this->build_int(X_inp,Y_inp,&n_inp,&m_inp,&sigma_inp);
}

template<typename Tt>
inline
void
model<mfe<Tt>>::build_int(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec* n_inp, const arma::vec* m_inp, const double* sigma_inp)
{
    need_norm = false;

    nbX = X_inp.n_rows;
    nbY = Y_inp.n_rows;

    dX = X_inp.n_cols;
    dY = Y_inp.n_cols;

    dim_theta = dX*dY;

    double sigma = (sigma_inp) ? *sigma_inp : 1.0;
    market_obj.sigma = sigma;
    //
    n = (n_inp) ? *n_inp : arma::ones(nbX,1);
    m = (m_inp) ? *m_inp : arma::ones(nbY,1);

    if (arma::accu(n) != arma::accu(m)) {
        printf("Unequal mass of individuals in an affinity model.\n");
    }
    //
    model_data = model_build_int(market_obj,X_inp,Y_inp);
}

//
// build markets

template<typename Tt>
void
model<mfe<Tt>>::model_to_market(const arma::mat& theta)
{
    model_to_market_int(market_obj,model_data,theta,n,m,nbX,nbY,dX,dY,need_norm);
}

//
// gradients

template<typename Tt>
void
model<mfe<Tt>>::dtheta(const arma::mat* delta_theta_inp, arma::mat& dtheta_M_out)
{
    arma::mat delta_theta_mat = (delta_theta_inp) ? *delta_theta_inp : arma::eye(dim_theta,dim_theta);

    dtheta_M_out = model_data * delta_theta_mat;
}

template<typename Tt>
arma::mat
model<mfe<Tt>>::dtheta(const arma::mat* delta_theta_inp)
{
    arma::mat dtheta_M_out;
    this->dtheta(delta_theta_inp,dtheta_M_out);

    return dtheta_M_out;
}

//
// MME

// template<typename Tt>
// bool 
// model<mfe<Tt>>::mme(const arma::mat& mu_hat, arma::mat& theta_hat, const arma::mat* theta_0_inp)
// {
//     bool success = false;
//     //
//     // double xtol_rel = (xtol_rel_inp) ? *xtol_rel_inp : 1E-04;
//     // int max_eval = (max_eval_inp) ? *max_eval_inp : 1E05;
//     double err_tol = 1E-06;
//     int max_iter = 5000;

//     arma::vec theta_0;
//     (theta_0_inp) ? theta_0 = *theta_0_inp : theta_0 = initial_theta();

//     arma::mat dtheta_Psi;
//     dtheta(NULL,dtheta_Psi);

//     model_to_market(theta_0);

//     arma::mat kron_term = dtheta_Psi;
//     arma::mat C_hat = arma::vectorise(arma::trans(arma::vectorise(mu_hat)) * kron_term);
//     //
//     // add optimization data
//     trame_model_mme_opt_data<mfe<Tt>> opt_data;
    
//     opt_data.market = market_obj;
//     opt_data.dim_theta = dim_theta;
//     opt_data.C_hat = C_hat;
//     opt_data.kron_term = kron_term;
//     //
//     arma::vec sol_vec = arma::join_cols(arma::vectorise(kron_term * theta_0)/2.0,theta_0);

//     double obj_val = 0;

//     success = model_mme_optim(sol_vec,model_mme_opt_objfn,&opt_data,&obj_val,&err_tol,&max_iter);
//     //
//     arma::mat U = arma::reshape(sol_vec.rows(0,nbX*nbY-1),nbX,nbY);
//     theta_hat = sol_vec.rows(nbX*nbY,nbX*nbY+dim_theta-1);

//     // double val_ret = obj_val;
//     //
//     return success;
// }

template<typename Tt>
bool
model<mfe<Tt>>::mme_woregul(const arma::mat& mu_hat, arma::mat& theta_hat, double& val_ret, double* xtol_rel_inp, int* max_eval_inp, double* tol_ipfp_inp, double* max_iter_ipfp_inp)
{
    bool success = false;
    //
    double xtol_rel = (xtol_rel_inp) ? *xtol_rel_inp : 1E-04;
    int max_eval = (max_eval_inp) ? *max_eval_inp : 1E05;

    double tol_ipfp = (tol_ipfp_inp) ? *tol_ipfp_inp : 1E-14;
    int max_iter_ipfp = (max_iter_ipfp_inp) ? *max_iter_ipfp_inp : 1E05;
    //
    arma::mat C_hat = Phi_k(mu_hat);
    //
    double total_mass = arma::accu(n);

    if (std::abs(arma::accu(mu_hat) - total_mass) > 1E-06) { // we use this instead of a raw '!=' to account for rounding error
        printf("Total number of couples does not coincide with margins.\n");
        return success;
    }
    //
    arma::vec p = n / total_mass;
    arma::vec q = m / total_mass;

    arma::mat IX = arma::ones(nbX,1);
    arma::mat tIY = arma::ones(1,nbY);

    arma::mat f = p * tIY;
    arma::mat g = IX * q.t();

    arma::mat Pi_hat = mu_hat / total_mass;
    arma::mat v = arma::zeros(1,nbY);

    // arma::mat phi_xy = arma::reshape(phi_xyk_aux,nbX*nbY,dim_theta);
    //
    // add optimization data
    trame_model_mfe_mme_opt_data opt_data;

    opt_data.max_iter_ipfp = max_iter_ipfp;
    opt_data.tol_ipfp = tol_ipfp;

    opt_data.nbX = nbX;
    opt_data.nbY = nbY;

    opt_data.dX = dX;
    opt_data.dY = dY;

    opt_data.sigma = market_obj.sigma;

    opt_data.p = p;
    opt_data.q = q;

    opt_data.IX = IX;
    opt_data.tIY = tIY;

    opt_data.f = f;
    opt_data.g = g;

    opt_data.v = v;
    opt_data.Pi_hat = Pi_hat;

    opt_data.phi_xy = model_data; // should be (nbX*nbY) x (dim_theta)
    //
    arma::vec opt_vec = arma::zeros(dX*dY,1);
    double obj_val = 0;

    success = model_mme_optim(opt_vec,model_mfe_mme_opt_objfn,&opt_data,&obj_val,&xtol_rel,&max_eval);
    //
    if (success) {
        theta_hat = opt_vec;
        val_ret = model_mfe_mme_opt_objfn(opt_vec,NULL,&opt_data);
    }
    //
    return success;
}

template<typename Tt>
bool
model<mfe<Tt>>::mme_regul(const arma::mat& mu_hat, const double& lambda, arma::mat& theta_hat, double& val_ret, double* xtol_rel_inp, int* max_eval_inp, double* tol_ipfp_inp, double* max_iter_ipfp_inp)
{
    bool success = false;
    //
    double xtol_rel = (xtol_rel_inp) ? *xtol_rel_inp : 1E-04;
    int max_eval = (max_eval_inp) ? *max_eval_inp : 1E05;

    double tol_ipfp = (tol_ipfp_inp) ? *tol_ipfp_inp : 1E-14;
    int max_iter_ipfp = (max_iter_ipfp_inp) ? *max_iter_ipfp_inp : 1E05;

    double sigma = market_obj.sigma;
    //
    arma::mat C_hat = Phi_k(mu_hat);
    // arma::mat C_hat = arma::vectorise(arma::trans(arma::vectorise(mu_hat)) * model_data);
    //
    double total_mass = arma::accu(n);

    if (std::abs(arma::accu(mu_hat) - total_mass) > 1E-06) { // we use this instead of a raw '!=' to account for rounding error
        printf("Total number of couples does not coincide with margins.\n");
        return success;
    }
    //
    arma::vec p = n / total_mass;
    arma::vec q = m / total_mass;

    arma::mat IX = arma::ones(nbX,1);
    arma::mat tIY = arma::ones(1,nbY);

    arma::mat f = p * tIY;
    arma::mat g = IX * q.t();

    arma::mat Pi_hat = mu_hat / total_mass;
    arma::mat v = arma::zeros(1,nbY);
    arma::mat A = arma::zeros(dX*dY,1);

    arma::mat Phi = arma::reshape(dtheta(&A),nbX,nbY);

    // arma::mat Phi = arma::reshape(Phi_xy(arma::vectorise(A)),nbX,nbY);
    //
    int iter_ipfp = 0, iter_count = 0;
    double err_ipfp = 2*tol_ipfp, err_val = 1.0;
    double t_k = 0.3; // step size for the prox grad algorithm (or grad descent when lambda=0)
    double alpha = 1.0; // for optimality check
    
    double the_val = 1.0, the_val_old = 1E04;
    arma::vec d, d_opt;
    arma::mat v_next = v, u, Pi, the_grad, A_mat, U, V, D, svd_mat, D_opt, opt_mat;

    while (err_val > xtol_rel && iter_count < max_eval) {
        iter_count++;

        Phi = arma::reshape(dtheta(&A),nbX,nbY);
        err_ipfp= 2*tol_ipfp;
        iter_ipfp = 0;

        while (err_ipfp > tol_ipfp && iter_ipfp < max_iter_ipfp) {
            iter_ipfp++;

            u = sigma * arma::log(arma::sum(g % arma::exp((Phi - IX * v)/sigma),1));
            v_next = sigma * arma::log(arma::sum(f % arma::exp((Phi - u * tIY)/sigma),0));

            err_ipfp = elem_max(arma::abs( arma::sum(g % arma::exp((Phi - IX * v_next - u*tIY)/sigma),1) - 1.0 ));
            v = v_next;
        }
        //
        Pi = f % g % arma::exp( (Phi - IX*v - u*tIY)/sigma );
        the_grad = Phi_k(Pi - Pi_hat);

        A -= t_k*the_grad;
        //
        if (lambda > 0) {
            // compute the proximal operator
            A_mat = arma::reshape(A,dX,dY);
            arma::svd(U,d,V,A_mat);
            D = arma::diagmat(elem_max(d - lambda*t_k,0.0));
            A = arma::vectorise(U * D * V.t());
        } // if lambda = 0 then we are just taking one step of gradient descent
        //
        if (iter_count % 10 == 0) {
            //alpha = 1.0;
            svd_mat = arma::reshape(A - alpha*the_grad,dX,dY);
            arma::svd(U,d_opt,V,svd_mat);
            D_opt = arma::diagmat(elem_max(d_opt - alpha*lambda, 0.0));

            opt_mat = arma::accu(arma::pow(A - arma::vectorise(U * D_opt *V.t()), 2));
        }
        //
        if (lambda > 0) {
            the_val = arma::accu(the_grad % arma::vectorise(A)) - sigma * arma::accu(Pi % arma::log(Pi)) + lambda * arma::accu(D);
        } else {
            the_val = arma::accu(the_grad % arma::vectorise(A)) - sigma * arma::accu(Pi % arma::log(Pi));
        }

        err_val = std::abs(the_val - the_val_old);
        //
        the_val_old = the_val;
    }

    if (err_val <= xtol_rel && iter_count < max_eval) {
        success = true;
    }
    //
    theta_hat = arma::vectorise(A);
    val_ret = the_val;
    //
    return success;
}

//
// solve wrappers

template<typename Tt>
bool
model<mfe<Tt>>::solve(arma::mat& mu_sol)
{
    bool res = market_obj.solve(mu_sol,NULL);
    //
    return res;
}

template<typename Tt>
bool
model<mfe<Tt>>::solve(arma::mat& mu_sol, const char* solver)
{
    bool res = market_obj.solve(mu_sol,solver);
    //
    return res;
}

template<typename Tt>
bool
model<mfe<Tt>>::solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver)
{
    bool res = market_obj.solve(mu_sol,U,V,solver);
    //
    return res;
}

//
// internal

template<typename Tt>
inline
void
model<mfe<Tt>>::initial_theta(arma::mat& theta_0)
{
    theta_0.zeros(dim_theta,1);
}

template<typename Tt>
inline
arma::mat
model<mfe<Tt>>::initial_theta()
{
    return arma::zeros(dim_theta,1);
}

template<typename Tt>
inline
arma::mat
model<mfe<Tt>>::Phi_k(const arma::mat& mu_hat)
{
    arma::mat ret = arma::vectorise(arma::trans(arma::vectorise(mu_hat))*model_data);
    //
    return ret;
}

//
// optimization-related functions

template<typename Tt>
bool
model<mfe<Tt>>::model_mme_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, double* value_out, double* err_tol_inp, int* max_iter_inp)
{
    bool success = generic_optim(init_out_vals,opt_objfn,opt_data,value_out,err_tol_inp,max_iter_inp);
    //
    return success;
}

template<typename Tt>
double
model<mfe<Tt>>::model_mfe_mme_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void* opt_data)
{
    trame_model_mfe_mme_opt_data *d = reinterpret_cast<trame_model_mfe_mme_opt_data*>(opt_data);
    //
    int nbX = d->nbX;
    int nbY = d->nbY;

    //int dX = d->dX;
    //int dY = d->dY;
    
    int max_iter_ipfp = d->max_iter_ipfp;
    double tol_ipfp = d->tol_ipfp;

    double sigma = d->sigma;

    arma::vec p = d->p;
    arma::vec q = d->q;

    arma::mat IX = d->IX;
    arma::mat tIY = d->tIY;

    arma::mat f = d->f;
    arma::mat g = d->g;

    arma::mat v = d->v;
    arma::mat Pi_hat = d->Pi_hat;

    arma::mat phi_xy = d->phi_xy; // should be (nbX*nbY) x (dim_theta), replaces a member function call
    //
    arma::mat Phi = arma::reshape(phi_xy * vals_inp,nbX,nbY);
    //
    int iter_ipfp = 0;
    double err_ipfp= 2*tol_ipfp;
    arma::mat v_next = v, u;

    while (err_ipfp > tol_ipfp && iter_ipfp < max_iter_ipfp) {
        iter_ipfp++;

        u = sigma * arma::log(arma::sum(g % arma::exp((Phi - IX * v)/sigma),1));
        v_next = sigma * arma::log(arma::sum( f % arma::exp((Phi - u * tIY)/sigma),0));

        err_ipfp = elem_max(arma::abs( arma::sum(g % arma::exp((Phi - IX * v_next - u*tIY)/sigma),1) - 1.0 ));
        v = v_next;
    }
    //
    arma::mat Pi = f % g % arma::exp( (Phi - IX*v - u*tIY)/sigma );
    arma::mat the_grad = phi_xy.t() * arma::vectorise(Pi - Pi_hat);

    if (grad) {
        *grad = the_grad;
    }
    //
    // update v for the next opt call
    d->v = v;
    opt_data = reinterpret_cast<void*>(d);
    //
    double ret = arma::accu(the_grad % arma::vectorise(vals_inp)) - sigma*arma::accu(Pi%arma::log(Pi));
    //
    return ret;
}

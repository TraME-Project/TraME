/*################################################################################
  ##
  ##   Copyright (C) 2015-2018 the TraME Team:
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
 * 05/07/2018
 */

// Notation: 'theta' refers to model parameters; 'params' refers to the structural parameters

//
// build

template<typename Tt>
void
model<mfe<Tt>>::build(const arma::mat& X_inp, const arma::mat& Y_inp)
{
    this->build_int(X_inp,Y_inp,nullptr,nullptr,nullptr);
}

template<typename Tt>
void
model<mfe<Tt>>::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    this->build_int(X_inp,Y_inp,&n_inp,&m_inp,nullptr);
}

template<typename Tt>
void
model<mfe<Tt>>::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp, const double sigma_inp)
{
    this->build_int(X_inp,Y_inp,&n_inp,&m_inp,&sigma_inp);
}

template<typename Tt>
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

    // types

    n = (n_inp) ? *n_inp : arma::ones(nbX,1);
    m = (m_inp) ? *m_inp : arma::ones(nbY,1);

    if (arma::accu(n) != arma::accu(m))
    {
        printf("Unequal mass of individuals in an affinity model.\n");
    }

    // construct model 'data'

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
    dtheta_M_out = (delta_theta_inp) ? model_data * (*delta_theta_inp) : model_data;
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

template<typename Tt>
bool
model<mfe<Tt>>::mme_woregul(const arma::mat& mu_hat, arma::mat& theta_hat)
{
    return this->mme_woregul(mu_hat,theta_hat,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tt>
bool
model<mfe<Tt>>::mme_woregul(const arma::mat& mu_hat, arma::mat& theta_hat, double& val_ret)
{
    return this->mme_woregul(mu_hat,theta_hat,&val_ret,nullptr,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tt>
bool
model<mfe<Tt>>::mme_woregul(const arma::mat& mu_hat, arma::mat& theta_hat, double* val_ret, 
                            double* xtol_rel_inp, int* max_eval_inp, 
                            double* tol_ipfp_inp, int* max_iter_ipfp_inp, 
                            const int* optim_method_inp)
{
    bool success = false;

    //

    const double xtol_rel = (xtol_rel_inp) ? *xtol_rel_inp : 1E-04;
    const int max_eval = (max_eval_inp) ? *max_eval_inp : 1E05;

    const double tol_ipfp = (tol_ipfp_inp) ? *tol_ipfp_inp : 1E-12;
    const int max_iter_ipfp = (max_iter_ipfp_inp) ? *max_iter_ipfp_inp : 1E05;

    //

    const double total_mass = arma::accu(n);

    if (std::abs(arma::accu(mu_hat) - total_mass) > 1E-06) { // we use this instead of a raw '!=' to account for rounding error
        printf("Total number of couples does not coincide with margins.\n");
        return success;
    }

    //

    arma::mat C_hat = Phi_k(mu_hat);

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

    opt_data.IX = IX;
    opt_data.tIY = tIY;

    opt_data.f = f;
    opt_data.g = g;

    opt_data.v = v;
    opt_data.Pi_hat = Pi_hat;

    opt_data.phi_xy = model_data; // should be (nbX*nbY) x (dim_theta)

    // optim setup

    const int optim_method = (optim_method_inp) ? *optim_method_inp : 1;

    optim::algo_settings_t settings;

    settings.err_tol = xtol_rel;
    settings.iter_max = max_eval;

    arma::vec opt_vec = arma::zeros(dX*dY,1);

    success = model_mme_optim(opt_vec,model_mfe_mme_opt_objfn,&opt_data,&settings,optim_method);

    //

    if (success)
    {
        theta_hat = opt_vec;

        if (val_ret) {
            *val_ret = settings.opt_value;
        }
    }

    return success;
}

template<typename Tt>
bool
model<mfe<Tt>>::mme_regul(const arma::mat& mu_hat, arma::mat& theta_hat, const double lambda)
{
    return this->mme_regul(mu_hat,theta_hat,lambda,nullptr,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tt>
bool
model<mfe<Tt>>::mme_regul(const arma::mat& mu_hat, arma::mat& theta_hat, const double lambda, double& val_ret)
{
    return this->mme_regul(mu_hat,theta_hat,lambda,&val_ret,nullptr,nullptr,nullptr,nullptr);
}

template<typename Tt>
bool
model<mfe<Tt>>::mme_regul(const arma::mat& mu_hat, arma::mat& theta_hat, const double lambda, double* val_ret, 
                          double* xtol_rel_inp, int* max_eval_inp, double* tol_ipfp_inp, int* max_iter_ipfp_inp)
{
    bool success = false;

    //

    const double xtol_rel = (xtol_rel_inp) ? *xtol_rel_inp : 1E-04;
    const int max_eval = (max_eval_inp) ? *max_eval_inp : 1E05;

    const double tol_ipfp = (tol_ipfp_inp) ? *tol_ipfp_inp : 1E-12;
    const int max_iter_ipfp = (max_iter_ipfp_inp) ? *max_iter_ipfp_inp : 1E05;

    const double sigma = market_obj.sigma;

    // check total mass condition

    double total_mass = arma::accu(n);

    if (std::abs(arma::accu(mu_hat) - total_mass) > 1E-06)
    {   // we use this instead of a raw '!=' to account for rounding error
        printf("Total number of couples does not coincide with margins.\n");
        return success;
    }

    //

    const arma::vec p = n / total_mass;
    const arma::vec q = m / total_mass;

    const arma::mat IX = arma::ones(nbX,1);
    const arma::mat tIY = arma::ones(1,nbY);

    const arma::mat f = p * tIY;
    const arma::mat g = IX * q.t();

    const arma::mat Pi_hat = mu_hat / total_mass;
    arma::rowvec v = arma::zeros(1,nbY);
    arma::vec A = arma::zeros(dX*dY,1);

    //

    int iter_ipfp = 0, iter_count = 0;
    double err_ipfp = 2*tol_ipfp, err_val = 1.0;

    const double t_k = 0.3;   // step size for the prox. grad algorithm (or grad descent when lambda=0)
    const double alpha = 1.0; // for optimality check

    double the_val = 1.0, the_val_old = 1E04;
    double opt_check = 1.0;

    arma::vec d, d_opt, u;
    arma::rowvec v_next = v;
    arma::mat U, V;

    while (err_val > xtol_rel && iter_count < max_eval)
    {
        iter_count++;

        arma::mat Phi = arma::reshape(dtheta(&A),nbX,nbY);

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

        arma::mat Pi = f % g % arma::exp( (Phi - IX*v - u*tIY)/sigma );
        arma::mat the_grad = Phi_k(Pi - Pi_hat);

        A -= t_k*the_grad;

        //

        if (lambda > 0.0)
        {   // compute the proximal operator
            arma::mat A_mat = arma::reshape(A,dX,dY);
            arma::svd(U,d,V,A_mat);

            arma::mat D = arma::diagmat(elem_max(d - lambda*t_k,0.0));

            if (dX != dY) 
            {   // R's SVD behaves differently
                const int d_l = d.n_elem;
                A = arma::vectorise(U.cols(0,d_l-1) * D * V.cols(0,d_l-1).t());
            }
            else
            {
                A = arma::vectorise(U * D * V.t());
            }

            // arma::cout << "A mfe up\n" << A << arma::endl;

            the_val = arma::accu(the_grad % arma::vectorise(A)) - sigma * arma::accu(Pi % arma::log(Pi)) + lambda * arma::accu(D);
        }
        else
        {   // if lambda = 0 then we are just taking one step of gradient descent
            the_val = arma::accu(the_grad % arma::vectorise(A)) - sigma * arma::accu(Pi % arma::log(Pi));
        }

        //

        if (iter_count % 10 == 0)
        {

            arma::mat svd_mat = arma::reshape(A - alpha*the_grad,dX,dY);
            arma::svd(U,d_opt,V,svd_mat);

            arma::mat D_opt = arma::diagmat(elem_max(d_opt - alpha*lambda, 0.0));
            const int d_l = d_opt.n_elem;

            // opt_check = arma::accu(arma::pow(A - arma::vectorise(U.cols(0,d_l-1) * D_opt *V.cols(0,d_l-1).t()), 2));
            opt_check = std::pow( arma::norm(A - arma::vectorise(U.cols(0,d_l-1) * D_opt * V.cols(0,d_l-1).t()),2) , 2);

            std::cout << "Regularized MME: testing optimality: " << opt_check << std::endl;
        }

        // conv check

        if (iter_count > 1) {
            err_val = std::abs(the_val - the_val_old);
        }

        the_val_old = the_val;
    }

    if (err_val <= xtol_rel && iter_count < max_eval)
    {
        success = true;
    }

    //

    // theta_hat = arma::vectorise(A);
    theta_hat = A;

    if (val_ret) {
        *val_ret = the_val;
    }

    //
    return success;
}

//
// solve wrappers

template<typename Tt>
bool
model<mfe<Tt>>::solve(arma::mat& mu_sol)
{
    return market_obj.solve(mu_sol,nullptr);
}

template<typename Tt>
bool
model<mfe<Tt>>::solve(arma::mat& mu_sol, const char* solver)
{
    return market_obj.solve(mu_sol,solver);
}

template<typename Tt>
bool
model<mfe<Tt>>::solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver)
{
    return market_obj.solve(mu_sol,U,V,solver);
}

//
// internal

template<typename Tt>
void
model<mfe<Tt>>::initial_theta(arma::mat& theta_0)
{
    theta_0.zeros(dim_theta,1);
}

template<typename Tt>
arma::mat
model<mfe<Tt>>::initial_theta()
{
    return arma::zeros(dim_theta,1);
}

template<typename Tt>
arma::mat
model<mfe<Tt>>::Phi_k(const arma::mat& mu_hat)
{
    return arma::vectorise(arma::trans(arma::vectorise(mu_hat))*model_data);
}

//
// optimization-related functions

template<typename Tt>
bool
model<mfe<Tt>>::model_mme_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn,
                                void* opt_data, optim::algo_settings_t* settings_inp, const int optim_method)
{
    if (optim_method == 1)
    {
        return optim::lbfgs_int(init_out_vals,opt_objfn,opt_data,settings_inp);
    }
    else if (optim_method == 2)
    {
        return optim::bfgs_int(init_out_vals,opt_objfn,opt_data,settings_inp);
    }
    else
    {
        printf("error: unrecognized optim_method choice.\n");
        return false;
    }
}

template<typename Tt>
double
model<mfe<Tt>>::model_mfe_mme_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void* opt_data)
{
    trame_model_mfe_mme_opt_data *d = reinterpret_cast<trame_model_mfe_mme_opt_data*>(opt_data);

    //

    const int nbX = d->nbX;
    const int nbY = d->nbY;

    const int max_iter_ipfp = d->max_iter_ipfp;
    const double tol_ipfp = d->tol_ipfp;

    const double sigma = d->sigma;

    const arma::mat f = d->f;
    const arma::mat g = d->g;

    arma::rowvec v = d->v;

    arma::mat phi_xy = d->phi_xy; // should be (nbX*nbY) x (dim_theta), replaces a member function call

    //

    const arma::mat Phi = arma::reshape(phi_xy * vals_inp,nbX,nbY);

    //

    int iter_ipfp = 0;
    double err_ipfp= 2*tol_ipfp;

    arma::rowvec v_next = v;
    arma::vec u;

    while (err_ipfp > tol_ipfp && iter_ipfp < max_iter_ipfp)
    {
        iter_ipfp++;

        u = sigma * arma::log(arma::sum(g % arma::exp((Phi - d->IX * v)/sigma),1));
        v_next = sigma * arma::log(arma::sum( f % arma::exp((Phi - u * d->tIY)/sigma),0));

        err_ipfp = elem_max(arma::abs( arma::sum(g % arma::exp((Phi - d->IX * v_next - u * d->tIY)/sigma),1) - 1.0 ));
        v = v_next;
    }

    //

    const arma::mat Pi = f % g % arma::exp( (Phi - d->IX * v - u * d->tIY)/sigma );
    const arma::mat the_grad = phi_xy.t() * arma::vectorise(Pi - d->Pi_hat);

    if (grad)
    {
        *grad = the_grad;
    }

    // update v for the next opt call

    d->v = v;
    // opt_data = reinterpret_cast<void*>(d);

    //

    return arma::accu(the_grad % arma::vectorise(vals_inp)) - sigma*arma::accu(Pi%arma::log(Pi));
}

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
 * model<dse> class
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * This version:
 * 07/26/2017
 */

// Note: 'theta' refers to model parameters; 'params' refers to structural parameters

//
// first method to build

template<typename Tg, typename Th, typename Tt>
inline
void
model<dse<Tg,Th,Tt>>::build(const arma::cube& phi_xyk_inp)
{
    this->build_int(phi_xyk_inp,nullptr,nullptr);
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

    model_data = cube_to_mat(phi_xyk_inp); // Phi_xy = matrix(phi_xyk)
}

//
// second method to build

template<typename Tg, typename Th, typename Tt>
inline
void
model<dse<Tg,Th,Tt>>::build(const arma::mat& X_inp, const arma::mat& Y_inp)
{
    this->build_int(X_inp,Y_inp,nullptr,nullptr);
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
}

//
// build markets

template<typename Tg, typename Th, typename Tt>
void
model<dse<Tg,Th,Tt>>::model_to_market(const arma::mat& theta)
{
    model_to_market_int(market_obj,model_data,theta,n,m,nbX,nbY,dX,dY,need_norm);
}

template<typename Tg, typename Th, typename Tt>
void
model<dse<Tg,Th,Tt>>::model_to_market(const arma::mat& theta, const Tg& arums_G_inp, const Th& arums_H_inp)
{
    model_to_market_int(market_obj,model_data,theta,arums_G_inp,arums_H_inp,n,m,nbX,nbY,dX,dY,need_norm);
}

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
    this->dtheta(delta_theta_inp,dtheta_Psi_out,nullptr,nullptr);
}

template<typename Tg, typename Th, typename Tt>
void
model<dse<Tg,Th,Tt>>::dtheta(const arma::mat* delta_theta_inp, arma::mat& dtheta_Psi_out, arma::mat* dtheta_G_out, arma::mat* dtheta_H_out)
{
    const arma::mat delta_theta_mat = (delta_theta_inp) ? *delta_theta_inp : arma::eye(dim_theta,dim_theta);

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
    dtheta(delta_theta,dtheta_Psi,nullptr,nullptr);
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
    const double err_tol = 1E-06;
    const int max_iter = 5000;

    arma::vec theta_0;
    (theta_0_inp) ? theta_0 = *theta_0_inp : theta_0 = initial_theta();

    model_to_market(theta_0);

    const bool by_individual = true;
    const double scale = std::max(arma::accu(n),arma::accu(m));

    const arma::vec mu_hat_x0 = n - arma::sum(mu_hat,1);
    const arma::vec mu_hat_0y = m - arma::trans(arma::sum(mu_hat,0));
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
model<dse<Tg,Th,Tt>>::mme(const arma::mat& mu_hat, arma::mat& theta_hat, const arma::mat* theta_0_inp)
{
    bool success = false;
    //
    const double err_tol = 1E-04;
    const int max_iter = 1000;

    arma::vec theta_0(dim_theta);
    (theta_0_inp) ? theta_0 = *theta_0_inp : theta_0 = initial_theta();

    arma::mat dtheta_Psi;
    dtheta(nullptr,dtheta_Psi);

    model_to_market(theta_0);

    const arma::mat kron_term = dtheta_Psi;
    const arma::mat C_hat = arma::vectorise(arma::trans(arma::vectorise(mu_hat)) * kron_term);
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
    // arma::mat U = arma::reshape(sol_vec.rows(0,nbX*nbY-1),nbX,nbY);
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
    return market_obj.solve(mu_sol,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
model<dse<Tg,Th,Tt>>::solve(arma::mat& mu_sol, const char* solver)
{
    return market_obj.solve(mu_sol,solver);
}

template<typename Tg, typename Th, typename Tt>
bool
model<dse<Tg,Th,Tt>>::solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver)
{
    return market_obj.solve(mu_sol,U,V,solver);
}

//
// internal

template<typename Tg, typename Th, typename Tt>
inline
void
model<dse<Tg,Th,Tt>>::initial_theta(arma::mat& theta_0)
{
    theta_0.zeros(dim_theta,1);
}

template<typename Tg, typename Th, typename Tt>
inline
arma::mat
model<dse<Tg,Th,Tt>>::initial_theta()
{
    return arma::zeros(dim_theta,1);
}

//
// optimization-related functions
//

//
// MLE functions

template<typename Tg, typename Th, typename Tt>
bool
model<dse<Tg,Th,Tt>>::model_mle_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, double* value_out, double* err_tol_inp, int* max_iter_inp)
{
    optim::optim_opt_settings opt_params;

    if (err_tol_inp) {
        opt_params.err_tol = *err_tol_inp;
    }

    if (max_iter_inp) {
        opt_params.iter_max = *max_iter_inp;
    }

    return optim::generic_optim_int(init_out_vals,opt_objfn,opt_data,value_out,&opt_params);
}

template<typename Tg, typename Th, typename Tt>
double
model<dse<Tg,Th,Tt>>::log_likelihood(const arma::vec& vals_inp, arma::vec* grad_vec, void* opt_data)
{
    trame_model_mle_opt_data<dse<Tg,Th,Tt>> *d = reinterpret_cast<trame_model_mle_opt_data<dse<Tg,Th,Tt>>*>(opt_data);

    const bool by_individual = d->by_individual;

    const int nbX = d->model_obj.nbX;
    const int nbY = d->model_obj.nbY;

    const double scale = d->scale;

    const arma::mat mu_hat = d->mu_hat;
    const arma::vec mu_hat_x0 = d->mu_hat_x0;
    const arma::vec mu_hat_0y = d->mu_hat_0y;
    //
    arma::vec mu_x0, mu_0y;
    arma::mat mu, dmu;

    d->model_obj.dtheta_mu(vals_inp,nullptr,mu,mu_x0,mu_0y,dmu);
    //
    double ret_val = 0.0;

    if (by_individual) {
        ret_val = - arma::accu(2.0*mu_hat % arma::log(mu)) - arma::accu(mu_hat_x0 % arma::log(mu_x0)) - arma::accu(mu_hat_0y % arma::log(mu_0y));

        if (grad_vec) {
            const arma::mat term_1 = arma::trans( elem_sub(2.0*mu_hat/arma::reshape(mu,nbX,nbY), mu_hat_x0/mu_x0) );
            const arma::mat term_2 = mu_hat_0y / mu_0y;

            const arma::mat term_grad = elem_prod(arma::vectorise(arma::trans(elem_sub(term_1,term_2))),dmu);

            *grad_vec = (- arma::trans(arma::sum(term_grad,0))) / scale;
        }
    } else {
        const double N = arma::accu(arma::join_cols(arma::vectorise(mu),arma::join_cols(mu_x0,mu_0y)));

        ret_val = - arma::accu(mu_hat % arma::log(mu / N)) - arma::accu(mu_hat_x0 % arma::log(mu_x0 / N)) - arma::accu(mu_hat_0y % arma::log(mu_0y / N));

        if (grad_vec) {
            const arma::mat term_1 = arma::trans( elem_sub(mu_hat/arma::reshape(mu,nbX,nbY), mu_hat_x0/mu_x0) );
            const arma::mat term_2 = mu_hat_0y / mu_0y;

            const arma::mat term_grad = elem_prod(arma::vectorise(arma::trans(elem_sub(term_1,term_2))),dmu);

            const arma::vec term_3 = (arma::accu(arma::join_cols(arma::vectorise(mu_hat),arma::join_cols(mu_hat_x0,mu_hat_0y))) / N)  * arma::trans(arma::sum(dmu,0));

            *grad_vec = (- arma::trans(arma::sum(term_grad,0)) - term_3) / scale;
        }
    }
    //
    return ret_val / scale;
}

//
// MME functions

template<typename Tg, typename Th, typename Tt>
bool
model<dse<Tg,Th,Tt>>::model_mme_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, double* value_out, double* err_tol_inp, int* max_iter_inp)
{
    optim::optim_opt_settings opt_params;

    if (err_tol_inp) {
        opt_params.err_tol = *err_tol_inp;
    }

    if (max_iter_inp) {
        opt_params.iter_max = *max_iter_inp;
    }

    return optim::generic_optim_int(init_out_vals,opt_objfn,opt_data,value_out,&opt_params);
}

template<typename Tg, typename Th, typename Tt>
double
model<dse<Tg,Th,Tt>>::model_mme_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void* opt_data)
{
    trame_model_mme_opt_data<dse<Tg,Th,Tt>> *d = reinterpret_cast<trame_model_mme_opt_data<dse<Tg,Th,Tt>>*>(opt_data);
    //
    const int nbX = d->market.nbX;
    const int nbY = d->market.nbY;
    const int dim_theta = d->dim_theta;

    const arma::mat C_hat = d->C_hat;
    const arma::mat kron_term = d->kron_term;
    //
    const arma::mat U = arma::reshape(vals_inp.rows(0,nbX*nbY-1),nbX,nbY);

    const arma::mat theta = vals_inp.rows(nbX*nbY,dim_theta + nbX*nbY - 1);
    const arma::mat phi_mat = arma::reshape(kron_term * theta,nbX,nbY);
    //
    arma::mat mu_G, mu_H;

    const double val_G = d->market.arums_G.G(d->market.n,U,mu_G);
    const double val_H = d->market.arums_H.G(d->market.m,arma::trans(phi_mat - U),mu_H);
    //
    const double ret = val_G + val_H - arma::accu(theta%C_hat);

    if (grad) {
        const arma::vec grad_U = arma::vectorise(mu_G - mu_H.t());
        const arma::vec grad_theta = arma::vectorise( arma::trans(arma::vectorise(mu_H.t())) * kron_term ) - C_hat;

        *grad = arma::join_cols(grad_U,grad_theta);
    }
    //
    return ret;
}

//
// MME for empirical class

template<>
inline
bool
model<dse<arums::empirical,arums::empirical,transfers::tu>>::mme(const arma::mat& mu_hat, arma::mat& theta_hat, double* val_out, arma::mat* mu_out, arma::mat* U_out, arma::mat* V_out)
{
    bool success = false;
    //
    const arma::mat kron_mat = model_data;
    const arma::mat kron_mat_2 = arma::reshape(kron_mat.t(),dim_theta*nbX,nbY);

    const arma::vec C_hat = arma::vectorise(arma::vectorise(mu_hat)*kron_mat);
    //
    arma::mat epsilon_iy, epsilon0_i, I_ix;
    arma::mat eta_xj, eta_0j, I_yj;

    const int nbDraws_1 = build_disaggregate_epsilon(n,market_obj.arums_G,epsilon_iy,epsilon0_i,I_ix);
    const int nbDraws_2 = build_disaggregate_epsilon(m,market_obj.arums_H,eta_xj,eta_0j,I_yj);

    epsilon0_i = arma::vectorise(epsilon0_i);

    eta_xj = eta_xj.t();
    eta_0j = arma::vectorise(eta_0j);
    I_yj = I_yj.t();
    
    const arma::vec n_i = arma::vectorise(I_ix * n) / (double) nbDraws_1;
    const arma::vec m_j = arma::vectorise(m.t() * I_yj) / (double) nbDraws_2;

    const int nbI = n_i.n_elem;
    const int nbJ = m_j.n_elem;

    arma::mat kron_data_mat = - arma::reshape(kron_mat_2*I_yj,dim_theta,nbX*nbJ);
    /*
     * use batch allocation to construct the sparse constraint matrix (A)
     *
     * A_sp_t size: (nbI + nbJ + nbX*nbY + dim_theta) x (nbI*nbY + nbJ*nbX)
     *
     * first block involves nbY blocks of nbI diagonal matrices, so nbI x (nbY*nbI) 
     *
     * second block begins on row nbI+1 and ends on row nbI+nbJ (inclusive); from column 
     * 1 to nbI*nbY nothing but zeros; then nbJ blocks of nbX-length row vectors of ones 
     *
     * third block begins on row nbI+nbJ+1 and ends on nbI+nbJ+nbX*nbY; 
     * there are nbX*nbY blocks of length nbDraws_1 from columns (1,nbI*nbY);
     * for columns (nbI+1)
     *
     * fourth block is filled from (nbI*nbY+1,nbI*nbY + nbJ*nbX) with kron_data_mat
     */
    int jj, kk, ll, count_val = 0;
    int num_nonzero_elem = nbI*nbY + nbJ*nbX + nbX*nbY*nbDraws_1 + nbY*nbDraws_2*nbX + dim_theta*nbJ*nbX;

    arma::umat location_mat_1(2,num_nonzero_elem);
    arma::rowvec vals_mat_1(num_nonzero_elem);

    for (jj=0; jj < nbY; jj++) { // first block
        for (kk=0; kk < nbI; kk++) {
            location_mat_1(0,count_val) = kk; // rows
            location_mat_1(1,count_val) = kk + jj*nbI; // columns

            vals_mat_1(count_val) = 1.0;

            ++count_val;
        }
    }

    for (jj=0; jj < nbJ; jj++) { // second block
        for (kk=0; kk < nbX; kk++) {
            location_mat_1(0,count_val) = nbI + jj;
            location_mat_1(1,count_val) = nbI*nbY + kk + jj*nbX;

            vals_mat_1(count_val) = 1.0;

            ++count_val;
        }
    }

    for (jj=0; jj < nbX*nbY; jj++) { // third block, part 1
        for (kk=0; kk < nbDraws_1; kk++) {
            location_mat_1(0,count_val) = nbI + nbJ + jj;
            location_mat_1(1,count_val) = kk + jj*nbDraws_1;

            vals_mat_1(count_val) = -1.0;

            ++count_val;
        }
    }

    for (jj=0; jj < nbY; jj++) { // third block, part 2
        for (kk=0; kk < nbDraws_2; kk++) {
            for (ll=0; ll < nbX; ll++) {
                location_mat_1(0,count_val) = nbI + nbJ + jj*nbX + ll;
                location_mat_1(1,count_val) = nbI*nbY + jj*(nbDraws_2*nbX) + kk*nbX + ll;

                vals_mat_1(count_val) = 1.0;

                ++count_val;
            }
        }
    }

    for (jj=0; jj < dim_theta; jj++) { // fourth block, data
        for (kk=0; kk < nbJ*nbX; kk++) {
            location_mat_1(0,count_val) = nbI + nbJ + nbX*nbY + jj;
            location_mat_1(1,count_val) = nbI*nbY + kk;

            vals_mat_1(count_val) = kron_data_mat(jj,kk);

            ++count_val;
        }
    }
    //
    // now proceed to solve the LP problem
    arma::sp_mat A_sp_t(location_mat_1,vals_mat_1); // transpose of A

    const int k_lp = A_sp_t.n_cols; // cols as we're working with the transpose
    const int n_lp = A_sp_t.n_rows; // rows as we're working with the transpose

    const arma::uword* row_vals = &(*A_sp_t.row_indices);
    const arma::uword* col_vals = &(*A_sp_t.col_ptrs);

    int* vind_lp = new int[num_nonzero_elem];
    int* vbeg_lp = new int[k_lp+1];
    double* vval_lp = new double[num_nonzero_elem];

    for (jj=0; jj<num_nonzero_elem; jj++) {
        vind_lp[jj] = row_vals[jj];
        vval_lp[jj] = A_sp_t.values[jj];
    }

    for (jj=0; jj<k_lp+1; jj++) {
        vbeg_lp[jj] = col_vals[jj];
    }
    //
    char* sense_lp = new char[k_lp];
    for (jj=0; jj<k_lp; jj++) {
        sense_lp[jj] = '>';
    }

    arma::vec lb_lp(epsilon0_i.n_elem + eta_0j.n_elem + nbX*nbY + dim_theta);
    lb_lp.rows(0,epsilon0_i.n_elem-1) = arma::vectorise(epsilon0_i);
    lb_lp.rows(epsilon0_i.n_elem,epsilon0_i.n_elem + eta_0j.n_elem - 1) = eta_0j;
    lb_lp.rows(epsilon0_i.n_elem + eta_0j.n_elem, lb_lp.n_rows - 1).fill(-arma::datum::inf);

    arma::vec rhs_lp(epsilon_iy.n_elem + eta_xj.n_elem);
    rhs_lp.rows(0,epsilon_iy.n_elem-1) = arma::vectorise(epsilon_iy);
    rhs_lp.rows(epsilon_iy.n_elem,rhs_lp.n_elem-1) = arma::vectorise(eta_xj);

    arma::vec obj_lp(nbI + nbJ + nbX*nbY + C_hat.n_elem);
    obj_lp.rows(0,nbI-1) = n_i;
    obj_lp.rows(nbI, nbI + nbJ - 1) = m_j;
    obj_lp.rows(nbI + nbJ, nbI + nbJ + nbX*nbY - 1).zeros();
    obj_lp.rows(nbI + nbJ + nbX*nbY, obj_lp.n_elem-1) = - C_hat;

    int modelSense = 0; // minimize
    
    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    bool lp_optimal = false;
    double val_lp = 0.0;
    //
    try {
        lp_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), num_nonzero_elem, vbeg_lp, vind_lp, vval_lp, modelSense, rhs_lp.memptr(), sense_lp, nullptr, lb_lp.memptr(), nullptr, nullptr, val_lp, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));
        
        if (lp_optimal) {

            theta_hat = sol_mat(arma::span(nbI+nbJ+nbX*nbY,nbI+nbJ+nbX*nbY+dim_theta-1),0);

            //
            // package up solution

            if (mu_out) {
                const arma::mat mu_iy = arma::reshape(dual_mat(arma::span(0,nbI*nbY-1),0),nbI,nbY);
                *mu_out = I_ix.t() * mu_iy;
            }

            if (U_out && V_out) {
                *U_out = arma::reshape(sol_mat(arma::span(nbI+nbJ,nbI+nbJ+nbX*nbY-1),0),nbX,nbY);
                *V_out = arma::reshape(kron_mat*theta_hat,nbX,nbY) - *U_out;
            } else if (U_out) {
                *U_out = arma::reshape(sol_mat(arma::span(nbI+nbJ,nbI+nbJ+nbX*nbY-1),0),nbX,nbY);
            } else if (V_out) {
                *V_out = arma::reshape(kron_mat*theta_hat,nbX,nbY) - arma::reshape(sol_mat(arma::span(nbI+nbJ,nbI+nbJ+nbX*nbY-1),0),nbX,nbY);
            }

            if (val_out) {
                *val_out = val_lp;
            }
            //
            success = true;
        } else {
            std::cout << "Non-optimal value found during optimization" << std::endl;
        }
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //
    return success;
}

//
// marp proj

template<>
inline
bool
model<dse<arums::none,arums::none,transfers::tu>>::mme(const arma::mat& mu_hat, arma::mat& theta_hat, const arma::mat* theta_0_inp)
{
    bool success = false;
    //
    const arma::mat kron_term = model_data;
    const arma::mat C_hat = arma::vectorise(arma::trans(arma::vectorise(mu_hat)) * kron_term);
    //
    arma::vec obj_lp = arma::join_cols(n,arma::join_cols(m,arma::zeros(dim_theta,1)));

    arma::mat A_11_lp = arma::kron(arma::ones(nbY,1), arma::eye(nbX,nbX));
    arma::mat A_12_lp = arma::kron(arma::eye(nbY,nbY), arma::ones(nbX,1));
    arma::mat A_1_lp = arma::join_rows(A_11_lp, arma::join_rows(A_12_lp,-kron_term));

    arma::mat A_2_lp = arma::join_rows(arma::zeros(1,nbX+nbY),arma::trans(arma::vectorise(C_hat)));

    arma::mat A_lp = arma::join_cols(A_1_lp, A_2_lp);

    arma::vec rhs_lp = arma::join_cols(arma::zeros(nbX*nbY,1),arma::ones(1,1));

    arma::vec lb_lp = arma::zeros(nbX + nbY + dim_theta,1);
    lb_lp.rows(nbX+nbY,nbX+nbY+dim_theta-1).fill(-arma::datum::inf);

    const int k_lp = A_lp.n_rows;
    const int n_lp = A_lp.n_cols;

    char* sense_lp = new char[k_lp];
    for (int jj=0; jj < k_lp - 1; jj++) {
        sense_lp[jj] = '>';
    }
    sense_lp[k_lp-1] = '=';

    const int modelSense = 0; // minimize

    arma::mat sol_mat(n_lp, 2);
    arma::mat dual_mat(k_lp, 2);

    bool LP_optimal = false;
    double val_lp = 0.0;
    //
    try {
        LP_optimal = generic_LP(k_lp, n_lp, obj_lp.memptr(), A_lp.memptr(), modelSense, rhs_lp.memptr(), sense_lp, nullptr, lb_lp.memptr(), nullptr, nullptr, val_lp, sol_mat.colptr(0), sol_mat.colptr(1), dual_mat.colptr(0), dual_mat.colptr(1));

        if (LP_optimal) {
            theta_hat = sol_mat(arma::span(nbX+nbY,nbX+nbY+dim_theta-1),0);
            //
            success = true;
        } else {
            std::cout << "Non-optimal value found during optimization" << std::endl;
        }
    } catch(...) {
        std::cout << "Exception during optimization" << std::endl;
    }
    //
    return success;
}

/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
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
 * general model class
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * This version:
 * 11/27/2016
 */

// first method to build
template<typename Ta>
void model<Ta>::build(const arma::cube& phi_xyk_inp)
{
    this->build_int(phi_xyk_inp,NULL,NULL);
}

template<typename Ta>
void model<Ta>::build(const arma::cube& phi_xyk_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    this->build_int(phi_xyk_inp,&n_inp,&m_inp);
}

template<typename Ta>
void model<Ta>::build_int(const arma::cube& phi_xyk_inp, const arma::vec* n_inp, const arma::vec* m_inp)
{
    need_norm = false;

    nbX = phi_xyk_inp.n_rows;
    nbY = phi_xyk_inp.n_cols;
    nbParams = phi_xyk_inp.n_slices;
    //
    if (n_inp) {
        n = *n_inp;
    } else {
        n.ones(nbX,1);
    }
    if (m_inp) {
        m = *m_inp;
    } else {
        m.ones(nbY,1);
    }

    phi_xyk = phi_xyk_inp;
    //
}

// second method to build
template<typename Ta>
void model<Ta>::build(const arma::mat& X_inp, const arma::mat& Y_inp)
{
    this->build_int(X_inp,Y_inp,NULL,NULL);
}

template<typename Ta>
void model<Ta>::build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp)
{
    this->build_int(X_inp,Y_inp,&n_inp,&m_inp);
}

template<typename Ta>
void model<Ta>::build_int(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec* n_inp, const arma::vec* m_inp)
{
    need_norm = false;

    nbX = X_inp.n_rows;
    nbY = Y_inp.n_rows;

    dX = X_inp.n_cols;
    dY = Y_inp.n_cols;

    nbParams = dX*dY;
    //
    if (n_inp) {
        n = *n_inp;
    } else {
        n.ones(nbX,1);
    }
    if (m_inp) {
        m = *m_inp;
    } else {
        m.ones(nbY,1);
    }

    arma::mat phi_xy_temp = arma::kron(Y_inp,X_inp);
    arma::cube phi_xyk_temp(phi_xy_temp.memptr(),nbX,nbY,nbParams,false); // share memory

    phi_xyk = phi_xyk_temp;

    //phi_xyk.set_size(nbX,nbY,nbParams);
    //phi_xyk.memptr() = phi_xy_temp.memptr();
    //phi_xyk(phi_xy_temp.memptr(),nbX,nbY,nbParams);
    //
}

//
// build markets
template<typename Ta>
void model<Ta>::build_market_TU(const arma::mat& theta)
{
    market_obj.build_TU(n,m,Phi_xy_theta(theta),need_norm);
}

// general arums input
template<typename Ta>
void model<Ta>::build_market_TU(const arma::mat& theta, const Ta& arums_G_inp, const Ta& arums_H_inp)
{
    market_obj.build_TU(n,m,Phi_xy_theta(theta),arums_G_inp,arums_H_inp,need_norm);
}

// empirical version
template<typename Ta> template<typename T>
void model<Ta>::build_market_TU(const arma::mat& theta, T arums_G_inp, T arums_H_inp, int nbDraws, int seed)
{
    market_obj.build_TU(n,m,Phi_xy_theta(theta),arums_G_inp,arums_H_inp,nbDraws,seed,need_norm);
}

template<typename Ta>
bool model<Ta>::mme(const arma::mat& mu_hat, arma::mat& theta_hat)
{
    bool success = false;
    //
    //double xtol_rel = (xtol_rel_inp) ? *xtol_rel_inp : 1E-04;
    //int max_eval = (max_eval_inp) ? *max_eval_inp : 1E05;

    arma::vec theta_0;
    init_param(theta_0);

    arma::mat dtheta_Psi;
    dparam(NULL,dtheta_Psi);

    build_market_TU(theta_0);

    arma::mat kron_term = dtheta_Psi;
    arma::mat C_hat = arma::vectorise(arma::trans(arma::vectorise(mu_hat)) * kron_term);
    //
    // add optimization data
    trame_model_opt_data<Ta> opt_data;
    
    opt_data.market = market_obj;
    opt_data.nbParams = nbParams;
    opt_data.C_hat = C_hat;
    opt_data.kron_term = kron_term;
    //
    arma::vec A_0 = arma::join_cols(arma::vectorise(kron_term * theta_0)/2.0,theta_0);
    int n_pars_opt = A_0.n_elem;

    std::vector<double> sol_vec = arma::conv_to< std::vector<double> >::from(A_0);
    double obj_val = 0;

    //std::vector<double> grad_vec;
    //double val_ret = model_mme_opt_objfn(sol_vec,grad_vec,&opt_data);

    success = model_mme_nlopt(n_pars_opt,sol_vec,obj_val,NULL,NULL,model_mme_opt_objfn,opt_data);
    //
    arma::mat sol_mat = arma::conv_to< arma::mat >::from(sol_vec);

    arma::mat U = arma::reshape(sol_mat.rows(0,nbX*nbY-1),nbX,nbY);
    theta_hat = sol_mat.rows(nbX*nbY,nbX*nbY+nbParams-1);
    double val_ret = obj_val;
    //
    return success;
}

template<typename Ta>
bool model<Ta>::mme_2(const arma::mat& mu_hat, arma::mat& theta_hat)
{
    bool success = false;
    //
    //double xtol_rel = (xtol_rel_inp) ? *xtol_rel_inp : 1E-04;
    //int max_eval = (max_eval_inp) ? *max_eval_inp : 1E05;

    arma::vec theta_0;
    init_param(theta_0);

    arma::mat dtheta_Psi;
    dparam(NULL,dtheta_Psi);

    build_market_TU(theta_0);

    arma::mat kron_term = dtheta_Psi;
    arma::mat C_hat = arma::vectorise(arma::trans(arma::vectorise(mu_hat)) * kron_term);
    //
    // add optimization data
    trame_model_opt_data<Ta> opt_data;
    
    opt_data.market = market_obj;
    opt_data.nbParams = nbParams;
    opt_data.C_hat = C_hat;
    opt_data.kron_term = kron_term;
    //
    arma::vec sol_vec = arma::join_cols(arma::vectorise(kron_term * theta_0)/2.0,theta_0);

    double obj_val = 0;

    success = model_mme_optim(sol_vec,model_mme_opt_objfn_2,&opt_data,&obj_val);
    //
    arma::mat U = arma::reshape(sol_vec.rows(0,nbX*nbY-1),nbX,nbY);
    theta_hat = sol_vec.rows(nbX*nbY,nbX*nbY+nbParams-1);

    double val_ret = obj_val;
    //
    return success;
}

template<typename Ta>
void model<Ta>::dparam(const arma::mat* dparams_inp, arma::mat& dparamsPsi_out)
{
    this->dparam(dparams_inp,dparamsPsi_out,NULL,NULL);
}

template<typename Ta>
void model<Ta>::dparam(const arma::mat* dparams_inp, arma::mat& dparamsPsi_out, arma::mat* dparamsG_out, arma::mat* dparamsH_out)
{
    arma::mat dparams_mat = (dparams_inp) ? *dparams_inp : arma::eye(nbParams,nbParams);

    dparamsPsi_out = Phi_xy() * dparams_mat;
    //
    if (dparamsG_out) {
        *dparamsG_out = arma::zeros(0,dparams_mat.n_cols);
    }
    if (dparamsH_out) {
        *dparamsH_out = arma::zeros(0,dparams_mat.n_cols);
    }
}

// Keith: should probably switch this to be a member variable
template<typename Ta>
arma::mat model<Ta>::Phi_xy()
{
    // mirror R's approach to creating a matrix from an array; take each slice and vectorise that matrix
    arma::mat phi_xy_mat(nbX*nbY,nbParams);
    for(int k = 0; k < nbParams; k++) {
        phi_xy_mat.col(k) = arma::vectorise(phi_xyk.slice(k));
    }

    return phi_xy_mat;
}

template<typename Ta>
arma::mat model<Ta>::Phi_xy_theta(const arma::mat& theta)
{
    arma::mat ret = arma::reshape(Phi_xy() * theta,nbX,nbY);
    return ret;
}

template<typename Ta>
void model<Ta>::init_param(arma::mat& params)
{
    params.zeros(nbParams,1);
}

// optimization-related functions

template<typename Ta>
bool model<Ta>::model_mme_nlopt(int n_pars, std::vector<double>& io_val, double& opt_val, double* lb, double* ub,
                                double (*opt_objfn)(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data),
                                trame_model_opt_data<Ta> opt_data)
{
    bool success = false;

    nlopt::opt opt_trame(nlopt::LD_LBFGS, n_pars);
    nlopt::result result;

    if (lb) {
        opt_trame.set_lower_bounds(*lb);
    }
    if (ub) {
        opt_trame.set_upper_bounds(*ub);
    }

    opt_trame.set_min_objective(*opt_objfn, &opt_data);
    
    opt_trame.set_xtol_rel(1E-04);
    //opt_trame.set_ftol_rel(1e-15);
    opt_trame.set_maxeval(1E05);

    double minf;
    try {
        result = opt_trame.optimize(io_val, minf);
    } catch(...) {
        printf("error in mme optimization (using LBFGS);\n");
        printf("retrying with MMA\n");

        nlopt::opt opt_trame_2(nlopt::LD_MMA, n_pars);

        if (lb) {
            opt_trame_2.set_lower_bounds(*lb);
        }
        if (ub) {
            opt_trame_2.set_upper_bounds(*ub);
        }

        opt_trame_2.set_min_objective(*opt_objfn, &opt_data);
        
        opt_trame_2.set_xtol_rel(1E-04);
        //opt_trame_2.set_ftol_rel(1e-15);
        opt_trame_2.set_maxeval(1E05);

        result = opt_trame_2.optimize(io_val, minf);
    }

    if (result > 0) {
        opt_val = minf;
        success = true;
    }

    return success;
}

template<typename Ta>
double model<Ta>::model_mme_opt_objfn(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data)
{
    //std::cout << "begin opt!" << std::endl;
    trame_model_opt_data<Ta> *d = reinterpret_cast<trame_model_opt_data<Ta>*>(opt_data);
    //
    int nbX = d->market.nbX;
    int nbY = d->market.nbY;
    int nbParams = d->nbParams;
    arma::mat C_hat = d->C_hat;
    arma::mat kron_term = d->kron_term;

    arma::mat x_mat = arma::conv_to< arma::mat >::from(x_inp);
    
    arma::mat U = arma::reshape(x_mat.rows(0,nbX*nbY-1),nbX,nbY);

    arma::mat theta = x_mat.rows(nbX*nbY,nbParams + nbX*nbY - 1);
    arma::mat phi_mat = arma::reshape(kron_term * theta,nbX,nbY);
    //
    arma::mat mu_G, mu_H;

    double val_G = d->market.arums_G.G(d->market.n,U,mu_G);
    double val_H = d->market.arums_H.G(d->market.m,arma::trans(phi_mat - U),mu_H);
    //
    if (!grad.empty()) {
        arma::vec grad_U = arma::vectorise(mu_G - mu_H.t());
        arma::vec grad_theta = arma::vectorise( arma::trans(arma::vectorise(mu_H.t())) * kron_term ) - C_hat;
        arma::vec grad_vec = arma::join_cols(grad_U,grad_theta);

        grad = arma::conv_to< std::vector<double> >::from(grad_vec);
    }
    //
    double ret = val_G + val_H - arma::accu(theta%C_hat);
    //
    return ret;
}

template<typename Ta>
bool model<Ta>::model_mme_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, double* value_out)
{
    bool success = generic_optim(init_out_vals,opt_objfn,opt_data,value_out);
    //
    return success;
}

template<typename Ta>
double model<Ta>::model_mme_opt_objfn_2(const arma::vec& vals_inp, arma::vec* grad, void* opt_data)
{
    trame_model_opt_data<Ta> *d = reinterpret_cast<trame_model_opt_data<Ta>*>(opt_data);
    //
    int nbX = d->market.nbX;
    int nbY = d->market.nbY;
    int nbParams = d->nbParams;
    arma::mat C_hat = d->C_hat;
    arma::mat kron_term = d->kron_term;
    //    
    arma::mat U = arma::reshape(vals_inp.rows(0,nbX*nbY-1),nbX,nbY);

    arma::mat theta = vals_inp.rows(nbX*nbY,nbParams + nbX*nbY - 1);
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

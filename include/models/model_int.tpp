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
 * internal functions
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * This version:
 * 08/20/2017
 */

//
// function declarations

template<typename Tm>
arma::mat model_build_int(const Tm& market_obj, const arma::mat& X_inp, const arma::mat& Y_inp);

template<typename Tm>
void model_to_market_int(Tm& market_obj, const arma::mat& model_data, const arma::mat& theta, const arma::vec& n, const arma::vec& m, int nbX, int nbY, int dX, int dY, bool need_norm);

template<typename Tg, typename Th, typename Tt>
void model_to_market_int(dse<Tg,Th,Tt>& market_obj, const arma::mat& model_data, const arma::mat& theta, const Tg& arums_G_inp, const Th& arums_H_inp, const arma::vec& n, const arma::vec& m, int nbX, int nbY, int dX, int dY, bool need_norm);

template<typename Tt>
void model_to_market_int(mfe<Tt>& market_obj, const arma::mat& model_data, const arma::mat& theta, const arma::vec& n, const arma::vec& m, int nbX, int nbY, int dX, int dY, double sigma, bool need_norm);

template<typename Tm>
void model_dmu(Tm& market_obj, const arma::mat& dtheta_Psi, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& dmu_out);

//
// structs

template<typename Tm>
struct trame_model_mme_opt_data {
    int dim_theta;

    arma::vec C_hat;
    arma::mat kron_term;

    Tm market;
};

template<typename Tm>
struct trame_model_mle_opt_data {
    bool by_individual;
    double scale;

    arma::mat mu_hat;
    arma::vec mu_hat_x0;
    arma::vec mu_hat_0y;

    model<Tm> model_obj;
};

struct trame_model_mfe_mme_opt_data {
    int max_iter_ipfp;
    double tol_ipfp;

    int nbX;
    int nbY;

    int dX;
    int dY;

    double sigma;

    arma::vec p;
    arma::vec q;

    arma::mat IX;
    arma::mat tIY;

    arma::mat f;
    arma::mat g;

    arma::mat v;
    arma::mat Pi_hat;

    arma::mat phi_xy; // should be (nbX*nbY) x (dim_params)
};

//
// functions with specializations
//

//
// build function

template<typename Tg, typename Th>
arma::mat
model_build_int(const dse<Tg,Th,transfers::etu>& market_obj, const arma::mat& X_inp, const arma::mat& Y_inp)
{
    int nbX = X_inp.n_rows;
    int nbY = Y_inp.n_rows;

    arma::mat model_data = arma::abs(arma::kron(arma::ones(nbY,1),X_inp) - arma::kron(Y_inp,arma::ones(nbX,1)));

    return model_data;
}

template<typename Tg, typename Th>
arma::mat
model_build_int(const dse<Tg,Th,transfers::tu>& market_obj, const arma::mat& X_inp, const arma::mat& Y_inp)
{
    int nbX = X_inp.n_rows;
    int nbY = Y_inp.n_rows;

    int dX = X_inp.n_cols;
    int dY = Y_inp.n_cols;

    int dim_theta = dX*dY;
    //
    arma::mat phi_xy_temp = arma::kron(Y_inp,X_inp);
    arma::cube phi_xyk_temp(phi_xy_temp.memptr(),nbX,nbY,dim_theta,false); // share memory

    arma::mat model_data = cube_to_mat(phi_xyk_temp);

    return model_data;
}

template<>
inline
arma::mat
model_build_int(const mfe<mmfs::geo>& market_obj, const arma::mat& X_inp, const arma::mat& Y_inp)
{
    int nbX = X_inp.n_rows;
    int nbY = Y_inp.n_rows;

    int dX = X_inp.n_cols;
    int dY = Y_inp.n_cols;

    int dim_theta = dX*dY;
    //
    arma::mat phi_xy_temp = arma::kron(Y_inp,X_inp);
    arma::cube phi_xyk_temp(phi_xy_temp.memptr(),nbX,nbY,dim_theta,false); // share memory

    arma::mat model_data = cube_to_mat(phi_xyk_temp);

    return model_data;
}

//
// model |-> market

template<typename Tg, typename Th>
void
model_to_market_int(dse<Tg,Th,transfers::etu>& market_obj, const arma::mat& model_data, const arma::mat& theta, const arma::vec& n, const arma::vec& m, int nbX, int nbY, int dX, int dY, bool need_norm)
{
    arma::mat theta_1 = theta.rows(0,dX-1);
    arma::mat theta_2 = theta.rows(dX,dX+dY-1);
    double theta_3 = theta(theta.n_rows-1);

    arma::mat alpha = arma::reshape(model_data*theta_1,nbX,nbY);
    arma::mat gamma = arma::reshape(model_data*theta_2,nbX,nbY);
    arma::mat tau(nbX,nbY);
    tau.fill(theta_3);

    market_obj.build(n,m,alpha,gamma,tau,need_norm);
}

template<typename Tg, typename Th>
void
model_to_market_int(dse<Tg,Th,transfers::etu>& market_obj, const arma::mat& model_data, const arma::mat& theta, const Tg& arums_G_inp, const Th& arums_H_inp, const arma::vec& n, const arma::vec& m, int nbX, int nbY, int dX, int dY, bool need_norm)
{
    arma::mat theta_1 = theta.rows(0,dX-1);
    arma::mat theta_2 = theta.rows(dX,dX+dY-1);
    double theta_3 = theta(theta.n_rows-1);

    arma::mat alpha = arma::reshape(model_data*theta_1,nbX,nbY);
    arma::mat gamma = arma::reshape(model_data*theta_2,nbX,nbY);
    arma::mat tau(nbX,nbY);
    tau.fill(theta_3);

    market_obj.build(n,m,alpha,gamma,tau,arums_G_inp,arums_H_inp,need_norm);
}

template<typename Tg, typename Th>
void
model_to_market_int(dse<Tg,Th,transfers::tu>& market_obj, const arma::mat& model_data, const arma::mat& theta, const arma::vec& n, const arma::vec& m, int nbX, int nbY, int dX, int dY, bool need_norm)
{
    arma::mat phi = arma::reshape(model_data*theta,nbX,nbY);
    market_obj.build(n,m,phi,need_norm);
}

template<typename Tg, typename Th>
void
model_to_market_int(dse<Tg,Th,transfers::tu>& market_obj, const arma::mat& model_data, const arma::mat& theta, const Tg& arums_G_inp, const Th& arums_H_inp, const arma::vec& n, const arma::vec& m, int nbX, int nbY, int dX, int dY, bool need_norm)
{
    arma::mat phi = arma::reshape(model_data*theta,nbX,nbY);
    market_obj.build(n,m,phi,arums_G_inp,arums_H_inp,need_norm);
}

template<>
inline
void
model_to_market_int(mfe<mmfs::geo>& market_obj, const arma::mat& model_data, const arma::mat& theta, const arma::vec& n, const arma::vec& m, int nbX, int nbY, int dX, int dY, double sigma, bool need_norm)
{
    arma::mat phi = arma::reshape(model_data*theta,nbX,nbY);

    market_obj.build(sigma,need_norm);
    market_obj.build(n,m,phi);
}

//
// gradient

template<typename Tg, typename Th, typename Tt>
void
model_dmu(dse<Tg,Th,Tt>& market_obj, const arma::mat& dtheta_Psi, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& dmu_out)
{
    arma::mat mu, U, V;
    char slv = 'j';
    // market_obj.solve(mu,U,V,nullptr);
    market_obj.solve(mu,U,V,&slv);

    arma::vec mu_x0 = market_obj.n - arma::sum(mu,1);
    arma::vec mu_0y = market_obj.m - arma::trans(arma::sum(mu,0));

    arma::vec du_Psi_vec = arma::vectorise(market_obj.trans_obj.du_Psi(U,V));
    arma::vec dv_Psi_vec = 1.0 - du_Psi_vec;
    //
    arma::mat HessGstar = market_obj.arums_G.D2Gstar(market_obj.n,mu,true);
    arma::mat HessHstar = market_obj.arums_H.D2Gstar(market_obj.m,mu.t(),false);
    //
    arma::mat denom = elem_prod(du_Psi_vec,HessGstar) + elem_prod(dv_Psi_vec,HessHstar);
    arma::mat term_1 = market_obj.trans_obj.dparams_Psi(U,V,&dtheta_Psi);

    arma::mat dmu = - arma::solve(denom,term_1);
    //
    mu_out = mu;
    mu_x0_out = mu_x0;
    mu_0y_out = mu_0y;
    dmu_out = dmu;
}

template<>
inline
void
model_dmu(dse<arums::logit,arums::logit,transfers::tu>& market_obj, const arma::mat& dtheta_Psi, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& dmu_out)
{ // dtheta_mu_logit of the R version

    const int nbX = market_obj.nbX;
    const int nbY = market_obj.nbY;
    const int dim_theta = dtheta_Psi.n_cols;

    const double sigma = market_obj.arums_G.sigma;

    arma::vec mu_x0, mu_0y, u_vec, v_vec;
    arma::mat mu, U, V;

    mfe<mmfs::geo> mfe_obj(market_obj.n,market_obj.m);
    mfe_obj.mmfs_obj = market_obj.trans_obj.gen_mmf();

    ipfp(mfe_obj,mu,mu_x0,mu_0y,U,V,u_vec,v_vec,nullptr,nullptr,nullptr);  // should be ipfp

    //

    // arma::mat u_mat = arma::reshape(u_vec,market_obj.nbX,market_obj.nbY); // us
    arma::mat u_mat = arma::repmat(u_vec,1,market_obj.nbY);        // us
    arma::mat v_mat = byrow(v_vec,market_obj.nbX,market_obj.nbY);  // vs

    arma::mat du_Psi_mat = market_obj.trans_obj.du_Psi(U,V);
    arma::mat dv_Psi_mat = 1.0 - du_Psi_mat;

    arma::mat delta_Psi = market_obj.trans_obj.dparams_Psi(u_mat,v_mat,&dtheta_Psi);

    // arma::mat mu_delta_Psi = arma::vectorise(mu) % arma::vectorise(delta_Psi);
    arma::mat mu_delta_Psi = elem_prod(arma::vectorise(mu), delta_Psi);

    arma::cube mu_delta(mu_delta_Psi.begin(), nbX, nbY, dim_theta, false);

    arma::mat d_1 = cube_sum(mu_delta,0) / sigma; // R: apply(mudeltaPsi, c(1,3), sum) / sigma
    arma::mat d_2 = cube_sum(mu_delta,1) / sigma;

    arma::mat numer = arma::join_cols(d_1,d_2);
    //
    const arma::mat Delta_11 = arma::diagmat(mu_x0 + arma::sum(elem_prod(mu,du_Psi_mat),1));
    const arma::mat Delta_22 = arma::diagmat(mu_0y + arma::trans(arma::sum(elem_prod(mu,dv_Psi_mat),0)));

    const arma::mat Delta_12 = elem_prod(mu,dv_Psi_mat);
    const arma::mat Delta_21 = arma::trans(elem_prod(mu,du_Psi_mat));

    const arma::mat Delta = arma::join_cols( arma::join_rows(Delta_11,Delta_12), arma::join_rows(Delta_21,Delta_22) );
    //
    arma::mat dlog_mu_sngls = arma::solve(Delta,numer);

    arma::mat dlog_mu_x0 = dlog_mu_sngls.rows(0,nbX-1);
    arma::mat dlog_mu_0y = dlog_mu_sngls.rows(nbX,nbX+nbY-1);

    // arma::cube dlog_mu_x0_full(nbX,nbY,dim_theta);
    // arma::cube dlog_mu_0y_full(nbX,nbY,dim_theta);

    arma::mat dlog_mu_x0_full_mat = arma::repmat(dlog_mu_x0,nbY,1); // R: matrix(dlogmux0full, ncol=rangeTheta)

    arma::mat dlog_mu_0y_full_mat(nbX*nbY,dim_theta);

    for (int y=0; y < nbY; y++) {
        dlog_mu_0y_full_mat.rows(y*nbX,(y+1)*nbX-1) = arma::repmat(dlog_mu_0y.row(y),nbX,1);
    }

    arma::mat dlogmu = elem_prod(arma::vectorise(du_Psi_mat),dlog_mu_x0_full_mat) + elem_prod(arma::vectorise(dv_Psi_mat),dlog_mu_0y_full_mat) - delta_Psi/sigma;
    //
    mu_out = mu;
    mu_x0_out = mu_x0;
    mu_0y_out = mu_0y;
    dmu_out = elem_prod(arma::vectorise(mu),dlogmu);
}

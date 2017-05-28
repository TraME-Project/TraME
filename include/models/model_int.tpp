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

//
// function declarations

template<typename Tm>
arma::mat model_build_int(const Tm& market_obj, const arma::mat& X_inp, const arma::mat& Y_inp);

template<typename Tm>
void model_dmu(const Tm& market_obj, const arma::mat& dtheta_Psi, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& dmu_out);

//
// structs

template<typename Tm>
struct trame_model_mme_opt_data {
    int dim_theta;

    arma::mat C_hat;
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

//
// functions with specializations

template<typename Tg, typename Th, typename Tt>
void
model_dmu(const dse<Tg,Th,Tt>& market_obj, const arma::mat& dtheta_Psi, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& dmu_out)
{
    arma::mat mu, U, V;
    market_obj.solve(mu,U,V,NULL);

    arma::vec mu_x0 = market_obj.n - arma::sum(mu,1);
    arma::vec mu_0y = market_obj.m - arma::trans(arma::sum(mu,0));

    arma::vec du_Psi_vec = arma::vectorise(market_obj.trans_obj.du_Psi(U,V));
    arma::vec dv_Psi_vec = 1.0 - du_Psi_vec;
    //
    arma::mat HessGstar = market_obj.arums_G.D2Gstar(market_obj.n,mu,true);
    arma::mat HessHstar = market_obj.arums_H.D2Gstar(market_obj.m,mu.t(),false);
    //
    arma::mat denom = elem_prod(du_Psi_vec,HessGstar) + elem_prod(dv_Psi_vec,HessHstar);
    arma::mat term_1 = market_obj.trans_obj.dparams_Psi(U,V,dtheta_Psi);

    arma::mat dmu = - arma::solve(denom,term_1);
    //
    mu_out = mu;
    mu_x0_out = mu_x0;
    mu_0y_out = mu_0y;
    dmu_out = dmu;
}

template<typename Tg, typename Th>
arma::mat 
model_build_int(const dse<Tg,Th,transfers::etu>& market_obj, const arma::mat& X_inp, const arma::mat& Y_inp)
{
    int nbX = X_inp.n_rows;
    int nbY = Y_inp.n_rows;

    // int dX = X_inp.n_cols;
    // int dY = Y_inp.n_cols;

    arma::mat eX = arma::ones(nbX,1);
    arma::mat eY = arma::ones(nbY,1);
    //
    arma::mat model_data = arma::abs(arma::kron(eY,X_inp) - arma::kron(Y_inp,eX));

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

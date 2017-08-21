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
 * 08/20/2017
 */

#ifndef _trame_model_HPP
#define _trame_model_HPP

class model_base
{
    public:
        // build objects
        bool need_norm;

        int nbX;
        int nbY;

        int dX;
        int dY;
        int dim_theta;

        arma::vec n;
        arma::vec m;
};

template<class Tm>
class model : public model_base
{
    Tm market_obj;
};

template<typename Tg, typename Th, typename Tt>
class model<dse<Tg,Th,Tt>> : public model_base
{
    public:
        // build objects
        arma::mat model_data;

        dse<Tg,Th,Tt> market_obj;
        // member functions
        void build(const arma::cube& phi_xyk_inp);
        void build(const arma::cube& phi_xyk_inp, const arma::vec& n_inp, const arma::vec& m_inp);
        void build(const arma::mat& X_inp, const arma::mat& Y_inp);
        void build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp);

        // template<typename Ta, typename Tb> void build_market_TU(const arma::mat& theta, const Ta& arums_G_inp, const Tb& arums_H_inp, int nb_draws, int seed);

        void model_to_market(const arma::mat& theta);
        void model_to_market(const arma::mat& theta, const Tg& arums_G_inp, const Th& arums_H_inp);

        void dtheta(const arma::mat* delta_theta_inp, arma::mat& dtheta_Psi_out);
        void dtheta(const arma::mat* delta_theta_inp, arma::mat& dtheta_Psi_out, arma::mat* dtheta_G_out, arma::mat* dtheta_H_out);

        void dtheta_mu(const arma::mat& theta, const arma::mat* dtheta, arma::mat& mu_out, arma::vec& mu_x0_out, arma::vec& mu_0y_out, arma::mat& dmu);

        bool mme(const arma::mat& mu_hat, arma::mat& theta_hat, const arma::mat* theta_0_inp);
        bool mme(const arma::mat& mu_hat, arma::mat& theta_hat, const arma::mat* theta_0_inp, const int* optim_method_inp);
        bool mme(const arma::mat& mu_hat, arma::mat& theta_hat, double* val_out, arma::mat* mu_out, arma::mat* U_out, arma::mat* V_out);

        bool mle(const arma::mat& mu_hat, arma::mat& theta_hat, const arma::mat* theta_0_inp);
        bool mle(const arma::mat& mu_hat, arma::mat& theta_hat, const arma::mat* theta_0_inp, const int* optim_method_inp);

        // solve wrappers
        bool solve(arma::mat& mu_sol);
        bool solve(arma::mat& mu_sol, const char* solver);
        bool solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver);

    // private:
        // internal build functions
        void build_int(const arma::cube& phi_xyk_inp, const arma::vec* n_inp, const arma::vec* m_inp);
        void build_int(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec* n_inp, const arma::vec* m_inp);

        void initial_theta(arma::mat& theta_0);
        arma::mat initial_theta();

        // optimization-related objects
        bool model_mle_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, optim::opt_settings* settings_inp, const int optim_method);
        static double log_likelihood(const arma::vec& vals_inp, arma::vec* grad, void* opt_data);

        bool model_mme_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, optim::opt_settings* settings_inp, const int optim_method);
        static double model_mme_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void* opt_data);
};

template<typename Tt>
class model<mfe<Tt>> : public model_base
{
    public:
        // build objects
        arma::mat model_data;

        mfe<Tt> market_obj;
        // member functions
        void build(const arma::mat& X_inp, const arma::mat& Y_inp);
        void build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp);
        void build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp, const double sigma_inp);

        void model_to_market(const arma::mat& theta);

        void dtheta(const arma::mat* delta_theta_inp, arma::mat& dtheta_M_out);
        arma::mat dtheta(const arma::mat* delta_theta_inp);

        bool mme_regul(const arma::mat& mu_hat, const double lambda, arma::mat& theta_hat, double& val_ret, double* xtol_rel_inp, int* max_eval_inp, double* tol_ipfp_inp, double* max_iter_ipfp_inp);
        bool mme_woregul(const arma::mat& mu_hat, arma::mat& theta_hat, double& val_ret, double* xtol_ret, int* max_iter, double* tol_ipfp, double* max_iter_ipfp);

        bool mme(const arma::mat& mu_hat, arma::mat& theta_hat);
        bool mme(const arma::mat& mu_hat, double lambda_inp, arma::mat& theta_hat);
        bool mme(const arma::mat& mu_hat, arma::mat& theta_hat, double* val_out, arma::mat* mu_out, arma::mat* U_out, arma::mat* V_out);

        // solve wrappers
        bool solve(arma::mat& mu_sol);
        bool solve(arma::mat& mu_sol, const char* solver);
        bool solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver);

    private:
        // internal build functions
        void build_int(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec* n_inp, const arma::vec* m_inp, const double* sigma_inp);

        void initial_theta(arma::mat& theta_0);
        arma::mat initial_theta();

        arma::mat Phi_k(const arma::mat& mu_hat);

        // optimization-related objects

        bool model_mme_optim(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad, void* opt_data)> opt_objfn, void* opt_data, double* value_out, const double* err_tol_inp, const int* max_iter_inp);
        static double model_mfe_mme_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void* opt_data);
};

#include "model_int.tpp"
#include "model_dse.tpp"
#include "model_mfe.tpp"

#endif

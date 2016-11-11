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
 * Derived classes to provide wrappers to the TraME library
 */

class transfers_R : public trame::transfers
{
    public:
        void build_ETU_R(arma::mat alpha_ETU, arma::mat gamma_ETU, arma::mat tau_ETU);
        void build_LTU_R(arma::mat lambda_LTU, arma::mat phi_LTU);
        void build_NTU_R(arma::mat alpha_NTU, arma::mat gamma_NTU);
        void build_TU_R(arma::mat phi_TU);

        void trans_R();

        SEXP Psi_R(arma::mat U, arma::mat V);
        SEXP Psi_R(arma::mat U, arma::mat V, Rcpp::IntegerVector x_ind, Rcpp::IntegerVector y_ind);

        SEXP du_Psi_R(arma::mat U, arma::mat V);
        SEXP du_Psi_R(arma::mat U, arma::mat V, Rcpp::IntegerVector x_ind, Rcpp::IntegerVector y_ind);

        SEXP Ucal_R(arma::mat vs);
        SEXP Ucal_R(arma::mat vs, Rcpp::IntegerVector x_ind, Rcpp::IntegerVector y_ind);

        SEXP Vcal_R(arma::mat us);
        SEXP Vcal_R(arma::mat us, Rcpp::IntegerVector x_ind, Rcpp::IntegerVector y_ind);
};

class mfe_mmf_R : public trame::mfe<trame::mmf>
{
    public:
        void build_ETU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, arma::mat tau_inp, bool need_norm_inp);
        void build_ETU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, arma::mat tau_inp, double sigma_inp, bool need_norm_inp);

        void build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, bool need_norm_inp);
        void build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, double sigma_inp, bool need_norm_inp);

        void build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp);
        void build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, double sigma_inp, bool need_norm_inp);
        
        void build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp);
        void build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, double sigma_inp, bool need_norm_inp);

        SEXP solve_R();
};

class dse_empirical_R : public trame::dse<trame::empirical>
{
    public:
        void build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, bool need_norm_inp);
        template<typename Ta> void build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, Ta arums_G_inp, Ta arums_H_inp, bool need_norm_inp);
        
        void build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp);
        template<typename Ta> void build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, Ta arums_G_inp, Ta arums_H_inp, bool need_norm_inp);
        
        void build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp);
        template<typename Ta> void build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, Ta arums_G_inp, Ta arums_H_inp, bool need_norm_inp);
        
        SEXP solve_R();
        SEXP solve_R(Rcpp::CharacterVector solver_inp);

        empirical_R get_arums_G();
        void set_arums_G(empirical_R arums_G_inp);
        empirical_R get_arums_H();
        void set_arums_H(empirical_R arums_H_inp);
        void set_arums(empirical_R arums_G_inp, empirical_R arums_H_inp);

        transfers_R get_transfers_R();
        void set_transfers_R(transfers_R trans_obj_inp);
};

class dse_logit_R : public trame::dse<trame::logit>
{
    public:
        void build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, bool need_norm_inp);
        void build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, logit_R arums_G_inp, logit_R arums_H_inp, bool need_norm_inp);
        
        void build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp);
        void build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, logit_R arums_G_inp, logit_R arums_H_inp, bool need_norm_inp);
        
        void build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp);
        void build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, logit_R arums_G_inp, logit_R arums_H_inp, bool need_norm_inp);
        
        SEXP solve_R();
        SEXP solve_R(Rcpp::CharacterVector solver_inp);

        logit_R get_arums_G();
        void set_arums_G(logit_R arums_G_inp);
        logit_R get_arums_H();
        void set_arums_H(logit_R arums_H_inp);
        void set_arums(logit_R arums_G_inp, logit_R arums_H_inp);

        transfers_R get_transfers_R();
        void set_transfers_R(transfers_R trans_obj_inp);
};

class dse_rsc_R : public trame::dse<trame::rsc>
{
    public:
        void build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, bool need_norm_inp);
        void build_LTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat lambda_inp, arma::mat phi_inp, rsc_R arums_G_inp, rsc_R arums_H_inp, bool need_norm_inp);
        
        void build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, bool need_norm_inp);
        void build_NTU_R(arma::vec n_inp, arma::vec m_inp, arma::mat alpha_inp, arma::mat gamma_inp, rsc_R arums_G_inp, rsc_R arums_H_inp, bool need_norm_inp);
        
        void build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, bool need_norm_inp);
        void build_TU_R(arma::vec n_inp, arma::vec m_inp, arma::mat phi_inp, rsc_R arums_G_inp, rsc_R arums_H_inp, bool need_norm_inp);
        
        SEXP solve_R();
        SEXP solve_R(Rcpp::CharacterVector solver_inp);

        rsc_R get_arums_G();
        void set_arums_G(rsc_R arums_G_inp);
        rsc_R get_arums_H();
        void set_arums_H(rsc_R arums_H_inp);
        void set_arums(rsc_R arums_G_inp, rsc_R arums_H_inp);

        transfers_R get_transfers_R();
        void set_transfers_R(transfers_R trans_obj_inp);
};

#include "trame_R_markets.tpp"

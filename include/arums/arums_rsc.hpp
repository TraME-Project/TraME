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
 * Random Scalar Coefficient (RSC) additive random utility model (ARUM) class
 *
 * Keith O'Hara
 * 08/08/2016
 * 
 * This version:
 * 07/03/2017
 */

#ifndef _trame_arums_rsc_HPP
#define _trame_arums_rsc_HPP

class rsc
{
    public:
        // build objects
        int nbX;
        int nbY;
        int dim_params;
        bool outside_option;
        
        double* dist_pars;
        
        arma::mat zeta;
        arma::mat aux_ord;
        
        arma::cube aux_Influence_lhs;
        arma::cube aux_Influence_rhs;
        
        arma::cube aux_DinvPsigma; 
        arma::cube aux_Psigma;
        
        double   (*aux_cdf_eps)(double x, double* dist_pars);
        double (*aux_quant_eps)(double x, double* dist_pars);
       
        double (*aux_pdf_eps)(double x, double* dist_pars);
        double (*aux_pot_eps)(double x, double* dist_pars);
        
        arma::vec   (*aux_cdf_eps_vec)(arma::vec x, double* dist_pars);
        arma::vec (*aux_quant_eps_vec)(arma::vec x, double* dist_pars);
       
        arma::vec (*aux_pdf_eps_vec)(arma::vec x, double* dist_pars);
        arma::vec (*aux_pot_eps_vec)(arma::vec x, double* dist_pars);
        
        // input objects
        arma::mat U;
        arma::mat mu;
        
        // equilibrium objects
        arma::mat U_sol;
        arma::mat mu_sol;
        
        // member functions
        ~rsc(){};
         rsc(){};
        explicit rsc(const int nbX_inp, const int nbY_inp);
        explicit rsc(const arma::mat& zeta_inp, const bool outside_option_inp);
        explicit rsc(const arma::mat& zeta_inp, const double alpha, const double beta);

        void build(const int nbX_inp, const int nbY_inp);
        void build(const arma::mat& zeta_inp, const bool outside_option_inp);

        void build_beta(const arma::mat& zeta_inp, const double alpha, const double beta);
        
        double G(const arma::vec& n);
        double G(const arma::vec& n, const arma::mat& U_inp, arma::mat& mu_out) const;
        double Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, const int x) const;
        
        double Gstar(const arma::vec& n);
        double Gstar(const arma::vec& n, const arma::mat& mu_inp, arma::mat& U_out) const;
        double Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, const int x) const;

        double Gbar(const arma::mat& Ubar, const arma::mat& mubar, const arma::vec& n, arma::mat& U_out, arma::mat& mu_out) const;
        double Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, int x) const;
        
        arma::mat D2Gstar(const arma::vec& n, const bool x_first) const;
        void D2Gstar(arma::mat &ret, const arma::vec& n, const bool x_first) const;
        arma::mat D2Gstar(const arma::vec& n, const arma::mat& mu_inp, const bool x_first) const;
        void D2Gstar(arma::mat &hess, const arma::vec& n, const arma::mat& mu_inp, const bool x_first) const;

        arma::mat dparams_NablaGstar(const arma::vec& n, const arma::mat* dparams_inp, const bool x_first) const;
        void dparams_NablaGstar(arma::mat &ret, const arma::vec& n, const arma::mat* dparams_inp, const bool x_first) const;
        arma::mat dparams_NablaGstar(const arma::vec& n, const arma::mat& mu_inp, const arma::mat* dparams_inp, const bool x_first) const;
        void dparams_NablaGstar(arma::mat &ret, const arma::vec& n, const arma::mat& mu_inp, const arma::mat* dparams_inp, const bool x_first) const;
        
        empirical simul() const;
        empirical simul(const int nbDraws, const int seed) const;
        void simul(empirical& obj_out) const;
        void simul(empirical& obj_out, const int nbDraws, const int seed) const;
        
    private:
        void simul_int(empirical& obj_out, const int* nbDraws, const int* seed) const;

        double cdf (double x) const;
        arma::vec cdf (arma::vec x) const;
        double pdf (double x) const;
        arma::vec pdf (arma::vec x) const;
        double quantile (double x) const;
        arma::vec quantile (arma::vec x) const;
        double pot (double x) const;
        arma::vec pot (arma::vec x) const;

        static double Gbar_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void* opt_data);
        static double Gbar_opt_constr(const arma::vec& vals_inp, arma::vec* grad, void* constr_data);

        static double Gstarx(arma::vec& U_x, const arma::vec& mu_x_inp, const arma::mat& zeta,
                             const arma::mat& aux_DinvPsigma, const arma::mat& aux_Psigma,
                             const arma::mat& aux_Influence_lhs, const arma::mat& aux_Influence_rhs,
                             arma::vec (*pot_eps_vec)(arma::vec pot_inp, double* dist_pars),
                             arma::vec (*quantile_eps_vec)(arma::vec quant_inp, double* dist_pars),
                             double* dist_pars, int nbY, int x);
};

struct trame_rsc_gbar_opt_data {
    int x;
    int nbY;

    arma::vec Ubar_x;
    arma::mat zeta;
    arma::mat aux_DinvPsigma;
    arma::mat aux_Psigma;
    arma::mat aux_Influence_lhs;
    arma::mat aux_Influence_rhs;

    double* dist_pars;
    arma::vec (*pot_eps_vec)(arma::vec x, double* dist_pars);
    arma::vec (*quantile_eps_vec)(arma::vec x, double* dist_pars);
};

#endif

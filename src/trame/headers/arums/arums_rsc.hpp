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
 * RSC class
 *
 * Keith O'Hara
 * 08/08/2016
 */

class rsc
{
    public:
        // build objects
        int nbX;
        int nbY;
        int nbParams;
        bool outsideOption;
        
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
        arma::mat mu;
        arma::mat U;
        
        // equilibrium objects
        arma::mat mu_sol;
        arma::mat U_sol;
        
        // member functions
        ~rsc(){};
         rsc(){};
        explicit rsc(arma::mat zeta_inp, double alpha, double beta);

        void build(arma::mat zeta_inp, bool outsideOption_inp);

        void build_beta(arma::mat zeta_inp, double alpha, double beta);
        
        double G(arma::vec n);
        double G(arma::vec n, const arma::mat& U_inp, arma::mat& mu_out);
        double Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, int x);
        
        double Gstar(arma::vec n);
        double Gstar(arma::vec n, const arma::mat& mu_inp, arma::mat& U_out);
        double Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x);
        static double Gstarx(arma::vec& U_x, arma::vec mu_x_inp, arma::mat zeta, 
                             arma::mat aux_DinvPsigma, arma::mat aux_Psigma, 
                             arma::mat aux_Influence_lhs, arma::mat aux_Influence_rhs,
                             arma::vec (*pot_eps_vec)(arma::vec pot_inp, double* dist_pars),
                             arma::vec (*quantile_eps_vec)(arma::vec quant_inp, double* dist_pars),
                             double* dist_pars, int nbY, int x);

        double Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_out, arma::mat& mu_out);
        double Gbarx(arma::vec Ubarx, arma::vec mubarx, arma::mat& U_x_out, arma::mat& mu_x_out, int x);
        
        void D2Gstar (arma::mat& hess, arma::vec n, bool x_first);
        void dtheta_NablaGstar (arma::mat& ret, arma::vec n, arma::mat* dtheta, bool x_first);
        
        void simul(empirical &ret, int nbDraws, int seed);
        
        double cdf (double x);
        arma::vec cdf (arma::vec x);
        double pdf (double x);
        arma::vec pdf (arma::vec x);
        double quantile (double x);
        arma::vec quantile (arma::vec x);
        double pot (double x);
        arma::vec pot (arma::vec x);
        
    private:
        static double Gbar_opt_objfn(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data);
        static double Gbar_opt_constr(const std::vector<double> &x_inp, std::vector<double> &grad, void *constr_data);
};

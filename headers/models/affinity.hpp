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
 * affinity model class
 *
 * Keith O'Hara
 * 09/20/2016
 *
 * This version:
 * 11/27/2016
 */

class affinity
{
    public:
        // build objects
        bool need_norm;

        int nbX;
        int nbY;

        int dX;
        int dY;
        int nbParams;

        double sigma;

        arma::vec n;
        arma::vec m;

        arma::mat phi_xyk_aux;
        // member functions
        void build(const arma::mat& X_inp, const arma::mat& Y_inp);
        void build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp);
        void build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp, double sigma_inp);

        mfe<mmf> build_market(const arma::mat& theta);

        void dparam(arma::mat* dparams_inp, arma::mat& dparamsPsi_out, arma::mat* dparamsG_out, arma::mat* dparamsH_out);

        bool mme_regul(const arma::mat& mu_hat, const double& lambda, arma::mat& theta_hat, double& val_ret, double* xtol_rel_inp, int* max_eval_inp, double* tol_ipfp_inp, double* max_iter_ipfp_inp);
        bool mme_woregul(const arma::mat& mu_hat, arma::mat& theta_hat, double& val_ret, double* xtol_ret, int* max_iter, double* tol_ipfp, double* max_iter_ipfp);

        bool mme(const arma::mat& mu_hat);
        bool mme(const arma::mat& mu_hat, const arma::mat& lambda);
        bool mme(const arma::mat& mu_hat, const arma::mat* lambda, double* xtol_rel_inp, int* max_eval_inp, double* tol_ipfp_inp, double* max_iter_ipfp_inp);
    
    private:
        void build_int(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec* n_inp, const arma::vec* m_inp, const double* sigma_inp);

        arma::mat Phi_xyk();
        arma::mat Phi_xy(const arma::mat& lambda);
        arma::mat Phi_k(const arma::mat& mu_hat);

        void init_param(arma::mat& params);
        static double mme_woregul_opt_objfn(const arma::vec& vals_inp, arma::vec* grad, void* opt_data);
};

struct trame_mme_woregal_opt_data {
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

    arma::mat phi_xy; // should be (nbX*nbY) x (nbParams)
};

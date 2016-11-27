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
 * general model class
 *
 * Keith O'Hara
 * 11/19/2016
 *
 * This version:
 * 11/26/2016
 */

template<class Ta>
class model
{
    public:
        // build objects
        bool need_norm;

        int nbX;
        int nbY;

        int dX;
        int dY;
        int nbParams;

        arma::vec n;
        arma::vec m;

        arma::mat phi_xyk;

        dse<Ta> market_obj;
        // member functions
        void build(const arma::mat& X_inp, const arma::mat& Y_inp);
        void build(const arma::mat& X_inp, const arma::mat& Y_inp, const arma::vec& n_inp, const arma::vec& m_inp);

        void build_market_TU(const arma::mat& theta);
        void build_market_TU(const arma::mat& theta, const Ta& arums_G_inp, const Ta& arums_H_inp);
        template<typename T> void build_market_TU(const arma::mat& theta, T arums_G_inp, T arums_H_inp, int nbDraws, int seed);

        void dparam(arma::mat* dparams_inp, arma::mat& dparamsPsi_out, arma::mat* dparamsG_out, arma::mat* dparamsH_out);
        
        bool mme(const arma::mat& mu_hat);
        
    private:
        void init_param(arma::mat& params);
        static double mme_woregul_opt_objfn(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data);
};

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
 * logit additive random utility model (ARUM) class
 *
 * Keith O'Hara
 * 08/08/2016
 *
 * This version:
 * 07/03/2017
 */

#ifndef _trame_arums_logit_HPP
#define _trame_arums_logit_HPP

class logit
{
    public:
        // build objects
        int nbX;
        int nbY;
        int dim_params;

        bool outside_option = true;

        double sigma = 1.0; // setting a default value obviates issues with general market construction (DSE)
        
        // input objects
        arma::mat U;
        arma::mat mu;
        
        // equilibrium objects
        arma::mat U_sol;
        arma::mat mu_sol;
        
        // member functions
        ~logit(){};
         logit(){};
        explicit logit(int nbX_inp, int nbY_inp);
        explicit logit(int nbX_inp, int nbY_inp, double sigma_inp, bool outside_option_inp);

        void build(int nbX_inp, int nbY_inp);
        void build(int nbX_inp, int nbY_inp, double sigma_inp, bool outside_option_inp);
        
        double G(const arma::vec& n);
        double G(const arma::vec& n, const arma::mat& U_inp, arma::mat& mu_out) const;
        //double Gx(const arma::mat& U_x_inp, arma::mat& mu_x_out, int x);
        
        double Gstar(const arma::vec& n);
        double Gstar(const arma::vec& n, const arma::mat& mu_inp, arma::mat& U_out) const;
        double Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out) const;
        double Gstarx(const arma::mat& mu_x_inp, arma::mat &U_x_out, int x) const;
        
        double Gbar(const arma::mat& Ubar, const arma::mat& mubar, const arma::vec& n, arma::mat& U_out, arma::mat& mu_out) const;
        double Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out) const;
        double Gbarx(const arma::vec& Ubar_x, const arma::vec& mubar_x, arma::mat& U_x_out, arma::mat& mu_x_out, int x) const;
        
        arma::mat D2G(const arma::vec& n, bool x_first) const;
        void D2G(arma::mat &H, const arma::vec& n, bool x_first) const;
        arma::mat D2G(const arma::vec& n, const arma::mat& U_inp, bool x_first) const;
        void D2G(arma::mat &H, const arma::vec& n, const arma::mat& U_inp, bool x_first) const;

        arma::mat D2Gstar(const arma::vec& n, bool x_first) const;
        void D2Gstar(arma::mat &H, const arma::vec& n, bool x_first) const;
        arma::mat D2Gstar(const arma::vec& n, const arma::mat& mu_inp, bool x_first) const;
        void D2Gstar(arma::mat &H, const arma::vec& n, const arma::mat& mu_inp, bool x_first) const;

        arma::mat dparams_NablaGstar(const arma::vec& n, arma::mat* dparams_inp, bool x_first) const;
        void dparams_NablaGstar(arma::mat &ret, const arma::vec& n, arma::mat* dparams_inp, bool x_first) const;
        arma::mat dparams_NablaGstar(const arma::vec& n, const arma::mat& mu_inp, arma::mat* dparams_inp, bool x_first) const;
        void dparams_NablaGstar(arma::mat &ret, const arma::vec& n, const arma::mat& mu_inp, arma::mat* dparams_inp, bool x_first) const;
        
        empirical simul() const;
        empirical simul(int nbDraws, int seed) const;
        void simul(empirical& obj_out) const;
        void simul(empirical& obj_out, int nbDraws, int seed) const;

    private:
        static double differMargX(double z, void* opt_data);

        void simul_int(empirical& obj_out, int* nbDraws, int* seed) const;
};

struct trame_logit_zeroin_data {
    arma::mat exp_Ubar_X;
    arma::mat mubar_X;
};

#endif

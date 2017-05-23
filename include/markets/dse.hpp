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
 * Demand-Supply Equilibrium (DSE) market
 *
 * Keith O'Hara
 * 08/17/2016
 *
 * This version:
 * 04/04/2017
 */

template<class Tg, class Th, class Tt>
class dse
{
    public:
        // build objects
        bool ETU = false;
        bool LTU = false;
        bool NTU = false;
        bool TU  = false;

        bool need_norm;
        bool outsideOption;

        int nbX;
        int nbY;

        arma::vec n;
        arma::vec m;

        Tg arums_G;
        Th arums_H;

        Tt trans_obj;

        // member functions
        void build(const arma::vec& n_inp, const arma::vec& m_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, bool need_norm_inp);

        void trans(dse<Th,Tg,Tt>& trans_market_obj) const;
        dse<Th,Tg,Tt> trans() const;

        bool solve(arma::mat& mu_sol);
        bool solve(arma::mat& mu_sol, const char* solver);
        bool solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver);
};

template<class Tg, class Th>
class dse<Tg,Th,transfers::etu>
{
    public:
        // build objects
        bool ETU = true;
        bool LTU = false;
        bool NTU = false;
        bool TU  = false;

        bool need_norm;
        bool outsideOption;

        int nbX;
        int nbY;

        arma::vec n;
        arma::vec m;

        Tg arums_G;
        Th arums_H;

        transfers::etu trans_obj;

        // member functions
        void build(const arma::vec& n_inp, const arma::vec& m_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, bool need_norm_inp);

        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, bool need_norm_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, const Tg& arums_G_inp, const Th& arums_H_inp, bool need_norm_inp);
        template<typename Ta, typename Tb> void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, Ta arums_G_inp, Tb arums_H_inp, int nbDraws, int seed, bool need_norm_inp);

        void trans(dse<Th,Tg,transfers::etu>& trans_market_obj);
        dse<Th,Tg,transfers::etu> trans();

        bool solve(arma::mat& mu_sol);
        bool solve(arma::mat& mu_sol, const char* solver);
        bool solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver);
};

template<class Tg, class Th>
class dse<Tg,Th,transfers::ltu>
{
    public:
        // build objects
        bool ETU = false;
        bool LTU = true;
        bool NTU = false;
        bool TU  = false;

        bool need_norm;
        bool outsideOption;

        int nbX;
        int nbY;

        arma::vec n;
        arma::vec m;

        Tg arums_G;
        Th arums_H;

        transfers::ltu trans_obj;

        // member functions
        void build(const arma::vec& n_inp, const arma::vec& m_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, bool need_norm_inp);

        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda, const arma::mat& phi, bool need_norm_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, const Tg& arums_G_inp, const Th& arums_H_inp, bool need_norm_inp);
        template<typename Ta, typename Tb> void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, Ta arums_G_inp, Tb arums_H_inp, int nbDraws, int seed, bool need_norm_inp);

        void trans(dse<Th,Tg,transfers::ltu>& trans_market_obj);
        dse<Th,Tg,transfers::ltu> trans();

        bool solve(arma::mat& mu_sol);
        bool solve(arma::mat& mu_sol, const char* solver);
        bool solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver);
};

template<class Tg, class Th>
class dse<Tg,Th,transfers::ntu>
{
    public:
        // build objects
        bool ETU = false;
        bool LTU = false;
        bool NTU = true;
        bool TU  = false;

        bool need_norm;
        bool outsideOption;

        int nbX;
        int nbY;

        arma::vec n;
        arma::vec m;

        Tg arums_G;
        Th arums_H;

        transfers::ntu trans_obj;

        // member functions
        void build(const arma::vec& n_inp, const arma::vec& m_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, bool need_norm_inp);

        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha, const arma::mat& gamma, bool need_norm_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const Tg& arums_G_inp, const Th& arums_H_inp, bool need_norm_inp);
        template<typename Ta, typename Tb> void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, Ta arums_G_inp, Tb arums_H_inp, int nbDraws, int seed, bool need_norm_inp);

        void trans(dse<Th,Tg,transfers::ntu>& trans_market_obj);
        dse<Th,Tg,transfers::ntu> trans();

        bool solve(arma::mat& mu_sol);
        bool solve(arma::mat& mu_sol, const char* solver);
        bool solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver);
};

template<class Tg, class Th>
class dse<Tg,Th,transfers::tu>
{
    public:
        // build objects
        bool ETU = false;
        bool LTU = false;
        bool NTU = false;
        bool TU  = true;

        bool need_norm;
        bool outsideOption;

        int nbX;
        int nbY;

        arma::vec n;
        arma::vec m;

        Tg arums_G;
        Th arums_H;

        transfers::tu trans_obj;

        // member functions
        void build(const arma::vec& n_inp, const arma::vec& m_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, bool need_norm_inp);

        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, bool need_norm_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, const Tg& arums_G_inp, const Th& arums_H_inp, bool need_norm_inp);
        template<typename Ta, typename Tb> void build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, Ta arums_G_inp, Tb arums_H_inp, int nbDraws, int seed, bool need_norm_inp);

        void trans(dse<Th,Tg,transfers::tu>& trans_market_obj);
        dse<Th,Tg,transfers::tu> trans();

        bool solve(arma::mat& mu_sol);
        bool solve(arma::mat& mu_sol, const char* solver);
        bool solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver);
};

#include "dse.tpp"

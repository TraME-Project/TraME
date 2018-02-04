/*################################################################################
  ##
  ##   Copyright (C) 2015-2018 the TraME Team:
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
 * 02/04/2018
 */

#ifndef _trame_dse_market_HPP
#define _trame_dse_market_HPP

class dse_base
{
    public:
        // build objects
        bool ETU = false;
        bool LTU = false;
        bool NTU = false;
        bool TU  = false;

        bool need_norm;
        bool outside_option;

        uint_t nbX;
        uint_t nbY;

        arma::vec n;
        arma::vec m;

        // member functions
        void build(const arma::vec& n_inp, const arma::vec& m_inp);
        void build(const arma::vec& n_inp, const arma::vec& m_inp, const bool need_norm_inp);
};

template<class Tg, class Th, class Tt>
class dse : public dse_base
{
    public:
        // build objects
        Tg arums_G;
        Th arums_H;

        Tt trans_obj;

        // member functions

        template<class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::etu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, const bool need_norm_inp);

        template<class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::etu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, const Tg& arums_G_inp, const Th& arums_H_inp, const bool need_norm_inp);

        template<typename Ta, typename Tb, class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::etu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, Ta arums_G_inp, Tb arums_H_inp, const uint_t n_draws, const uint_t seed, const bool need_norm_inp);

        //

        template<class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::ltu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda, const arma::mat& phi, const bool need_norm_inp);

        template<class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::ltu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, const Tg& arums_G_inp, const Th& arums_H_inp, const bool need_norm_inp);

        template<typename Ta, typename Tb, class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::ltu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, Ta arums_G_inp, Tb arums_H_inp, const uint_t n_draws, const uint_t seed, const bool need_norm_inp);

        //

        template<class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::ntu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha, const arma::mat& gamma, const bool need_norm_inp);

        template<class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::ntu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const Tg& arums_G_inp, const Th& arums_H_inp, const bool need_norm_inp);

        template<typename Ta, typename Tb, class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::ntu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, Ta arums_G_inp, Tb arums_H_inp, const uint_t n_draws, const uint_t seed, const bool need_norm_inp);

        //

        template<class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::tu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, const bool need_norm_inp);

        template<class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::tu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, const Tg& arums_G_inp, const Th& arums_H_inp, const bool need_norm_inp);

        template<typename Ta, typename Tb, class Tq = Tt>
        typename std::enable_if<std::is_same<Tq,transfers::tu>::value>::type
        build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, Ta arums_G_inp, Tb arums_H_inp, const uint_t n_draws, const uint_t seed, const bool need_norm_inp);

        //

        void trans(dse<Th,Tg,Tt>& trans_market_obj) const;
        dse<Th,Tg,Tt> trans() const;

        bool solve(arma::mat& mu_sol);
        bool solve(arma::mat& mu_sol, const char* solver);
        bool solve(arma::mat& mu_sol, arma::mat& U, arma::mat& V, const char* solver);

        //
    private:
        void build_simple(const arma::vec& n_inp, const arma::vec& m_inp, const bool need_norm_inp);

        void build_with_arums(const arma::vec& n_inp, const arma::vec& m_inp, const bool need_norm_inp);
        void build_with_arums(const arma::vec& n_inp, const arma::vec& m_inp, const Tg& arums_G_inp, const Th& arums_H_inp, const bool need_norm_inp);

        template<typename Ta, typename Tb>
        void build_with_arums(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, 
                              Ta arums_G_inp, Tb arums_H_inp, const uint_t n_draws, const uint_t seed, const bool need_norm_inp);
};

template<class Tg, class Th, class Tt>
void trans_market(const dse<Tg,Th,Tt>& market_obj, dse<Th,Tg,Tt>& trans_market_obj);

// template<class Tg, class Th, class Tt>
// dse<Th,Tg,Tt> trans_market(const dse<Tg,Th,Tt>& market_obj);

#include "dse.tpp"

#endif

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

//
// build functions
//

//
// short build functions

inline
void
dse_base::build(const arma::vec& n_inp, const arma::vec& m_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;

    need_norm = false;

    n = n_inp;
    m = m_inp;

    outside_option = true;
}

inline
void
dse_base::build(const arma::vec& n_inp, const arma::vec& m_inp, const bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;

    need_norm = need_norm_inp;

    n = n_inp;
    m = m_inp;

    outside_option = (need_norm_inp) ? false : true;
}

//

template<typename Tg, typename Th, typename Tt>
void
dse<Tg,Th,Tt>::build_simple(const arma::vec& n_inp, const arma::vec& m_inp, const bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;

    need_norm = need_norm_inp;

    n = n_inp;
    m = m_inp;

    outside_option = (need_norm_inp) ? false : true;
}

//

template<typename Tg, typename Th, typename Tt>
void
dse<Tg,Th,Tt>::build_with_arums(const arma::vec& n_inp, const arma::vec& m_inp, const bool need_norm_inp)
{
    build_simple(n_inp,m_inp,need_norm_inp);
    //
    arums_G.build(nbX,nbY);
    arums_H.build(nbY,nbX);
}

template<typename Tg, typename Th, typename Tt>
void
dse<Tg,Th,Tt>::build_with_arums(const arma::vec& n_inp, const arma::vec& m_inp, const Tg& arums_G_inp, const Th& arums_H_inp, const bool need_norm_inp)
{
    build_simple(n_inp,m_inp,need_norm_inp);
    //
    arums_G = arums_G_inp;
    arums_H = arums_H_inp;
}

template<typename Tg, typename Th, typename Tt>
template<typename Ta, typename Tb>
void
dse<Tg,Th,Tt>::build_with_arums(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, 
                                Ta arums_G_inp, Tb arums_H_inp, const uint_t n_draws, const uint_t seed, const bool need_norm_inp)
{
    build_simple(n_inp,m_inp,need_norm_inp);
    //
    arums_G_inp.simul(arums_G,n_draws,seed);
    arums_H_inp.simul(arums_H,n_draws,seed);
}

//
// etu case

template<typename Tg, typename Th, typename Tt>
template<typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::etu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,need_norm_inp);
    //
    trans_obj.build(alpha_inp,gamma_inp,tau_inp,need_norm_inp);
    //
    ETU = true;
}

template<typename Tg, typename Th, typename Tt>
template<typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::etu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, const Tg& arums_G_inp, const Th& arums_H_inp, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,arums_G_inp,arums_H_inp,need_norm_inp);
    //
    trans_obj.build(alpha_inp,gamma_inp,tau_inp,need_norm_inp);
    //
    ETU = true;
}

template<typename Tg, typename Th, typename Tt>
template<typename Ta, typename Tb, typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::etu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp, Ta arums_G_inp, Tb arums_H_inp, const uint_t n_draws, const uint_t seed, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,arums_G_inp,arums_H_inp,n_draws,seed,need_norm_inp);
    //
    trans_obj.build(alpha_inp,gamma_inp,tau_inp,need_norm_inp);
    //
    ETU = true;
}

//
// ltu case

template<typename Tg, typename Th, typename Tt>
template<typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::ltu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,need_norm_inp);
    //
    trans_obj.build(lambda_inp,phi_inp,need_norm_inp);
    //
    LTU = true;
}

template<typename Tg, typename Th, typename Tt>
template<typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::ltu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, const Tg& arums_G_inp, const Th& arums_H_inp, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,arums_G_inp,arums_H_inp,need_norm_inp);
    //
    trans_obj.build(lambda_inp,phi_inp,need_norm_inp);
    //
    LTU = true;
}

template<typename Tg, typename Th, typename Tt>
template<typename Ta, typename Tb, typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::ltu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, Ta arums_G_inp, Tb arums_H_inp, const uint_t n_draws, const uint_t seed, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,arums_G_inp,arums_H_inp,n_draws,seed,need_norm_inp);
    //
    trans_obj.build(lambda_inp,phi_inp,need_norm_inp);
    //
    LTU = true;
}

//
// ntu case

template<typename Tg, typename Th, typename Tt>
template<typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::ntu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,need_norm_inp);
    //
    trans_obj.build(alpha_inp,gamma_inp,need_norm_inp);
    //
    NTU = true;
}

template<typename Tg, typename Th, typename Tt>
template<typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::ntu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const Tg& arums_G_inp, const Th& arums_H_inp, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,arums_G_inp,arums_H_inp,need_norm_inp);
    //
    trans_obj.build(alpha_inp,gamma_inp,need_norm_inp);
    //
    NTU = true;
}

template<typename Tg, typename Th, typename Tt>
template<typename Ta, typename Tb, typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::ntu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, Ta arums_G_inp, Tb arums_H_inp, const uint_t n_draws, const uint_t seed, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,arums_G_inp,arums_H_inp,n_draws,seed,need_norm_inp);
    //
    trans_obj.build(alpha_inp,gamma_inp,need_norm_inp);
    //
    NTU = true;
}

//
// tu case

template<typename Tg, typename Th, typename Tt>
template<typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::tu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,need_norm_inp);
    //
    trans_obj.build(phi_inp,need_norm_inp);
    //
    TU = true;
}

template<typename Tg, typename Th, typename Tt>
template<typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::tu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, const Tg& arums_G_inp, const Th& arums_H_inp, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,arums_G_inp,arums_H_inp,need_norm_inp);
    //
    trans_obj.build(phi_inp,need_norm_inp);
    //
    TU = true;
}

template<typename Tg, typename Th, typename Tt>
template<typename Ta, typename Tb, typename Tq>
typename std::enable_if<std::is_same<Tq,transfers::tu>::value>::type
dse<Tg,Th,Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, Ta arums_G_inp, Tb arums_H_inp, const uint_t n_draws, const uint_t seed, const bool need_norm_inp)
{
    build_with_arums(n_inp,m_inp,arums_G_inp,arums_H_inp,n_draws,seed,need_norm_inp);
    //
    trans_obj.build(phi_inp,need_norm_inp);
    //
    TU = true;
}

//
// market transpose

template<typename Tg, typename Th, typename Tt>
void
trans_market(const dse<Tg,Th,Tt>& market_obj, dse<Th,Tg,Tt>& trans_market_obj)
{
    trans_market_obj.nbX = market_obj.nbY;
    trans_market_obj.nbY = market_obj.nbX;

    trans_market_obj.n = market_obj.m;
    trans_market_obj.m = market_obj.n;

    trans_market_obj.trans_obj = market_obj.trans_obj;
    trans_market_obj.trans_obj.trans();

    trans_market_obj.arums_G = market_obj.arums_H;
    trans_market_obj.arums_H = market_obj.arums_G;
}

template<typename Tg, typename Th, typename Tt>
void
dse<Tg,Th,Tt>::trans(dse<Th,Tg,Tt>& trans_market_obj)
const
{
    trans_market(*this,trans_market_obj);
}

template<typename Tg, typename Th, typename Tt>
dse<Th,Tg,Tt> 
dse<Tg,Th,Tt>::trans()
const
{
    dse<Th,Tg,Tt> new_market;

    trans_market(*this,new_market);
    //
    return new_market;
}

//
// solve functions

template<typename Tg, typename Th, typename Tt>
bool
dse<Tg,Th,Tt>::solve(arma::mat& mu_sol)
{
    return equil_solve(*this,mu_sol);
}

template<typename Tg, typename Th, typename Tt>
bool
dse<Tg,Th,Tt>::solve(arma::mat& mu_sol, const char* solver)
{
    return equil_solve(*this,mu_sol,solver);
}

template<typename Tg, typename Th, typename Tt>
bool
dse<Tg,Th,Tt>::solve(arma::mat& mu_sol, arma::mat& U_out, arma::mat& V_out, const char* solver)
{
    return equil_solve(*this,mu_sol,U_out,V_out,solver);
}

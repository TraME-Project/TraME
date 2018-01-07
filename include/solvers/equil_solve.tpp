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
 * generic solve function for DSE markets
 *
 * Keith O'Hara
 * 08/01/2016
 *
 * This version:
 * 07/26/2017
 */

template<typename Tg, typename Th, typename Tt>
bool
equil_solve(const dse<Tg,Th,Tt>& market_obj, arma::mat& mu_sol)
{
    return equil_solve(market_obj,mu_sol,nullptr);
}

template<typename Tg, typename Th, typename Tt>
bool
equil_solve(const dse<Tg,Th,Tt>& market_obj, arma::mat& mu_sol, const char* solver)
{
    bool res = false;
    const char sig = (solver != nullptr) ? solver[0] : char(0);

    if (solver) { // not nullptr
        if (sig=='a') {
            res = arc_newton(market_obj,mu_sol);
        }
        // if (sig=='c') { // only works with empirical case
        //     res = cupids_lp(market_obj,mu_sol);
        // }
        if (sig=='d') {
            res = darum(market_obj,mu_sol);
        }
        if (sig=='e') {
            res = eap_nash(market_obj,mu_sol);
        }
        if (sig=='j') {
            res = jacobi(market_obj,mu_sol);
        }
        if (sig=='m') {
            res = max_welfare(market_obj,mu_sol);
        }
        if (sig=='o') {
            res = oap_lp(market_obj,mu_sol);
        }
    } else { // default
        if (market_obj.NTU) {
            res = darum(market_obj,mu_sol);
        } else if (market_obj.TU) {
            res = max_welfare(market_obj,mu_sol);
        } else {
            res = jacobi(market_obj,mu_sol);
        }
    }
    //
    return res;
}

template<typename Tg, typename Th, typename Tt>
bool
equil_solve(const dse<Tg,Th,Tt>& market_obj, arma::mat& mu_sol, arma::mat& U_out, arma::mat& V_out, const char* solver)
{
    bool res = false;
    const char sig = (solver != nullptr) ? solver[0] : char(0);

    if (solver) { // not nullptr
        if (sig=='a') {
            res = arc_newton(market_obj,mu_sol,U_out,V_out);
        }
        // if (sig=='c') { // only works with empirical case
        //     res = cupids_lp(market_obj,mu_sol);
        // }
        if (sig=='d') {
            res = darum(market_obj,mu_sol,U_out,V_out);
        }
        if (sig=='e') {
            res = eap_nash(market_obj,mu_sol,U_out,V_out);
        }
        if (sig=='j') {
            res = jacobi(market_obj,mu_sol,U_out,V_out);
        }
        if (sig=='m') {
            res = max_welfare(market_obj,mu_sol,U_out,V_out);
        }
        if (sig=='o') {
            // res = oap_lp(market_obj,mu_sol,U_out,V_out);
        }
    } else { // default
        if (market_obj.NTU) {
            res = darum(market_obj,mu_sol,U_out,V_out);
        } else if (market_obj.TU) {
            res = max_welfare(market_obj,mu_sol,U_out,V_out);
        } else {
            res = jacobi(market_obj,mu_sol,U_out,V_out);
        }
    }
    //
    return res;
}

//
// specializations

template<typename Tg, typename Th>
bool
equil_solve(const dse<Tg,Th,transfers::etu>& market_obj, arma::mat& mu_sol, const char* solver)
{
    bool res = false;
    const char sig = (solver != nullptr) ? solver[0] : char(0);

    if (solver) { // not nullptr
        if (sig=='j') {
            res = jacobi(market_obj,mu_sol);
        }
    } else { // default
        res = jacobi(market_obj,mu_sol);
    }
    //
    return res;
}

template<typename Tg, typename Th>
bool
equil_solve(const dse<Tg,Th,transfers::ltu>& market_obj, arma::mat& mu_sol, const char* solver)
{
    bool res = false;
    const char sig = (solver != nullptr) ? solver[0] : char(0);

    if (solver) { // not nullptr
        if (sig=='j') {
            res = jacobi(market_obj,mu_sol);
        }
    } else { // default
        res = jacobi(market_obj,mu_sol);
    }
    //
    return res;
}

template<typename Tg, typename Th>
bool
equil_solve(const dse<Tg,Th,transfers::ntu>& market_obj, arma::mat& mu_sol, const char* solver)
{
    bool res = false;
    const char sig = (solver != nullptr) ? solver[0] : char(0);

    if (solver) { // not nullptr
        if (sig=='d') {
            res = darum(market_obj,mu_sol);
        }
        if (sig=='j') {
            res = jacobi(market_obj,mu_sol);
        }
    } else { // default
        res = darum(market_obj,mu_sol);
    }
    //
    return res;
}

template<>
inline
bool
equil_solve(const dse<arums::logit,arums::logit,transfers::tu>& market_obj, arma::mat& mu_sol, const char* solver)
{
    mfe<mmfs::geo> mfe_obj(market_obj.n,market_obj.m);
    mfe_obj.mmfs_obj = market_obj.trans_obj.gen_mmf();

    return ipfp(mfe_obj,mu_sol);
}

template<>
inline
bool
equil_solve(const dse<arums::logit,arums::logit,transfers::tu>& market_obj, arma::mat& mu_sol, arma::mat& U_out, arma::mat& V_out, const char* solver)
{
    std::cout << "yep" << std::endl;
    mfe<mmfs::geo> mfe_obj(market_obj.n,market_obj.m);
    mfe_obj.mmfs_obj = market_obj.trans_obj.gen_mmf();

    return ipfp(mfe_obj,mu_sol,U_out,V_out);
}

template<typename Tg, typename Th>
bool
equil_solve(const dse<Tg,Th,transfers::tu>& market_obj, arma::mat& mu_sol, const char* solver)
{
    bool res = false;
    const char sig = (solver != nullptr) ? solver[0] : char(0);

    if (solver) { // not nullptr
        if (sig=='a') {
            res = arc_newton(market_obj,mu_sol);
        }
        if (sig=='j') {
            res = jacobi(market_obj,mu_sol);
        }
        if (sig=='m') {
            res = max_welfare(market_obj,mu_sol);
        }
    } else { // default
        res = max_welfare(market_obj,mu_sol);
    }
    //
    return res;
}

// template<typename Tg, typename Th>
// bool
// equil_solve(const dse<Tg,Th,transfers::tu>& market_obj, arma::mat& mu_sol, arma::mat& U_out, arma::mat& V_out, const char* solver)
// {
//     bool res = false;
//     const char sig = (solver != nullptr) ? solver[0] : char(0);

//     if (solver) { // not nullptr
//         if (sig=='a') {
//             std::cout << "arc_newton" << std::endl;
//             res = arc_newton(market_obj,mu_sol,U_out,V_out);
//         }
//         if (sig=='j') {
//             std::cout << "jacobi" << std::endl;
//             res = jacobi(market_obj,mu_sol,U_out,V_out);
//         }
//         if (sig=='m') {
//             std::cout << "max_welfare" << std::endl;
//             res = max_welfare(market_obj,mu_sol,U_out,V_out);
//         }
//     } else { // default
//         res = max_welfare(market_obj,mu_sol,U_out,V_out);
//     }
//     //
//     return res;
// }

// we specialize because cupids_lp is only defined for empirical classes
template<>
inline
bool
equil_solve(const dse<arums::empirical,arums::empirical,transfers::tu>& market_obj, arma::mat& mu_sol, const char* solver)
{
    return cupids_lp(market_obj,mu_sol);
}

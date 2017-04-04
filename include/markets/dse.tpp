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
 * 03/22/2017
 */

template<typename Tg, typename Th, typename Tm>
void 
dse<Tg,Th,Tm>::build_TU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    outsideOption = (need_norm_inp) ? false : true;

    trans_obj.build(phi_inp,need_norm_inp);

    arums_G.build(nbX,nbY); // this avoids nbX and nbY not being set in arums
    arums_H.build(nbY,nbX);
    //
    TU = true;
}

// general arums input
template<typename Tg, typename Th, typename Tm>
void 
dse<Tg,Th,Tm>::build_TU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, const Tg& arums_G_inp, const Th& arums_H_inp, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    outsideOption = (need_norm_inp) ? false : true;

    trans_obj.build(phi_inp,need_norm_inp);

    arums_G = arums_G_inp;
    arums_H = arums_H_inp;
    //
    TU = true;
}

// empirical version
template<typename Tg, typename Th, typename Tm> template<typename Ta, typename Tb>
void 
dse<Tg,Th,Tm>::build_TU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp, Ta arums_G_inp, Tb arums_H_inp, int nbDraws, int seed, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    outsideOption = (need_norm_inp) ? false : true;

    trans_obj.build(phi_inp,need_norm_inp);

    arums_G_inp.simul(arums_G,nbDraws,seed);
    arums_H_inp.simul(arums_H,nbDraws,seed);
    //
    TU = true;
}

// arums none
template<typename Tg, typename Th, typename Tm>
void 
dse<Tg,Th,Tm>::build_NTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    outsideOption = (need_norm_inp) ? false : true;

    trans_obj.build(alpha_inp,gamma_inp,need_norm_inp);

    arums_G.build(nbX,nbY);
    arums_H.build(nbY,nbX);
    //
    NTU = true;
}

// general arums input
template<typename Tg, typename Th, typename Tm>
void 
dse<Tg,Th,Tm>::build_NTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const Tg& arums_G_inp, const Th& arums_H_inp, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    outsideOption = (need_norm_inp) ? false : true;

    trans_obj.build(alpha_inp,gamma_inp,need_norm_inp);

    arums_G = arums_G_inp;
    arums_H = arums_H_inp;
    //
    NTU = true;
}

// general simulation
template<typename Tg, typename Th, typename Tm> template<typename Ta, typename Tb>
void 
dse<Tg,Th,Tm>::build_NTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, Ta arums_G_inp, Tb arums_H_inp, int nbDraws, int seed, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    outsideOption = (need_norm_inp) ? false : true;

    trans_obj.build(alpha_inp,gamma_inp,need_norm_inp);

    arums_G_inp.simul(arums_G,nbDraws,seed);
    arums_H_inp.simul(arums_H,nbDraws,seed);
    //
    NTU = true;
}

// arums none
template<typename Tg, typename Th, typename Tm>
void 
dse<Tg,Th,Tm>::build_LTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    outsideOption = (need_norm_inp) ? false : true;

    trans_obj.build_LTU(lambda_inp,phi_inp);

    arums_G.build(nbX,nbY);
    arums_H.build(nbY,nbX);
    //
    LTU = true;
}

// general arums input
template<typename Tg, typename Th, typename Tm>
void 
dse<Tg,Th,Tm>::build_LTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, const Tg& arums_G_inp, const Th& arums_H_inp, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    outsideOption = (need_norm_inp) ? false : true;

    trans_obj.build_LTU(lambda_inp,phi_inp);

    arums_G = arums_G_inp;
    arums_H = arums_H_inp;
    //
    LTU = true;
}

// general simulation
template<typename Tg, typename Th, typename Tm> template<typename Ta, typename Tb>
void 
dse<Tg,Th,Tm>::build_LTU(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp, Ta arums_G_inp, Tb arums_H_inp, int nbDraws, int seed, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    outsideOption = (need_norm_inp) ? false : true;

    trans_obj.build_LTU(lambda_inp,phi_inp);

    arums_G_inp.simul(arums_G,nbDraws,seed);
    arums_H_inp.simul(arums_H,nbDraws,seed);
    //
    LTU = true;
}

template<typename Tg, typename Th, typename Tm>
void 
dse<Tg,Th,Tm>::trans()
{
    int nbX_temp = nbX;
    nbX = nbY;
    nbY = nbX_temp;

    arma::vec n_temp = n;
    n = m;
    m = n_temp;
    // Keith: fill in normalization later

    //mmf_obj.trans();
    trans_obj.trans();

    Tg arums_G_temp = arums_G;
    arums_G = arums_H;
    arums_H = arums_G_temp;
    //
}

template<typename Tg, typename Th, typename Tm>
bool 
dse<Tg,Th,Tm>::solve(arma::mat& mu_sol)
{
    bool res = this->solve(mu_sol,NULL);
    //
    return res;
}

template<typename Tg, typename Th, typename Tm>
bool 
dse<Tg,Th,Tm>::solve(arma::mat& mu_sol, const char* solver)
{
    bool res = false;
    const char sig = (solver != NULL) ? solver[0] : char(0);
    
    if (solver) { // not NULL
        if (sig=='a') {
            res = arc_newton(*this,mu_sol);
        }
        /*if (sig=='c') { // only works with empirical case
            res = cupids_lp(*this,mu_sol);
        }*/
        if (sig=='d') {
            res = darum(*this,mu_sol);
        }
        if (sig=='e') {
            res = eap_nash(*this,mu_sol);
        }
        if (sig=='j') {
            res = jacobi(*this,mu_sol);
        }
        if (sig=='m') {
            res = max_welfare(*this,mu_sol);
        }
        if (sig=='o') {
            res = oap_lp(*this,mu_sol);
        }
    }/* else { // default
        if (NTU) {
            res = darum(*this,mu_sol);
        } else if (TU) {
            res = max_welfare(*this,mu_sol);
        } else {
            res = jacobi(*this,mu_sol);
        }
    }*/
    //
    return res;
}

template<typename Tg, typename Th, typename Tm>
bool 
dse<Tg,Th,Tm>::solve(arma::mat& mu_sol, arma::mat& U_out, arma::mat& V_out, const char* solver)
{
    bool res = false;
    const char sig = (solver != NULL) ? solver[0] : char(0);
    
    if (solver) { // not NULL
        if (sig=='a') {
            res = arc_newton(*this,mu_sol,U_out,V_out);
        }
        /*if (sig=='c') { // only works with empirical case
            res = cupids_lp(*this,mu_sol);
        }*/
        if (sig=='d') {
            res = darum(*this,mu_sol,U_out,V_out);
        }
        if (sig=='e') {
            res = eap_nash(*this,mu_sol,U_out,V_out);
        }
        if (sig=='j') {
            res = jacobi(*this,mu_sol,U_out,V_out);
        }
        if (sig=='m') {
            res = max_welfare(*this,mu_sol,U_out,V_out);
        }
        if (sig=='o') {
            res = oap_lp(*this,mu_sol,U_out,V_out);
        }
    }/* else { // default
        if (NTU) {
            res = darum(*this,mu_sol,U_out,V_out);
        } else if (TU) {
            res = max_welfare(*this,mu_sol,U_out,V_out);
        } else {
            res = jacobi(*this,mu_sol,U_out,V_out);
        }
    }*/
    //
    return res;
}

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
 * Demand-Supply Equilibrium (DSE) market
 *
 * Keith O'Hara
 * 08/08/2016
 */

template <class Ta>
class DSE
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

        double sigma;

        arma::vec n;
        arma::vec m;

        mmf mmf_obj;
        transfers trans_obj;

        //none arums_G_none;
        //none arums_H_none;

        //empirical arums_G_empirical;
        //empirical arums_H_empirical;

        Ta arums_G;
        Ta arums_H;

        // member functions
        void build_TU(arma::vec n_inp, arma::vec m_inp, arma::mat phi, bool need_norm_inp);
        template <typename T> void build_TU(arma::vec n_inp, arma::vec m_inp, arma::mat phi, T arums_G, T arums_H, int nbDraws, int seed, bool need_norm_inp);

        void build_NTU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha, arma::mat gamma, bool need_norm_inp);
        template <typename T> void build_NTU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha, arma::mat gamma, T arums_G, T arums_H, int nbDraws, int seed, bool need_norm_inp);

        void build_LTU(arma::vec n_inp, arma::vec m_inp, arma::mat lambda, arma::mat phi, bool need_norm_inp);
        template <typename T> void build_LTU(arma::vec n_inp, arma::vec m_inp, arma::mat lambda, arma::mat phi, T arums_G, T arums_H, int nbDraws, int seed, bool need_norm_inp);

        void trans();

    private:
        bool arum_none;
        bool arum_empirical;
        bool arum_general; // need to finish this later
};

// arums none
void DSE::build_TU(arma::vec n_inp, arma::vec m_inp, arma::mat phi, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = true;
    }

    trans_obj.build_TU(phi);

    arums_G.build(nbX,nbY);
    arums_H.build(nbY,nbX);
    //
    arum_none = true;
}

// need templates to allow for general arums type
template <typename T>
void DSE::build_TU(arma::vec n_inp, arma::vec m_inp, arma::mat phi, T arums_G_inp, T arums_H_inp, int nbDraws, int seed, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = true;
    }

    trans_obj.build_TU(phi);

    arums_G_inp.simul(arums_G,nbDraws,seed);
    arums_H_inp.simul(arums_H,nbDraws,seed);
    //
    arum_empirical = true;
}

// 'general' class output requires complicated polymorphism for class structure; finish later

// arums none
void DSE::build_NTU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha, arma::mat gamma, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = true;
    }

    trans_obj.build_NTU(alpha,gamma);

    arums_G.build(nbX,nbY);
    arums_H.build(nbY,nbX);
    //
    arum_none = true;
}

// general simulation
template <typename T>
void DSE::build_NTU(arma::vec n_inp, arma::vec m_inp, arma::mat alpha, arma::mat gamma, T arums_G, T arums_H, int nbDraws, int seed, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = true;
    }

    trans_obj.build_NTU(alpha,gamma);

    arums_G.simul(arums_G_empirical,nbDraws,seed);
    arums_H.simul(arums_H_empirical,nbDraws,seed);
    //
    arum_empirical = true;
}

// arums none
void DSE::build_LTU(arma::vec n_inp, arma::vec m_inp, arma::mat lambda, arma::mat phi, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = true;
    }

    trans_obj.build_LTU(lambda,phi);

    arums_G.build(nbX,nbY);
    arums_H.build(nbY,nbX);
    //
    arum_none = true;
}

// general simulation
template <typename T>
void DSE::build_LTU(arma::vec n_inp, arma::vec m_inp, arma::mat lambda, arma::mat phi, T arums_G_inp, T arums_H_inp, int nbDraws, int seed, bool need_norm_inp)
{
    need_norm = need_norm_inp;

    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    if (need_norm) {
        outsideOption = false;
    } else {
        outsideOption = true;
    }

    trans_obj.build_LTU(lambda,phi);

    arums_G_inp.simul(arums_G,nbDraws,seed);
    arums_H_inp.simul(arums_H,nbDraws,seed);
    //
    arum_empirical = true;
}

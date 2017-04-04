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
 * Matching Function Equilibrium (MFE) market
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 03/14/2017
 */

// short build function, mmf_obj not touched
template<typename Tt>
void 
mfe<Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    need_norm = false;

    outsideOption = true;
    sigma = 1.0;
}

template<typename Tt>
void 
mfe<Tt>::trans()
{
    int nbX_temp = nbX;

    nbX = nbY;
    nbY = nbX_temp;
    //
    arma::vec n_temp = n;
    n = m;
    m = n_temp;
    // Keith: fill in normalization later

    mmfs_obj.trans();
    //
}

template<typename Tt>
arma::vec 
mfe<Tt>::marg_x_inv(const arma::mat& B_ys)
const
{
    arma::vec ret = this->marg_x_inv(B_ys,NULL);
    //
    return ret;
}

template<typename Tt>
arma::vec 
mfe<Tt>::marg_x_inv(const arma::mat& B_ys, arma::uvec* xs)
const
{
    arma::uvec temp_ind = (xs) ? *xs : uvec_linspace(0, (int) nbX-1);
    bool coeff = (need_norm) ? false : true;
    arma::vec ubs(nbX);

    if (need_norm) {
        ubs.fill(1E10);
    } else {
        ubs = n;
    }
    //
    trame_mfe_zeroin_data<Tt> root_data;

    root_data.mfe_obj = *this;
    root_data.coeff = coeff;
    root_data.B_ys  = B_ys;
    //
    arma::vec the_a_xs(temp_ind.n_elem);

    for (int j=0; j < (int) temp_ind.n_elem; j++) {
        int x = temp_ind(j);
        root_data.x_ind = x;

        the_a_xs(j) = zeroin(0.0, ubs(x), marg_x_inv_fn, &root_data, NULL, NULL);
    }
    //
    return the_a_xs;
}

template<typename Tt>
arma::vec 
mfe<Tt>::marg_y_inv(const arma::mat& A_xs)
const
{
    arma::vec ret = this->marg_y_inv(A_xs,NULL);
    //
    return ret;
}

template<typename Tt>
arma::vec 
mfe<Tt>::marg_y_inv(const arma::mat& A_xs, arma::uvec* ys)
const
{
    arma::uvec temp_ind = (ys) ? *ys : uvec_linspace(0, nbY-1);
    bool coeff = (need_norm) ? false : true;
    arma::vec ubs(nbY);

    if (need_norm) {
        ubs.fill(1E10);
    } else {
        ubs = m;
    }
    //
    trame_mfe_zeroin_data<Tt> root_data;

    root_data.mfe_obj = *this;
    root_data.coeff = coeff;
    root_data.A_xs  = A_xs;
    //
    arma::vec the_b_ys(temp_ind.n_elem);
        
    for (int j=0; j < (int) temp_ind.n_elem; j++) {
        int y = temp_ind(j);
        root_data.y_ind = y;

        the_b_ys(j) = zeroin(0.0, ubs(y), marg_y_inv_fn, &root_data, NULL, NULL);
    }
    //
    return the_b_ys;
}

// solve

template<typename Tt>
bool 
mfe<Tt>::solve(arma::mat& mu_sol)
{
    bool res = ipfp(*this,mu_sol);
    //
    return res;
}

template<typename Tt>
bool 
mfe<Tt>::solve(arma::mat& mu_sol, const char* solver)
{
    bool res = false;
    const char sig = (solver != NULL) ? solver[0] : char(0);
    
    if (solver) { // not NULL
        if (sig=='i') {
            res = ipfp(*this,mu_sol);
        }
        if (sig=='n') {
            //res = nodal_newton(*this,mu_sol);
        }
    } else {
        res = ipfp(*this,mu_sol);
    }
    //
    return res;
}

template<typename Tt>
bool 
mfe<Tt>::solve(arma::mat& mu_sol, arma::mat& U_out, arma::mat& V_out, const char* solver)
{
    bool res = false;
    const char sig = (solver != NULL) ? solver[0] : char(0);
    
    if (solver) { // not NULL
        if (sig=='i') {
            res = ipfp(*this,mu_sol,U_out,V_out);
        }
        if (sig=='n') {
            //res = nodal_newton(*this,mu_sol,U_out,V_out);
        }
    } else {
        res = ipfp(*this,mu_sol,U_out,V_out);
    }
    //
    return res;
}

// root finding functions

template<typename Tt>
double
mfe<Tt>::marg_x_inv_fn(double z, void* opt_data)
{
    trame_mfe_zeroin_data<Tt> *d = reinterpret_cast<trame_mfe_zeroin_data<Tt>*>(opt_data);
    //
    arma::uvec x_ind_temp(1);
    x_ind_temp(0) = d->x_ind;
    
    double term_1 = (d->coeff) ? z : 0;
    //
    double ret = term_1 - d->mfe_obj.n(d->x_ind) + arma::accu(d->mfe_obj.mmfs_obj.M(z,d->B_ys,&x_ind_temp,NULL));
    //
    return ret;
}

template<typename Tt>
double
mfe<Tt>::marg_y_inv_fn(double z, void* opt_data)
{
    trame_mfe_zeroin_data<Tt> *d = reinterpret_cast<trame_mfe_zeroin_data<Tt>*>(opt_data);
    //
    arma::uvec y_ind_temp(1);
    y_ind_temp(0) = d->y_ind;
    
    double term_1 = (d->coeff) ? z : 0;
    //
    double ret = term_1 - d->mfe_obj.m(d->y_ind) + arma::accu(d->mfe_obj.mmfs_obj.M(d->A_xs,z,NULL,&y_ind_temp));
    //
    return ret;
}

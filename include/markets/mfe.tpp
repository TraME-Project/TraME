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

//
// market transpose

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
inline
arma::vec
mfe<Tt>::marg_x_inv(const arma::mat& B_ys)
const
{
    return this->marg_x_inv(B_ys,NULL);
}

template<typename Tt>
inline
arma::vec
mfe<Tt>::marg_x_inv(const arma::mat& B_ys, const arma::uvec* xs)
const
{
    arma::uvec index_vec = (xs) ? *xs : uvec_linspace(0, (int) nbX-1);
    bool coeff = (need_norm) ? false : true;

    arma::vec ubs(nbX);
    (need_norm) ? ubs.fill(1E10) : ubs = n;
    //
    trame_mfe_zeroin_data<Tt> root_data;

    root_data.mfe_obj = *this;
    root_data.coeff = coeff;
    root_data.B_ys  = B_ys;
    //
    arma::vec the_a_xs(index_vec.n_elem);

    for (int j=0; j < (int) index_vec.n_elem; j++) {
        int x = index_vec(j);
        root_data.x_ind = x;

        the_a_xs(j) = zeroin(0.0, ubs(x), marg_x_inv_fn, &root_data, NULL, NULL);
    }
    //
    return the_a_xs;
}

//
// ipfp

template<typename Tt>
inline
arma::vec
mfe<Tt>::marg_y_inv(const arma::mat& A_xs)
const
{
    return this->marg_y_inv(A_xs,NULL);
}

template<typename Tt>
inline
arma::vec
mfe<Tt>::marg_y_inv(const arma::mat& A_xs, const arma::uvec* ys)
const
{
    arma::uvec index_vec = (ys) ? *ys : uvec_linspace(0, nbY-1);
    bool coeff = (need_norm) ? false : true;
    arma::vec ubs(nbY);

    (need_norm) ? ubs.fill(1E10) : ubs = m;
    //
    trame_mfe_zeroin_data<Tt> root_data;

    root_data.mfe_obj = *this;
    root_data.coeff = coeff;
    root_data.A_xs  = A_xs;
    //
    arma::vec the_b_ys(index_vec.n_elem);
        
    for (int j=0; j < (int) index_vec.n_elem; j++) {
        int y = index_vec(j);
        root_data.y_ind = y;

        the_b_ys(j) = zeroin(0.0, ubs(y), marg_y_inv_fn, &root_data, NULL, NULL);
    }
    //
    return the_b_ys;
}

template<>
inline
arma::vec
mfe<mmfs::min>::marg_x_inv(const arma::mat& B_ys, const arma::uvec* xs)
const
{
    arma::uvec index_vec = (xs) ? *xs : uvec_linspace(0, nbX-1);
    //
    arma::vec a_NTU = n.elem(index_vec);
    arma::mat B_NTU = arma::trans( elem_prod(arma::trans(mmfs_obj.aux_gamma_exp.rows(index_vec)/mmfs_obj.aux_alpha_exp.rows(index_vec)), B_ys) );
    arma::mat C_NTU = mmfs_obj.aux_alpha_exp.rows(index_vec);

    arma::vec the_a_xs = inv_pwa(a_NTU, B_NTU, C_NTU, 1.0);
    //
    return the_a_xs;
}

template<>
inline
arma::vec
mfe<mmfs::geo>::marg_x_inv(const arma::mat& B_ys, const arma::uvec* xs)
const
{
    arma::uvec index_vec = (xs) ? *xs : uvec_linspace(0, nbX-1);
    //
    arma::mat sqrt_A_xs;
    arma::mat sqrt_B_ys = arma::sqrt(B_ys);

    if (!need_norm) {
        arma::mat b = (mmfs_obj.aux_phi_exp.rows(index_vec) * sqrt_B_ys) / 2;
        sqrt_A_xs = arma::sqrt(n.rows(index_vec) + b%b) - b;
    } else{
        sqrt_A_xs = n.elem(index_vec) / arma::vectorise(mmfs_obj.aux_phi_exp.rows(index_vec) * sqrt_B_ys);
    }
        
    arma::vec the_a_xs = arma::vectorise(sqrt_A_xs % sqrt_A_xs);
    //
    return the_a_xs;
}

template<>
inline
arma::vec
mfe<mmfs::min>::marg_y_inv(const arma::mat& A_xs, const arma::uvec* ys)
const
{
    arma::uvec index_vec = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::vec a_NTU = m.elem(index_vec);
    arma::mat B_NTU = arma::trans( elem_prod(mmfs_obj.aux_alpha_exp.cols(index_vec)/mmfs_obj.aux_gamma_exp.cols(index_vec), A_xs) );
    arma::mat C_NTU = arma::trans(mmfs_obj.aux_gamma_exp.cols(index_vec));

    arma::vec the_b_ys = inv_pwa(a_NTU, B_NTU, C_NTU, 1.0);
    //
    return the_b_ys;
}

template<>
inline
arma::vec
mfe<mmfs::geo>::marg_y_inv(const arma::mat& A_xs, const arma::uvec* ys)
const
{
    arma::uvec index_vec = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat sqrt_B_ys;
    arma::mat sqrt_A_xs = arma::sqrt(A_xs);

    if (!need_norm) {
        arma::mat b = arma::trans(sqrt_A_xs.t() * mmfs_obj.aux_phi_exp.cols(index_vec)) / 2; // not sure about this
        sqrt_B_ys = arma::sqrt(m.rows(index_vec) + b%b) - b;
    } else {
        sqrt_B_ys = m.elem(index_vec) / arma::vectorise(arma::trans(sqrt_A_xs.t() * mmfs_obj.aux_phi_exp.cols(index_vec))); // not sure about this
    }
        
    arma::vec the_b_ys = arma::vectorise(sqrt_B_ys % sqrt_B_ys);
    //
    return the_b_ys;
}

//
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

//
// root finding functions

template<typename Tt>
inline
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
inline
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

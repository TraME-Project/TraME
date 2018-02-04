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
 * Matching Function Equilibrium (MFE) market
 *
 * Keith O'Hara
 * 08/16/2016
 *
 * This version:
 * 02/04/2018
 */

template<typename Tt>
mfe<Tt>::mfe(const arma::vec& n_inp, const arma::vec& m_inp)
{   
    this->build(n_inp,m_inp);
}

template<typename Tt>
mfe<Tt>::mfe(const double sigma_inp, const bool need_norm_inp)
{   
    this->build(sigma_inp,need_norm_inp);
}

template<typename Tt>
mfe<Tt>::mfe(const arma::vec& n_inp, const arma::vec& m_inp, const double sigma_inp, const bool need_norm_inp)
{   
    this->build(n_inp,m_inp,sigma_inp,need_norm_inp);
}

// short build functions, mmf_obj not touched

template<typename Tt>
void
mfe<Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;
}

template<typename Tt>
void
mfe<Tt>::build(const double sigma_inp, const bool need_norm_inp)
{
    nbX = 0;
    nbY = 0;

    sigma = sigma_inp;

    need_norm = need_norm_inp;
    outside_option = (need_norm_inp) ? false : true;
}

template<typename Tt>
void
mfe<Tt>::build(const arma::vec& n_inp, const arma::vec& m_inp, const double sigma_inp, const bool need_norm_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    sigma = sigma_inp;

    need_norm = need_norm_inp;
    outside_option = (need_norm_inp) ? false : true;
}

//
// specialized builds

template<>
inline
void
mfe<mmfs::ces>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp, const arma::mat& tau_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    mmfs_obj.build(alpha_inp,gamma_inp,tau_inp,need_norm);
}

template<>
inline
void
mfe<mmfs::cd>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& lambda_inp, const arma::mat& phi_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;
    
    mmfs_obj.build(lambda_inp,phi_inp/sigma,need_norm);
}

template<>
inline
void
mfe<mmfs::min>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& alpha_inp, const arma::mat& gamma_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    mmfs_obj.build(alpha_inp/sigma,gamma_inp/sigma,need_norm);
}

template<>
inline
void
mfe<mmfs::geo>::build(const arma::vec& n_inp, const arma::vec& m_inp, const arma::mat& phi_inp)
{
    nbX = n_inp.n_elem;
    nbY = m_inp.n_elem;
    
    n = n_inp;
    m = m_inp;

    mmfs_obj.build(phi_inp/sigma,need_norm);
}

//
// market transpose

template<typename Tt>
void
mfe<Tt>::trans()
{
    std::swap(nbX,nbY); 
    //
    n.swap(m);

    mmfs_obj.trans();
}

//
// ipfp

template<typename Tt>
arma::vec
mfe<Tt>::marg_x_inv(const arma::mat& B_ys)
const
{
    return this->marg_x_inv(B_ys,nullptr);
}

template<typename Tt>
arma::vec
mfe<Tt>::marg_x_inv(const arma::mat& B_ys, const arma::uvec* xs)
const
{
    const arma::uvec index_vec = (xs) ? *xs : uvec_linspace(0, static_cast<int>(nbX - 1));
    const bool coeff = (need_norm) ? false : true;

    arma::vec ubs(nbX);
    (need_norm) ? ubs.fill(1E10) : ubs = n;

    //

    trame_mfe_zeroin_data<Tt> root_data;

    root_data.mfe_obj = *this;
    root_data.coeff = coeff;
    root_data.B_ys  = B_ys;

    //

    arma::vec the_a_xs(index_vec.n_elem);

#ifdef TRAME_USE_OMP
    #pragma omp parallel for firstprivate(root_data)
#endif
    for (uint_t j=0; j < index_vec.n_elem; j++)
    {
        uint_t x = index_vec(j);
        root_data.x_ind = x;

        the_a_xs(j) = zeroin(0.0, ubs(x), marg_x_inv_fn, &root_data, nullptr, nullptr);
    }
    //
    return the_a_xs;
}

template<typename Tt>
arma::vec
mfe<Tt>::marg_y_inv(const arma::mat& A_xs)
const
{
    return this->marg_y_inv(A_xs,nullptr);
}

template<typename Tt>
arma::vec
mfe<Tt>::marg_y_inv(const arma::mat& A_xs, const arma::uvec* ys)
const
{
    const arma::uvec index_vec = (ys) ? *ys : uvec_linspace(0, nbY-1);
    const bool coeff = (need_norm) ? false : true;

    arma::vec ubs(nbY);
    (need_norm) ? ubs.fill(1E10) : ubs = m;
    //
    trame_mfe_zeroin_data<Tt> root_data;

    root_data.mfe_obj = *this;
    root_data.coeff = coeff;
    root_data.A_xs  = A_xs;
    //
    arma::vec the_b_ys(index_vec.n_elem);
    
#ifdef TRAME_USE_OMP
    #pragma omp parallel for firstprivate(root_data)
#endif
    for (uint_t j=0; j < index_vec.n_elem; j++)
    {
        uint_t y = index_vec(j);
        root_data.y_ind = y;

        the_b_ys(j) = zeroin(0.0, ubs(y), marg_y_inv_fn, &root_data, nullptr, nullptr);
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
    const arma::uvec index_vec = (xs) ? *xs : uvec_linspace(0, nbX-1);
    //
    const arma::vec a_NTU = n.elem(index_vec);
    const arma::mat B_NTU = arma::trans( elem_prod(arma::trans(mmfs_obj.aux_gamma_exp.rows(index_vec)/mmfs_obj.aux_alpha_exp.rows(index_vec)), B_ys) );
    const arma::mat C_NTU = mmfs_obj.aux_alpha_exp.rows(index_vec);

    const arma::vec the_a_xs = inv_pwa(a_NTU, B_NTU, C_NTU, 1.0);
    //
    return the_a_xs;
}

template<>
inline
arma::vec
mfe<mmfs::geo>::marg_x_inv(const arma::mat& B_ys, const arma::uvec* xs)
const
{
    const arma::uvec index_vec = (xs) ? *xs : uvec_linspace(0, nbX-1);
    //
    arma::mat sqrt_A_xs;
    const arma::mat sqrt_B_ys = arma::sqrt(B_ys);

    if (!need_norm) {
        const arma::mat b = (mmfs_obj.aux_phi_exp.rows(index_vec) * sqrt_B_ys) / 2.0;
        sqrt_A_xs = arma::sqrt(n.rows(index_vec) + b%b) - b;
    } else{
        sqrt_A_xs = n.elem(index_vec) / arma::vectorise(mmfs_obj.aux_phi_exp.rows(index_vec) * sqrt_B_ys);
    }
        
    return arma::vectorise(sqrt_A_xs % sqrt_A_xs);
}

template<>
inline
arma::vec
mfe<mmfs::min>::marg_y_inv(const arma::mat& A_xs, const arma::uvec* ys)
const
{
    const arma::uvec index_vec = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    const arma::vec a_NTU = m.elem(index_vec);
    const arma::mat B_NTU = arma::trans( elem_prod(mmfs_obj.aux_alpha_exp.cols(index_vec)/mmfs_obj.aux_gamma_exp.cols(index_vec), A_xs) );
    const arma::mat C_NTU = arma::trans(mmfs_obj.aux_gamma_exp.cols(index_vec));

    const arma::vec the_b_ys = inv_pwa(a_NTU, B_NTU, C_NTU, 1.0);
    //
    return the_b_ys;
}

template<>
inline
arma::vec
mfe<mmfs::geo>::marg_y_inv(const arma::mat& A_xs, const arma::uvec* ys)
const
{
    const arma::uvec index_vec = (ys) ? *ys : uvec_linspace(0, nbY-1);
    //
    arma::mat sqrt_B_ys;
    const arma::mat sqrt_A_xs = arma::sqrt(A_xs);

    if (!need_norm) {
        const arma::mat b = arma::trans(sqrt_A_xs.t() * mmfs_obj.aux_phi_exp.cols(index_vec)) / 2.0; // not sure about this
        sqrt_B_ys = arma::sqrt(m.rows(index_vec) + b%b) - b;
    } else {
        sqrt_B_ys = m.elem(index_vec) / arma::vectorise(arma::trans(sqrt_A_xs.t() * mmfs_obj.aux_phi_exp.cols(index_vec))); // not sure about this
    }
        
    return arma::vectorise(sqrt_B_ys % sqrt_B_ys);
}

//
// solve

template<typename Tt>
bool
mfe<Tt>::solve(arma::mat& mu_sol)
{
    return ipfp(*this,mu_sol);
}

template<typename Tt>
bool
mfe<Tt>::solve(arma::mat& mu_sol, const char* solver)
{
    bool res = false;
    const char sig = (solver) ? solver[0] : char(0);
    
    if (solver) { // not nullptr
        if (sig=='i') {
            res = ipfp(*this,mu_sol);
        }
        if (sig=='n') {
            res = nodal_newton(*this,mu_sol);
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
    const char sig = (solver) ? solver[0] : char(0);
    
    if (solver) { // not nullptr
        if (sig=='i') {
            res = ipfp(*this,mu_sol,U_out,V_out);
        }
        if (sig=='n') {
            res = nodal_newton(*this,mu_sol,U_out,V_out);
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
double
mfe<Tt>::marg_x_inv_fn(const double z, void* opt_data)
{
    trame_mfe_zeroin_data<Tt> *d = reinterpret_cast<trame_mfe_zeroin_data<Tt>*>(opt_data);
    //
    arma::uvec x_ind_temp(1);
    x_ind_temp(0) = d->x_ind;
    //
    const double term_1 = (d->coeff) ? z : 0;
    
    const double ret = term_1 - d->mfe_obj.n(d->x_ind) + arma::accu(d->mfe_obj.mmfs_obj.M(z,d->B_ys,&x_ind_temp,nullptr));
    //
    return ret;
}

template<typename Tt>
double
mfe<Tt>::marg_y_inv_fn(const double z, void* opt_data)
{
    trame_mfe_zeroin_data<Tt> *d = reinterpret_cast<trame_mfe_zeroin_data<Tt>*>(opt_data);
    //
    arma::uvec y_ind_temp(1);
    y_ind_temp(0) = d->y_ind;
    //
    const double term_1 = (d->coeff) ? z : 0;
    
    const double ret = term_1 - d->mfe_obj.m(d->y_ind) + arma::accu(d->mfe_obj.mmfs_obj.M(d->A_xs,z,nullptr,&y_ind_temp));
    //
    return ret;
}

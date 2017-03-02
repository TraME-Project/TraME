/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
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
 * mmf class
 *
 * Keith O'Hara
 * 08/15/2016
 */

#include "trame.hpp"

void trame::mmf::build_ETU(const arma::vec& n_ETU, const arma::vec& m_ETU, const arma::mat& C_ETU, const arma::mat& D_ETU, const arma::mat& kappa_ETU, bool need_norm_ETU)
{
    n = n_ETU;
    m = m_ETU;

    C = C_ETU;
    D = D_ETU;

    aux_log_C = arma::log(C);
    aux_log_D = arma::log(D);

    kappa = kappa_ETU;
    
    need_norm = need_norm_ETU;

    ETU = true;
}

void trame::mmf::build_LTU(const arma::vec& n_LTU, const arma::vec& m_LTU, const arma::mat& lambda_LTU, const arma::mat& K_LTU, bool need_norm_LTU)
{
    n = n_LTU;
    m = m_LTU;
    K = K_LTU;

    lambda = lambda_LTU;
    aux_zeta = 1 - lambda_LTU;
    
    need_norm = need_norm_LTU;

    LTU = true;
}

void trame::mmf::build_NTU(const arma::vec& n_NTU, const arma::vec& m_NTU, const arma::mat& A_NTU, const arma::mat& B_NTU, bool need_norm_NTU)
{
    n = n_NTU;
    m = m_NTU;
    
    A = A_NTU;
    B = B_NTU;
    
    need_norm = need_norm_NTU;

    NTU = true;
}

void trame::mmf::build_TU(const arma::vec& n_TU, const arma::vec& m_TU, const arma::mat& K_TU, bool need_norm_TU)
{
    n = n_TU;
    m = m_TU;
    K = K_TU;

    need_norm = need_norm_TU;

    TU = true;
}

arma::mat trame::mmf::M(const arma::mat& a_xs, const arma::mat& b_ys)
const
{
    arma::mat ret = this->M(a_xs,b_ys,NULL,NULL);
    //
    return ret;
}

arma::mat trame::mmf::M(const arma::mat& a_xs, const arma::mat& b_ys, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::mat ret;
    
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, (int) n.n_elem-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, (int) m.n_elem-1);
    //
    if (ETU) {
        arma::mat term_1 = arma::exp(aux_log_C(x_ind,y_ind) + elem_prod(kappa(x_ind,y_ind), arma::log(a_xs)));
        arma::mat term_2 = arma::exp(aux_log_D(x_ind,y_ind) + arma::trans(arma::trans(kappa(x_ind,y_ind)) % arma::log(b_ys)));

        ret = arma::exp((1/kappa(x_ind,y_ind)) % arma::log((term_1+term_2)/2)); 
    }
    //
    if (LTU) {
        arma::mat term_1 = arma::exp(elem_prod(lambda(x_ind,y_ind), arma::log(a_xs)));
        arma::mat term_2 = arma::trans(arma::exp( elem_prod(arma::trans(aux_zeta(x_ind,y_ind)), arma::log(b_ys)) ));
        arma::mat term_3 = K(x_ind,y_ind);

        ret = term_1 % term_2 % term_3;
    }
    //
    if (NTU) {
        arma::mat term_1 = elem_prod(a_xs,A(x_ind,y_ind));
        arma::mat term_2 = arma::trans( elem_prod(b_ys, arma::trans(A(x_ind,y_ind))) );

        ret = arma::min(term_1, term_2);
    }
    //
    if (TU) {
        arma::mat term_1 = K(x_ind,y_ind);
        arma::mat term_2 = arma::sqrt(a_xs * b_ys.t());

        ret = term_1 % term_2;
    }
    //
    return ret;
}

arma::mat trame::mmf::M(const double& a_xs, const arma::mat& b_ys, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::mat ret;

    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, (int) n.n_elem-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, (int) m.n_elem-1);
    //
    if (ETU) {
        arma::mat term_1 = arma::exp(aux_log_C(x_ind,y_ind) + kappa(x_ind,y_ind) * std::log(a_xs));
        arma::mat term_2 = arma::exp(aux_log_D(x_ind,y_ind) + arma::trans(arma::trans(kappa(x_ind,y_ind)) % arma::log(b_ys)));

        ret = arma::exp((1/kappa(x_ind,y_ind)) % arma::log((term_1+term_2)/2)); 
    }
    //
    if (LTU) {
        arma::mat term_1 = arma::exp(lambda(x_ind,y_ind) * std::log(a_xs));
        arma::mat term_2 = arma::trans(arma::exp(arma::trans(aux_zeta(x_ind,y_ind)) % arma::log(b_ys)));
        arma::mat term_3 = K(x_ind,y_ind);

        ret = term_1 % term_2 % term_3;
    }
    //
    if (NTU) {
        arma::mat term_1 = a_xs * A(x_ind,y_ind);
        arma::mat term_2 = arma::trans(b_ys % arma::trans(A(x_ind,y_ind)));

        ret = arma::min(term_1, term_2);
    }
    //
    if (TU) {
        arma::mat term_1 = K(x_ind,y_ind);
        arma::mat term_2 = arma::sqrt(a_xs * b_ys.t());

        ret = term_1 % term_2;
    }
    //
    return ret;
}

arma::mat trame::mmf::M(const arma::mat& a_xs, const double& b_ys, arma::uvec* xs, arma::uvec* ys)
const
{
    arma::mat ret;
    
    arma::uvec x_ind = (xs) ? *xs : uvec_linspace(0, (int) n.n_elem-1); 
    arma::uvec y_ind = (ys) ? *ys : uvec_linspace(0, (int) m.n_elem-1);
    //
    if (ETU) {
        arma::mat term_1 = arma::exp(aux_log_C(x_ind,y_ind) + kappa(x_ind,y_ind) % arma::log(a_xs));
        arma::mat term_2 = arma::exp(aux_log_D(x_ind,y_ind) + arma::trans(arma::trans(kappa(x_ind,y_ind)) * std::log(b_ys)));

        ret = arma::exp((1/kappa(x_ind,y_ind)) % arma::log((term_1+term_2)/2)); 
    }
    //
    if (LTU) {
        arma::mat term_1 = arma::exp(lambda(x_ind,y_ind) % arma::log(a_xs));
        arma::mat term_2 = arma::trans(arma::exp(arma::trans(aux_zeta(x_ind,y_ind)) * std::log(b_ys)));
        arma::mat term_3 = K(x_ind,y_ind);
        
        ret = term_1 % term_2 % term_3;
    }
    //
    if (NTU) {
        arma::mat term_1 = a_xs % A(x_ind,y_ind);
        arma::mat term_2 = arma::trans(b_ys * arma::trans(A(x_ind,y_ind)));

        ret = arma::min(term_1, term_2);
    }
    //
    if (TU) {
        arma::mat term_1 = K(x_ind,y_ind);
        arma::mat term_2 = arma::sqrt(a_xs * b_ys);

        ret = term_1 % term_2;
    }
    //
    return ret;
}

arma::mat trame::mmf::Mx0(const arma::mat& a_x)
const
{
    return a_x;
}

arma::mat trame::mmf::M0y(const arma::mat& b_y)
const
{
    return b_y;
}

void trame::mmf::trans()
{
    arma::vec n_temp = n;

    n = m;
    m = n_temp;
    //
    if (ETU) {
        arma::mat C_temp = C;

        C = D.t();
        D = C_temp.t();
        kappa = kappa.t();
    }
    
    if (LTU) {
        arma::mat lambda_temp;

        K = K.t();
        lambda = aux_zeta.t();
        aux_zeta = lambda_temp.t();
    }
    
    if (NTU) {
        arma::mat A_temp = A;

        A = B.t();
        B = A_temp.t();
    }
    
    if (TU) {
        K = K.t();
    }
    //
}

arma::vec trame::mmf::marg_x_inv(const arma::mat& B_ys)
{
    arma::vec ret = this->marg_x_inv(B_ys,NULL);
    //
    return ret;
}

arma::vec trame::mmf::marg_x_inv(const arma::mat& B_ys, arma::uvec* xs)
{
    arma::uvec temp_ind = (xs) ? *xs : uvec_linspace(0, (int) n.n_elem-1);
    //
    if (NTU) {
        arma::vec a_NTU = n.elem(temp_ind);
        arma::mat B_NTU = arma::trans( elem_prod(arma::trans(B.rows(temp_ind)/A.rows(temp_ind)), B_ys) );
        arma::mat C_NTU = A.rows(temp_ind);

        arma::vec the_a_xs = inv_pwa(a_NTU, B_NTU, C_NTU, 1.0);
        //
        return the_a_xs;
    } else if (TU) {
        arma::mat sqrt_A_xs;
        arma::mat sqrt_B_ys = arma::sqrt(B_ys);

        if (!need_norm) {
            arma::mat b = (K.rows(temp_ind) * sqrt_B_ys) / 2;
            sqrt_A_xs = arma::sqrt(n.rows(temp_ind) + b%b) - b;
        } else{
            sqrt_A_xs = n.elem(temp_ind) / arma::vectorise(K.rows(temp_ind) * sqrt_B_ys);
        }
        
        arma::vec the_a_xs = arma::vectorise(sqrt_A_xs % sqrt_A_xs);
        //
        return the_a_xs;
    } else { // 'default'
        int x;
        int nb_X = n.n_elem;

        bool coeff;
        arma::vec ubs(nb_X);

        if (need_norm) {
            coeff = false;
            ubs.fill(1E10);
        } else {
            coeff = true;
            ubs = n;
        }
        //
        trame_mmf_zeroin_data root_data;

        root_data.mmf_obj = *this;

        root_data.coeff = coeff;
        root_data.B_ys  = B_ys;
        //
        arma::vec the_a_xs(temp_ind.n_elem);

        for (int j=0; j < (int) temp_ind.n_elem; j++) {
            x = temp_ind(j);
            root_data.x_ind = x;

            the_a_xs(j) = zeroin(0.0, ubs(x), marg_x_inv_fn, &root_data, NULL, NULL);
        }
        //
        return the_a_xs;
    }
}

double trame::mmf::marg_x_inv_fn(double z, void* opt_data)
{
    trame_mmf_zeroin_data *d = reinterpret_cast<trame_mmf_zeroin_data*>(opt_data);
    //
    arma::uvec x_ind_temp(1);
    x_ind_temp(0) = d->x_ind;
    
    double term_1 = (d->coeff) ? z : 0;
    //
    double ret = term_1 - d->mmf_obj.n(d->x_ind) + arma::accu(d->mmf_obj.M(z,d->B_ys,&x_ind_temp,NULL));
    //
    return ret;
}

arma::vec trame::mmf::marg_y_inv(const arma::mat& A_xs)
{
    arma::vec ret = this->marg_y_inv(A_xs,NULL);
    //
    return ret;
}

arma::vec trame::mmf::marg_y_inv(const arma::mat& A_xs, arma::uvec* ys)
{
    arma::uvec temp_ind;

    if (ys) {
        temp_ind = *ys;
    } else {
        temp_ind = uvec_linspace(0, (int) m.n_elem-1);
    }
    //
    if (NTU) {
        arma::vec a_NTU = m.elem(temp_ind);
        arma::mat B_NTU = arma::trans( elem_prod(A.cols(temp_ind)/B.cols(temp_ind), A_xs) );
        arma::mat C_NTU = arma::trans(B.cols(temp_ind));

        arma::vec the_b_ys = inv_pwa(a_NTU, B_NTU, C_NTU, 1.0);
        //
        return the_b_ys;
    } else if (TU) {
        arma::mat sqrt_B_ys;
        arma::mat sqrt_A_xs = arma::sqrt(A_xs);

        if (!need_norm) {
            arma::mat b = arma::trans(sqrt_A_xs.t() * K.cols(temp_ind)) / 2; // not sure about this
            sqrt_B_ys = arma::sqrt(m.rows(temp_ind) + b%b) - b;
        } else {
            sqrt_B_ys = m.elem(temp_ind) / arma::vectorise(arma::trans(sqrt_A_xs.t() * K.cols(temp_ind))); // not sure about this
        }
        
        arma::vec the_b_ys = arma::vectorise(sqrt_B_ys % sqrt_B_ys);
        //
        return the_b_ys;
    } else { // 'default'
        int nb_Y = m.n_elem;

        bool coeff;
        arma::vec ubs(nb_Y);

        if (need_norm) {
            coeff = false;
            ubs.fill(1E10);
        } else {
            coeff = true;
            ubs = m;
        }
        //
        trame_mmf_zeroin_data root_data;

        root_data.mmf_obj = *this;

        root_data.coeff = coeff;
        root_data.A_xs  = A_xs;
        //
        int j, y;
        arma::vec the_b_ys(temp_ind.n_elem);
        
        for (j=0; j < (int) temp_ind.n_elem; j++) {
            y = temp_ind(j);
            root_data.y_ind = y;

            the_b_ys(j) = zeroin(0.0, ubs(y), marg_y_inv_fn, &root_data, NULL, NULL);
        }
        //
        return the_b_ys;
    }
}

double trame::mmf::marg_y_inv_fn(double z, void* opt_data)
{
    trame_mmf_zeroin_data *d = reinterpret_cast<trame_mmf_zeroin_data*>(opt_data);
    //
    arma::uvec y_ind_temp(1);
    y_ind_temp(0) = d->y_ind;
    
    double term_1 = (d->coeff) ? z : 0;
    //
    double ret = term_1 - d->mmf_obj.m(d->y_ind) + arma::accu(d->mmf_obj.M(d->A_xs,z,NULL,&y_ind_temp));
    //
    return ret;
}

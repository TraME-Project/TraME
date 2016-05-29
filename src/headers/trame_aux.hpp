/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2016 Alfred Galichon and the TraME Team
  ##
  ##   This file is part of the R package TraME.
  ##
  ##   The R package TraME is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package TraME is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

extern "C" {
    void beta_cdf_(double* x, double* p, double* q, double* result);
	void beta_cdf_inv_(double* x, double* p, double* q, double* result);
    void beta_cdf_inv_int_(double* x, double* p, double* q, double* result);
    void beta_pdf_(double* x, double* p, double* q, double* result);
}

double pbeta (double x, double* fn_pars)
{
    double p, q, res;
    p = fn_pars[0];
    q = fn_pars[1];
    
    beta_cdf_(&x, &p, &q, &res);
    
    return res;
}

arma::vec pbeta (arma::vec x, double* fn_pars)
{
    int i;
    int N = x.n_elem;
    double p, q, x_i, res_temp;
    arma::vec res(N);
    
    p = fn_pars[0];
    q = fn_pars[1];
    
    for(i=0; i<N; i++){
        x_i = x(i);
        beta_cdf_(&x_i, &p, &q, &res_temp);
        res(i) = res_temp;
    }
    
    return res;
}

double qbeta (double x, double* fn_pars)
{
    double p, q, res;
    p = fn_pars[0];
    q = fn_pars[1];
    
    beta_cdf_inv_(&x, &p, &q, &res);
    
    return res;
}

arma::vec qbeta (arma::vec x, double* fn_pars)
{
    int i;
    int N = x.n_elem;
    double p, q, x_i, res_temp;
    arma::vec res(N);
    
    p = fn_pars[0];
    q = fn_pars[1];
    
    for(i=0; i<N; i++){
        x_i = x(i);
        beta_cdf_inv_(&x_i, &p, &q, &res_temp);
        res(i) = res_temp;
    }
    
    return res;
}

double iqbeta (double x, double* fn_pars)
{
    double p, q, res;
    p = fn_pars[0];
    q = fn_pars[1];
    
    beta_cdf_inv_int_(&x, &p, &q, &res);
    
    return res;
}

arma::vec iqbeta (arma::vec x, double* fn_pars)
{
    int i;
    int N = x.n_elem;
    double p, q, x_i, res_temp;
    arma::vec res(N);
    
    p = fn_pars[0];
    q = fn_pars[1];
    
    for(i=0; i<N; i++){
        x_i = x(i);
        beta_cdf_inv_int_(&x_i, &p, &q, &res_temp);
        res(i) = res_temp;
    }
    
    return res;
}

double dbeta (double x, double* fn_pars)
{
    double p, q, res;
    p = fn_pars[0];
    q = fn_pars[1];
    
    beta_pdf_(&x, &p, &q, &res);
    
    return res;
}

arma::vec dbeta (arma::vec x, double* fn_pars)
{
    int i;
    int N = x.n_elem;
    double p, q, x_i, res_temp;
    arma::vec res(N);
    
    p = fn_pars[0];
    q = fn_pars[1];
    
    for(i=0; i<N; i++){
        x_i = x(i);
        beta_pdf_(&x_i, &p, &q, &res_temp);
        res(i) = res_temp;
    }
    
    return res;
}
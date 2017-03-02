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
 * auxiliary functions for equilibrium solvers
 *
 * Keith O'Hara
 * 08/23/2016
 */

/*
 * w upper bound
 * used in 'jacobi' solver and arc_newton
 *
 * Keith O'Hara
 * 08/25/2016
 */

template<typename Ta>
arma::mat w_upper_bound(const dse<Ta>& market)
{
    int nbX = market.nbX;
    int nbY = market.nbY;

    arma::vec n = market.n;
    arma::vec m = market.m;

    Ta arums_G = market.arums_G;
    Ta arums_H = market.arums_H;

    transfers trans_obj = market.trans_obj;
    int transfers_type = trans_obj.transfers_type;
    //
    //bool success = false;
    int iter = 0;
    int max_iter = 1000;
    double Z_min_val = -10.0;
    
    int x;
    double k = 1.0;
    double mu_fill;
    arma::uvec x_ind(1);
    arma::vec mu_cond_x(nbY,1);
    arma::mat U_star_x, Z;

    arma::mat w(nbX,nbY);
    arma::mat U(nbX,nbY);
    arma::mat V(nbX,nbY);

    while (Z_min_val < 0 && iter < max_iter) {
        iter ++;

        for (x=0; x < nbX; x++) {
            x_ind(0) = x;

            if (transfers_type==1) {
                mu_fill = 1.0 / (std::pow(2.0,-k) + (double) nbY);
                mu_cond_x.fill(mu_fill);

                arums_G.Gstarx(mu_cond_x,U_star_x,x);
                U.row(x) = arma::trans(U_star_x);

                w.row(x) = trans_obj.WU(U_star_x.t(),&x_ind,NULL);
                V.row(x) = trans_obj.VW(w.row(x),&x_ind,NULL);
            } else if (transfers_type==2) {
                w.row(x).fill(std::pow(2,k));

                U.row(x) = trans_obj.UW(w.row(x),&x_ind,NULL);
                V.row(x) = trans_obj.VW(w.row(x),&x_ind,NULL);
            } else {
                printf("w_upper_bound error: unrecognized transfers_type");
                return w;
            }
        }
        //
        arums_G.U = U;
        arums_H.U = V.t();

        arums_G.G(n);
        arums_H.G(m);

        Z = arums_G.mu_sol - arma::trans(arums_H.mu_sol);
        Z_min_val = elem_min(Z);
        //
        k *= 2;
    }
    //
    return w;
}

/*
template<typename Ta>
bool max_welfare_nlopt(int n_pars, std::vector<double>& io_val, double& opt_val, double* lb, double* ub,
                       double (*opt_objfn)(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data),
                       trame_market_opt_data<Ta> opt_data)
{
    bool success = false;

    nlopt::opt opt_trame(nlopt::LD_LBFGS, n_pars);
    nlopt::result result;

    if (lb) {
        opt_trame.set_lower_bounds(*lb);
    }
    if (ub) {
        opt_trame.set_upper_bounds(*ub);
    }

    opt_trame.set_min_objective(*opt_objfn, &opt_data);
    
    opt_trame.set_xtol_rel(1e-8);
    opt_trame.set_ftol_rel(1e-15);
    opt_trame.set_maxeval(5000);

    double minf;
    try {
        result = opt_trame.optimize(io_val, minf);
    } catch(...) {
        printf("error in max_welfare optimization (using LBFGS);\n");
        printf("retrying with MMA\n");

        nlopt::opt opt_trame_2(nlopt::LD_MMA, n_pars);

        if (lb) {
            opt_trame_2.set_lower_bounds(*lb);
        }
        if (ub) {
            opt_trame_2.set_upper_bounds(*ub);
        }

        opt_trame_2.set_min_objective(*opt_objfn, &opt_data);
        
        opt_trame_2.set_xtol_rel(1e-8);
        opt_trame_2.set_ftol_rel(1e-15);
        opt_trame_2.set_maxeval(5000);

        result = opt_trame_2.optimize(io_val, minf);
    }

    if (result > 0) {
        opt_val = minf;
        success = true;
    }

    return success;
}
*/
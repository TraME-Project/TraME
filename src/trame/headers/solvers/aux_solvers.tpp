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
 * auxiliary functions for equilibrium solvers
 *
 * Keith O'Hara
 * 08/23/2016
 */

/*
 * w upper bound
 * used in 'jacobi' solver
 *
 * Keith O'Hara
 * 08/25/2016
 */

template<typename Ta>
arma::mat w_upper_bound(dse<Ta> market)
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

// for use in jacobi
template<typename Ta>
class trame_zeroin_data_tmp
{
    public:
        int x_ind;
        int y_ind;

        arma::vec n;
        arma::vec m;

        arma::mat U;
        arma::mat V;

        Ta arums_G;
        Ta arums_H;

        transfers trans_obj;
};

template<typename Ta>
double jacobi_zeroin_fn(double z, const trame_zeroin_data_tmp<Ta>& opt_data)
{
    double ret = 1.0;

    int x_ind = opt_data.x_ind;
    int y_ind = opt_data.y_ind;

    arma::mat U = opt_data.U;
    arma::mat V = opt_data.V;

    Ta arums_G = opt_data.arums_G;
    Ta arums_H = opt_data.arums_H;
    
    transfers trans_obj = opt_data.trans_obj;
    //
    U(x_ind,y_ind) = trans_obj.UW(z,x_ind,y_ind);
    V(x_ind,y_ind) = trans_obj.VW(z,x_ind,y_ind);
    
    arums_G.U = U;
    arums_H.U = V.t();

    arums_G.G(opt_data.n);
    arums_H.G(opt_data.m);

    arma::mat Z = arums_G.mu_sol - arma::trans(arums_H.mu_sol);
    ret = Z(x_ind,y_ind);
    //
    return ret;
}

template<typename Ta>
double zeroin_tmp(double ax, double bx, double (*f)(double x, const trame_zeroin_data_tmp<Ta>& opt_data), const trame_zeroin_data_tmp<Ta>& zeroin_data, double* tol_inp, int* max_iter_inp)
{
	double a,b,c;
	double fa;
	double fb;
	double fc;

	double tol;
    int max_iter;
    
    if (tol_inp) {
        tol = *tol_inp;
    } else {
        tol = 1E-12;
    }

    if (max_iter_inp) {
        max_iter = *max_iter_inp;
    } else {
        max_iter = 10000;
    }
		
	a = ax;  b = bx;  fa = (*f)(a,zeroin_data);  fb = (*f)(b,zeroin_data);
	c = a;   fc = fa;
	
	// check endpoints
	if(fa == 0.0){
		return b;
	}
	if(fb == 0.0){
		return b;
	}
	
	// otherwise begin iterations
	int iter = 0;
	
	double eps_temp = std::numeric_limits<double>::epsilon();
	double tol_act = 2*eps_temp*fabs(b) + tol/2;

	double p, q, prev_step, new_step;
	//register double t1,cb,t2;
	double t1,cb,t2; // Keith: register is deprecated as of C++-11
	
	new_step = (c - b)/2;
	
	while(fabs(new_step) > tol_act && iter < max_iter){
		iter++;
		prev_step = b-a;
			
		if( fabs(fc) < fabs(fb) ){
			a = b;  b = c;  c = a;
			fa=fb;  fb=fc;  fc=fa;
		}
		
		new_step = (c-b)/2;

		if( fabs(prev_step) >= tol_act	&& fabs(fa) > fabs(fb) ){
			
			cb = c - b;
						
			if( a==c ){
				t1 = fb/fa;
				p = cb*t1;
				q = 1.0 - t1;
			}else{
				q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
				p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
				q = (q-1.0) * (t1-1.0) * (t2-1.0);
			}
						
			if( p > 0.0 ){
				q = -q;
			}else{
				p = -p;
			}

			if( p < (0.75*cb*q-fabs(tol_act*q)/2) && p < fabs(prev_step*q/2) ){
				new_step = p/q;
			}
		}
				
		if( fabs(new_step) < tol_act ){
			if( new_step > 0.0 ){
				new_step = tol_act;
			}else{
				new_step = -tol_act;
			}
		}
						
		a = b;  fa = fb;
		b += new_step;  fb = (*f)(b,zeroin_data);
		if( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) ){
			c = a;  fc = fa;
		}
	}
	//
	return b;
}

template<typename Ta>
bool max_welfare_nlopt(int n_pars, std::vector<double>& io_val, double& opt_val, double* lb, double* ub,
                       double (*opt_objfn)(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data),
                       trame_market_opt_data<Ta> opt_data)
{
    bool success = false;

    nlopt::opt opt_trame(nlopt::LD_LBFGS, n_pars);

    if (lb) {
        opt_trame.set_lower_bounds(*lb);
    }
    if (ub) {
        opt_trame.set_upper_bounds(*ub);
    }

    opt_trame.set_min_objective(*opt_objfn, &opt_data);
    
    opt_trame.set_xtol_rel(1e-7);
    opt_trame.set_maxeval(5000);

    double minf;
    nlopt::result result = opt_trame.optimize(io_val, minf);

    if (result > 0) {
        opt_val = minf;
        success = true;
    }

    return success;
}

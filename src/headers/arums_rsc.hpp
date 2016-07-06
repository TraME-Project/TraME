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

//#include "trame_aux.hpp"

#include "omp.h" //  OpenMP

// RSC class
class RSC
{
    public:
        // build_logit objects
        int nbX;
        int nbY;
        int nbParams;
        bool outsideOption;
        
        double* dist_pars;
        
        arma::mat zeta;
        arma::mat aux_ord;
        
        arma::cube aux_Influence_lhs;
        arma::cube aux_Influence_rhs;
        
        arma::cube aux_DinvPsigma; 
        arma::cube aux_Psigma;
        
        double   (*aux_cdf_eps)(double x, double* dist_pars);
        double (*aux_quant_eps)(double x, double* dist_pars);
       
        double (*aux_pdf_eps)(double x, double* dist_pars);
        double (*aux_pot_eps)(double x, double* dist_pars);
        
        arma::vec   (*aux_cdf_eps_vec)(arma::vec x, double* dist_pars);
        arma::vec (*aux_quant_eps_vec)(arma::vec x, double* dist_pars);
       
        arma::vec (*aux_pdf_eps_vec)(arma::vec x, double* dist_pars);
        arma::vec (*aux_pot_eps_vec)(arma::vec x, double* dist_pars);
        
        // input objects
        arma::mat mu;
        arma::mat U;
        
        // equilibrium objects
        arma::mat mu_sol;
        arma::mat U_sol;
        
        // member functions
        void build(arma::mat zeta_b, bool outsideOption_b);

        void build_beta(arma::mat zeta_b, double alpha, double beta);
        
        double G(arma::vec n);
        double Gx(arma::vec& mu_x, int x);
        
        double Gstar(arma::vec n);
        double Gstar(arma::mat& U_inp, arma::vec n);
        double Gstarx(arma::vec& U_x, double n_x, int x);
        double Gstarx(arma::vec& U_x, arma::vec mu_x_inp, int x);
        static double Gstarx(arma::vec& U_x, arma::vec mu_x_inp, arma::mat zeta, 
                             arma::mat aux_DinvPsigma, arma::mat aux_Psigma, 
                             arma::mat aux_Influence_lhs, arma::mat aux_Influence_rhs,
                             arma::vec (*pot_eps_vec)(arma::vec pot_inp, double* dist_pars),
                             arma::vec (*quantile_eps_vec)(arma::vec quant_inp, double* dist_pars),
                             double* dist_pars, int nbY, int x);

        double Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_inp, arma::mat& mu_inp);
        double Gbarx(arma::vec Ubarx, arma::vec mubarx, arma::mat& Ux_inp, arma::mat& mu_x_inp, int x);
        
        void D2Gstar (arma::mat& hess, arma::vec n, bool x_first);
        void dtheta_NablaGstar (arma::mat& ret, arma::vec n, arma::mat* dtheta, bool x_first);
        
        void simul(empirical &ret, int nbDraws, int seed);
        
        double cdf (double x);
        arma::vec cdf (arma::vec x);
        double pdf (double x);
        arma::vec pdf (arma::vec x);
        double quantile (double x);
        arma::vec quantile (arma::vec x);
        double pot (double x);
        arma::vec pot (arma::vec x);
        
    private:
        static double Gbar_opt_objfn(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data);
        static double Gbar_opt_constr(const std::vector<double> &x_inp, std::vector<double> &grad, void *constr_data);
};
/*
void RSC::build(arma::mat zeta_b, bool outsideOption_b, 
                double   (*aux_cdf_eps_b)(double x, double* dist_pars),
                double (*aux_quant_eps_b)(double x, double* dist_pars),
                double   (*aux_pdf_eps_b)(double x, double* dist_pars),
                double   (*aux_pot_eps_b)(double x, double* dist_pars))
*/
void RSC::build(arma::mat zeta_b, bool outsideOption_b)
{
    if (!outsideOption_b) {
        return;
    }
    //
    int i,j;
    //
    nbX = zeta_b.n_rows;
    nbY = zeta_b.n_cols - 1;
    nbParams = zeta_b.n_elem;
    
    zeta = zeta_b;
    //
    aux_ord = arma::zeros(nbX,nbY+1);
    
    arma::mat D =  arma::eye(nbY+1,nbY+1) - arma::join_cols(arma::zeros(1,nbY+1), arma::join_rows(arma::eye(nbY,nbY),arma::zeros(nbY,1)));
    arma::mat D_inv = arma::inv(D);
    
    arma::mat neg_ones(nbY,1);
    neg_ones.fill(-1);
    arma::mat N_temp = arma::join_rows(arma::eye(nbY,nbY),neg_ones);
    //
    arma::uvec ordx_temp;
    arma::mat Psigmax(nbY+1,nbY+1);
    
    aux_Influence_lhs.set_size(nbY,nbY+1,nbX);
    aux_Influence_lhs.zeros();
    aux_Influence_rhs.set_size(nbY+1,nbY+1,nbX);
    aux_Influence_rhs.zeros();
    aux_Psigma.set_size(nbY+1,nbY+1,nbX);
    aux_Psigma.zeros();
    aux_DinvPsigma.set_size(nbY+1,nbY+1,nbX);
    aux_DinvPsigma.zeros();
    //
    for (i=0; i<nbX; i++) {
        ordx_temp = arma::sort_index(zeta_b.row(i));
        aux_ord.row(i) = arma::conv_to< arma::rowvec >::from(ordx_temp);
        
        Psigmax.zeros();
        for (j=0; j<nbY+1; j++) {
            Psigmax(j,ordx_temp(j)) = 1;
        }
        
        aux_Influence_lhs.slice(i) = N_temp * (Psigmax.t() * D_inv);
        aux_Influence_rhs.slice(i) = D * Psigmax;
        
        aux_Psigma.slice(i) = arma::eye(nbY+1,nbY+1) * Psigmax;
        aux_DinvPsigma.slice(i) = D_inv * Psigmax;
    }
    //
    outsideOption = true;
}

// epsilon is a beta(alpha,beta) distribution
void RSC::build_beta(arma::mat zeta_b, double alpha, double beta)
{
    dist_pars = new double[2];
    dist_pars[0] = alpha;
    dist_pars[1] = beta;
    
    aux_cdf_eps   = pbeta;
    aux_quant_eps = qbeta;
    aux_pdf_eps   = dbeta;
    aux_pot_eps   = iqbeta;
    
    aux_cdf_eps_vec   = pbeta;
    aux_quant_eps_vec = qbeta;
    aux_pdf_eps_vec   = dbeta;
    aux_pot_eps_vec   = iqbeta;
    //
    RSC::build(zeta_b,true);
}

double RSC::G(arma::vec n)
{   
    int i;
    double val=0.0, val_x_temp;
    
    mu_sol.set_size(nbX,nbY);
    arma::vec mu_x;
    //
    //#pragma omp parallel for private(mu_x,val_x_temp)
        for (i=0; i<nbX; i++) {
            val_x_temp = Gx(mu_x,i);
            //
            val += n(i)*val_x_temp;
            mu_sol.row(i) = arma::trans(n(i)*mu_x);
        }
    //
    return val;
}

double RSC::Gx(arma::vec& mu_x, int x)
{
    int nbAlt = nbY + 1;
    int i,j,y,z;
    
    double val_x=0, E_eps_temp=0, E_eps_temp_next=0, cumul_temp=0; 
    double run_max=0, run_min=0, run_temp=0;
    double mu_x_tilde_y, e_y;
    
    arma::vec U_x = U.row(x).t();
    
    arma::vec mu_x_tilde = arma::zeros(nbAlt,1);
    arma::vec U_x_tilde = arma::join_cols(arma::vectorise(U_x),arma::zeros(1,1));
    //
    for (i=0; i<nbAlt; i++) {
        y = aux_ord(x,i);
        run_max = quantile(0);
        //
        j = 0;
        
        while (j < i) {
            z = aux_ord(x,j);
            
            if (zeta(x,z) != zeta(x,y)) {
                run_temp = (U_x_tilde(y) - U_x_tilde(z)) / (zeta(x,z) - zeta(x,y)); 
                run_max = std::max(run_max,run_temp);
            } else {
                run_max = INFINITY;
            }
            
            j++;
        }
        //
        run_min = quantile(1);
        //
        j = nbAlt-1;
        
        while (j > i) {
            z = aux_ord(x,j);
            
            if (zeta(x,z) != zeta(x,y)) {
                run_temp = (U_x_tilde(y) - U_x_tilde(z)) / (zeta(x,z) - zeta(x,y)); 
                run_min = std::min(run_min,run_temp);
            }
            
            j--;
        }
        //
        if (run_min > run_max) {
            mu_x_tilde_y = std::max(cdf(run_min) - cdf(run_max),0.0);
        } else {
            mu_x_tilde_y = 0;
        }
        
        mu_x_tilde(y) = mu_x_tilde_y;
        //
        if (mu_x_tilde_y > 0) {
            cumul_temp += mu_x_tilde_y;
            e_y = quantile(cumul_temp);
            
            E_eps_temp_next = e_y * cumul_temp - pot(cumul_temp);
            
            val_x += mu_x_tilde_y*U_x_tilde(y) + zeta(x,y)*(E_eps_temp_next - E_eps_temp);
        }
    }
    //
    mu_x = mu_x_tilde.rows(0,nbAlt-2);
    //
    return val_x;
}

double RSC::Gstar(arma::vec n)
{   
    int i;
    double val=0.0, val_x_temp;
    
    U_sol.set_size(nbX,nbY);
    arma::vec U_x_temp;
    //
    for (i=0; i<nbX; i++) {
        //val_x_temp = Gstarx((mu.row(i).t())/n(i),U_x_temp,i);
        val_x_temp = Gstarx(U_x_temp,n(i),i);
        //
        val += n(i)*val_x_temp;
        U_sol.row(i) = arma::trans(U_x_temp);
    }
    //
    return val;
}

double RSC::Gstar(arma::mat& U_inp, arma::vec n)
{   
    int i;
    double val=0.0, val_x_temp;
    
    U_inp.set_size(nbX,nbY);
    arma::vec U_x_temp;
    //
    for (i=0; i<nbX; i++) {
        //val_x_temp = Gstarx((mu.row(i).t())/n(i),U_x_temp,i);
        val_x_temp = Gstarx(U_x_temp,n(i),i);
        //
        val += n(i)*val_x_temp;
        U_inp.row(i) = arma::trans(U_x_temp);
    }
    //
    return val;
}

double RSC::Gstarx(arma::vec& U_x, double n_x, int x)
{
    double val_x = 0;
    arma::vec mu_x = (mu_sol.row(x).t())/n_x; // we divide by n(x)
    
    arma::vec ts_temp(mu_x.n_elem+1);
    ts_temp.rows(0,mu_x.n_elem-1) = mu_x;
    ts_temp(mu_x.n_elem) = 1-arma::accu(mu_x);
    
    arma::vec ts_full = aux_DinvPsigma.slice(x) * ts_temp;
    arma::vec ts = ts_full.rows(0,nbY-1);
    //
    arma::vec pots_inp = arma::join_cols(arma::zeros(1,1),ts_full);
    arma::vec pots = pot(pots_inp);
    
    arma::vec diff_pots = pots.rows(1,nbY+1) - pots.rows(0,nbY);
    //
    val_x = - arma::accu( (aux_Psigma.slice(x)*zeta.row(x).t()) % diff_pots );
    
    arma::mat e_mat = arma::diagmat( arma::join_cols(arma::zeros(1,1),quantile(ts)) );
    U_x = - aux_Influence_lhs.slice(x) * e_mat * aux_Influence_rhs.slice(x) * zeta.row(x).t();
    //
    return val_x;
}

double RSC::Gstarx(arma::vec& U_x, arma::vec mu_x_inp, int x)
{
    double val_x = 0;
    
    arma::vec ts_temp(mu_x_inp.n_elem+1);
    ts_temp.rows(0,mu_x_inp.n_elem-1) = mu_x_inp;
    ts_temp(mu_x_inp.n_elem) = 1-arma::accu(mu_x_inp);
    
    arma::vec ts_full = aux_DinvPsigma.slice(x) * ts_temp;
    arma::vec ts = ts_full.rows(0,nbY-1);
    //
    arma::vec pots_inp = arma::join_cols(arma::zeros(1,1),ts_full);
    arma::vec pots = pot(pots_inp);
    
    arma::vec diff_pots = pots.rows(1,nbY+1) - pots.rows(0,nbY);
    //
    val_x = - arma::accu( (aux_Psigma.slice(x)*zeta.row(x).t()) % diff_pots );
    
    arma::mat e_mat = arma::diagmat( arma::join_cols(arma::zeros(1,1),quantile(ts)) );
    U_x = - aux_Influence_lhs.slice(x) * e_mat * aux_Influence_rhs.slice(x) * zeta.row(x).t();
    //
    return val_x;
}

double RSC::Gstarx(arma::vec& U_x, arma::vec mu_x_inp, arma::mat zeta, 
                   arma::mat aux_DinvPsigma, arma::mat aux_Psigma, 
                   arma::mat aux_Influence_lhs, arma::mat aux_Influence_rhs,
                   arma::vec (*pot_eps_vec)(arma::vec pot_inp, double* dist_pars),
                   arma::vec (*quantile_eps_vec)(arma::vec quant_inp, double* dist_pars),
                   double* dist_pars, int nbY, int x)
{
    double val_x = 0;
    
    arma::vec ts_temp(mu_x_inp.n_elem+1);
    ts_temp.rows(0,mu_x_inp.n_elem-1) = mu_x_inp;
    ts_temp(mu_x_inp.n_elem) = 1-arma::accu(mu_x_inp);
    
    arma::vec ts_full = aux_DinvPsigma * ts_temp;
    arma::vec ts = ts_full.rows(0,nbY-1);
    //
    arma::vec pots_inp = arma::join_cols(arma::zeros(1,1),ts_full);
    arma::vec pots = (*pot_eps_vec)(pots_inp, dist_pars);
    
    arma::vec diff_pots = pots.rows(1,nbY+1) - pots.rows(0,nbY);
    //
    val_x = - arma::accu( (aux_Psigma*zeta.row(x).t()) % diff_pots );
    
    arma::mat e_mat = arma::diagmat( arma::join_cols(arma::zeros(1,1),(*quantile_eps_vec)(ts,dist_pars)) );
    U_x = - aux_Influence_lhs * e_mat * aux_Influence_rhs * zeta.row(x).t();
    //
    return val_x;
}

double RSC::Gbar(arma::mat Ubar, arma::mat mubar, arma::vec n, arma::mat& U_inp, arma::mat& mu_inp)
{   
    int i;
    double val=0.0, val_temp;

    U_inp.set_size(nbX,nbY);
    mu_inp.set_size(nbX,nbY);
    arma::mat Ux_temp, mu_x_temp;
    //
    for (i=0; i<nbX; i++) {
        val_temp = Gbarx(Ubar.row(i).t(),(mubar.row(i).t())/n(i),Ux_temp,mu_x_temp,i);
        
        val += n(i)*val_temp;
        U_inp.row(i) = arma::trans(Ux_temp);
        mu_inp.row(i) = arma::trans(n(i)*mu_x_temp);
    }
    //
    return val;
}

double RSC::Gbarx(arma::vec Ubarx, arma::vec mubarx, arma::mat& Ux_inp, arma::mat& mu_x_inp, int x)
{
    if (!outsideOption) {
        printf("Gbarx not implemented yet when outsideOption==false");
        return 0.0;
    }
    //
    arma::vec lb = arma::zeros(nbY,1);
    arma::vec ub = mubarx;

    std::vector<double> sol_vec = arma::conv_to< std::vector<double> >::from(mubarx/2.0);
    double obj_val = 0;
    double ret = 0;
    //
    // opt data
    trame_nlopt_opt_data opt_data;

    opt_data.x = x;
    opt_data.nbY = nbY;
    opt_data.Ubar_x = Ubarx;
    opt_data.zeta = zeta;

    opt_data.aux_DinvPsigma = aux_DinvPsigma.slice(x);
    opt_data.aux_Psigma = aux_Psigma.slice(x);
    opt_data.aux_Influence_lhs = aux_Influence_lhs.slice(x);
    opt_data.aux_Influence_rhs = aux_Influence_rhs.slice(x);

    opt_data.pot_eps_vec = aux_pot_eps_vec;
    opt_data.quantile_eps_vec = aux_quant_eps_vec;

    opt_data.dist_pars = dist_pars;

    //
    // constr data
    trame_nlopt_constr_data constr_data;
    constr_data.nbY = nbY;
    //
    bool success = generic_nlopt(nbY,sol_vec,obj_val,lb.memptr(),ub.memptr(),RSC::Gbar_opt_objfn,RSC::Gbar_opt_constr,opt_data,constr_data);
    //
    arma::vec Ux_temp;
    arma::vec sol_temp = arma::conv_to<arma::vec>::from(sol_vec);
    
    if (success) {
        mu_x_inp = sol_temp;
        Gstarx(Ux_temp,sol_temp,x);
        Ux_inp = Ux_temp;
        ret = -obj_val;
    }
    return ret;
}

void RSC::D2Gstar(arma::mat& hess, arma::vec n, bool x_first)
{
    int i,j;
    
    hess.set_size(nbX*nbY,nbX*nbY);
    hess.zeros();
    
    arma::vec mu_x_0 = n - arma::sum(mu,1);
    
    arma::umat mat_inds(nbY,nbX); // indices to place the results (complicated)
    for (i=0; i<nbX; i++) {
        for (j=0; j<nbY; j++) {
            if (x_first) {
                mat_inds(j,i) = i + j*nbX;
            } else {
                mat_inds(j,i) = i*nbY + j;
            }
        }
    }
    //
    arma::mat C, d_mu_e_temp, d_mu_e;
    arma::vec ts_temp(nbY+1), ts_full, erestr_temp, erestr, erestr_pdf;
    
    for (i=0; i<nbX; i++) {
        C = - arma::trans( arma::repmat(aux_Influence_rhs.slice(i) * zeta.row(i).t(),1,nbY) % aux_Influence_lhs.slice(i).t() );
        
        ts_temp.rows(0,nbY-1) = mu.row(i).t();
        ts_temp(nbY) = mu_x_0(i);
        
        ts_full = aux_DinvPsigma.slice(i) * ts_temp / n(i);
        //
        erestr_temp = quantile(ts_full);
        erestr = erestr_temp.rows(0,nbY-1);
        erestr_pdf = pdf(erestr);
        //
        d_mu_e_temp = arma::join_cols(arma::zeros(1,nbY),arma::diagmat(1 / erestr_pdf)) * aux_DinvPsigma.slice(i).rows(0,nbY-1);
        d_mu_e = d_mu_e_temp.cols(0,nbY-1) - arma::repmat(d_mu_e_temp.col(nbY),1,nbY);
        //
        hess(mat_inds.col(i),mat_inds.col(i)) = C * d_mu_e / n(i);
    }
}

void RSC::dtheta_NablaGstar(arma::mat& ret, arma::vec n, arma::mat* dtheta, bool x_first)
{
    int i,j;
    arma::mat dtheta_mat;
    
    if (dtheta==NULL) {
        dtheta_mat = arma::eye(nbParams,nbParams);
    } else {
        dtheta_mat = *dtheta;
    }
    
    int nbDirs = std::floor(dtheta_mat.n_elem / (nbX*nbX*(nbY+1)));
    
    ret.set_size(nbX*nbY,nbX*nbDirs);
    ret.zeros();
    
    arma::umat mat_inds_1(nbDirs,nbX); // indices to place the results (complicated)
    arma::umat mat_inds_2(nbDirs,nbX);
    arma::umat mat_inds_3(nbDirs,nbX);
    
    for (i=0; i<nbX; i++) {
        for (j=0; j<nbDirs; j++) {
            if (x_first) {
                mat_inds_1(j,i) = i + j*nbX;
                mat_inds_2(j,i) = i + j*nbX;
                mat_inds_3(j,i) = i + j*nbX;
            } else {
                mat_inds_1(j,i) = i + j*nbX;
                mat_inds_2(j,i) = i + j*nbX;
                mat_inds_3(j,i) = i*nbY + j;
            }
        }
    }
    
    if (nbDirs > nbY) {
        mat_inds_3.shed_row(nbDirs-1);
    }
    //
    arma::mat e_mat;
    arma::vec ts_temp(nbY+1), ts_full, ts;
    
    for (i=0; i<nbX; i++) {
        ts_temp.rows(0,nbY-1) = mu.row(i).t();
        ts_temp(nbY) = n(i) - arma::accu(mu.row(i));
        
        arma::vec ts_full = aux_DinvPsigma.slice(i) * ts_temp / n(i);
        arma::vec ts = ts_full.rows(0,nbY-1);

        e_mat = arma::diagmat( arma::join_cols(arma::zeros(1,1),quantile(ts)) );
        
        ret(mat_inds_3.col(i),mat_inds_1.col(i)) = - aux_Influence_lhs.slice(i) * e_mat * aux_Influence_rhs.slice(i) * dtheta_mat(mat_inds_2.col(i),mat_inds_2.col(i));
    }
}

void RSC::simul(empirical &ret, int nbDraws, int seed_val)
{
    int i;
    arma::arma_rng::set_seed(seed_val);
    //
    arma::cube atoms(nbDraws,nbY+1,nbX);
    
    for (i=0; i<nbX; i++) {
        atoms.slice(i) = quantile(arma::randu(nbDraws,1)) * zeta.row(i);
    }
    //
    ret.nbX = nbX;
    ret.nbY = nbY;
    ret.nbParams = atoms.n_elem;
    ret.atoms = atoms;
    ret.aux_nbDraws = nbDraws;
    ret.xHomogenous = false;
    ret.outsideOption = outsideOption;
    if (outsideOption) {
        ret.nbOptions = nbY + 1;
    } else {
        ret.nbOptions = nbY;
    }
    //
    arma::arma_rng::set_seed_random(); // need to reset the seed
}

/*
 * optimization-related functions
 */

double RSC::Gbar_opt_objfn(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data)
{
    trame_nlopt_opt_data *d = reinterpret_cast<trame_nlopt_opt_data*>(opt_data);

    int x = d->x;
    arma::mat Ubar_x = d->Ubar_x;
    arma::mat zeta = d->zeta;
    arma::mat aux_DinvPsigma = d->aux_DinvPsigma;
    arma::mat aux_Psigma = d->aux_Psigma;
    arma::mat aux_Influence_lhs = d->aux_Influence_lhs;
    arma::mat aux_Influence_rhs = d->aux_Influence_rhs;

    double* dist_pars = d->dist_pars;

    int nbY = d->nbY;

    arma::vec U_x_temp;
    arma::vec mu_x_inp = arma::conv_to<arma::vec>::from(x_inp);

	double val_x = Gstarx(U_x_temp, mu_x_inp, zeta,
                          aux_DinvPsigma,aux_Psigma,
                          aux_Influence_lhs,aux_Influence_rhs,
                          d->pot_eps_vec,d->quantile_eps_vec,
                          dist_pars,nbY,x);
	
    double ret = val_x - arma::accu(mu_x_inp % Ubar_x);
    //
    if (!grad.empty()) {
        grad = arma::conv_to< std::vector<double> >::from(U_x_temp - Ubar_x);
    }
    //
    return ret;
}

double RSC::Gbar_opt_constr(const std::vector<double> &x_inp, std::vector<double> &grad, void *constr_data)
{
    trame_nlopt_constr_data *d = reinterpret_cast<trame_nlopt_constr_data*>(constr_data);

    int nbY = d->nbY;

    arma::vec mu_x_inp = arma::conv_to<arma::vec>::from(x_inp);

    double ret = arma::accu(mu_x_inp) - 1;
    //
    if (!grad.empty()) {
        grad = arma::conv_to< std::vector<double> >::from(arma::ones<arma::vec>(nbY,1));
    }
    //
    return ret;
}

/*
 * Distribution-related functions
 */

double RSC::cdf(double x)
{
    double res;
    
    res = (*aux_cdf_eps)(x, dist_pars);
    
    return res;
}

arma::vec RSC::cdf(arma::vec x)
{
    arma::vec res;
    
    res = (*aux_cdf_eps_vec)(x, dist_pars);
    
    return res;
}

double RSC::pdf(double x)
{
    double res;
    
    res = (*aux_pdf_eps)(x, dist_pars);
    
    return res;
}

arma::vec RSC::pdf(arma::vec x)
{
    arma::vec res;
    
    res = (*aux_pdf_eps_vec)(x, dist_pars);
    
    return res;
}

double RSC::quantile(double x)
{
    double res;
    
    res = (*aux_quant_eps)(x, dist_pars);
    
    return res;
}

arma::vec RSC::quantile(arma::vec x)
{
    arma::vec res;
    
    res = (*aux_quant_eps_vec)(x, dist_pars);
    
    return res;
}

double RSC::pot(double x)
{
    double res;
    
    res = (*aux_pot_eps)(x, dist_pars);
    
    return res;
}

arma::vec RSC::pot(arma::vec x)
{
    arma::vec res;
    
    res = (*aux_pot_eps_vec)(x, dist_pars);
    
    return res;
}

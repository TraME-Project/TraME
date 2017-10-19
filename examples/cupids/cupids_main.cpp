
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
 * Cupids example
 *
 * Keith O'Hara
 * 07/17/2017
 *
 * This version:
 * 08/21/2017
 */

// g++-mp-7 -O2 -Wall -std=c++11 -I/../../include cupids_1.cpp -o cupids_1.test -framework Accelerate -L/../.. -ltrame
// g++ -O2 -Wall -std=c++11 -I/usr/local/include/trame cupids_1.cpp -o cupids_1.test -L/usr/local/lib -ltrame -lopenblas

#include "trame.hpp"
#include "headers/misc.hpp"
#include "headers/build_model.hpp"

int main()
{
    arma::mat nSingles, marr, nAvail;
    nSingles.load("data/n_singles.txt",arma::auto_detect);
    marr.load("data/marr.txt",arma::auto_detect);
    nAvail.load("data/n_avail.txt",arma::auto_detect);
    //
    int spec_id = 10;
    // bool short_version = true;

    const int n_draws = 200;
    const int sim_seed = 777;
    //
    int nb_categ = 25;
    arma::mat mu_hat_x0 = nSingles(arma::span(0,nb_categ-1),0);
    arma::mat mu_hat_0y = nSingles(arma::span(0,nb_categ-1),1);

    arma::mat mu_hat_xy = marr(arma::span(0,nb_categ-1),arma::span(0,nb_categ-1));

    arma::mat n = nAvail(arma::span(0,nb_categ-1),0);
    arma::mat m = nAvail(arma::span(0,nb_categ-1),1);
    //
    int nb_couples = static_cast<int>(arma::accu(mu_hat_xy));

    mu_hat_x0 /= nb_couples;
    mu_hat_0y /= nb_couples;
    mu_hat_xy /= nb_couples;

    n /= nb_couples;
    m /= nb_couples;

    int nbX = n.n_elem;
    int nbY = m.n_elem;
    //
    if (nbX != nbY) {
        printf("error: nbX != nbY\n");
        return 1;
    }

    const int optim_method = 2; 

    //
    // build the phi_xyk cube

    printf("*** Estimations for specification #%u ***\n",spec_id);
    arma::cube phi_xyk = build_tu_model_inp(spec_id,nb_categ);
    
    const int n_theta = phi_xyk.n_slices;
    //
    // none model

    printf("none model:\n");
    trame::model< trame::dse<trame::arums::none, trame::arums::none, trame::transfers::tu> > cupids_none_model;

    cupids_none_model.build(phi_xyk,n,m);

    arma::mat theta_hat_1;

    cupids_none_model.mme(mu_hat_xy,theta_hat_1,nullptr);
    
    arma::cout << "none_model theta_hat:\n" << theta_hat_1 << arma::endl;
    
    //
    // logit model

    printf("logit model:\n");
    trame::model< trame::dse<trame::arums::logit, trame::arums::logit, trame::transfers::tu> > cupids_logit_model;

    cupids_logit_model.build(phi_xyk,n,m);

    arma::mat theta_hat_2;
    // arma::mat initial_theta_2 = arma::zeros(n_theta,1);

    cupids_logit_model.mme(mu_hat_xy,theta_hat_2,nullptr,&optim_method);
    // cupids_logit_model.mle(mu_hat_xy,theta_hat_2,nullptr);
    
    arma::cout << "logit_model theta_hat:\n" << theta_hat_2 << arma::endl;

    //
    // logit model via simulation

    printf("logit_sim model:\n");

    trame::arums::logit logit_1(nbX,nbY), logit_2(nbY,nbX);

    trame::arums::empirical logit_sim_1(nbX,nbY), logit_sim_2(nbY,nbX);
    logit_1.simul(logit_sim_1, n_draws, sim_seed);
    logit_2.simul(logit_sim_2, n_draws, sim_seed);

    trame::model< trame::dse<trame::arums::empirical, trame::arums::empirical, trame::transfers::tu> > cupids_logit_sim_model;

    cupids_logit_sim_model.build(phi_xyk,n,m);
    cupids_logit_sim_model.market_obj.arums_G = logit_sim_1;
    cupids_logit_sim_model.market_obj.arums_H = logit_sim_2;

    arma::mat theta_hat_3;

    cupids_logit_sim_model.mme(mu_hat_xy,theta_hat_3,nullptr);
    
    arma::cout << "logit_sim_model theta_hat:\n" << theta_hat_3 << arma::endl;

    return 0;
}

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
 * Derived classes to provide wrappers to the TraME library
 */

class empirical_R : public trame::empirical
{
    public:
        SEXP G_R(arma::vec n);
        SEXP G_R(arma::vec n, arma::mat U_inp);
        SEXP Gx_R(arma::mat U_x_inp, int x);
        SEXP Gstar_R(arma::vec n);
        SEXP Gstar_R(arma::vec n, arma::mat mu_inp);
        SEXP Gstarx_R(arma::mat mu_x_inp, int x);
        SEXP Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);
};

class logit_R : public trame::logit
{
    public:
        SEXP G_R(arma::vec n);
        SEXP G_R(arma::vec n, arma::mat U_inp);
        SEXP Gstar_R(arma::vec n);
        SEXP Gstar_R(arma::vec n, arma::mat mu_inp);
        SEXP Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);

        empirical_R simul_R(int nbDraws);
};

class none_R : public trame::none
{
    public:
        SEXP G_R(arma::vec n);
        SEXP G_R(arma::vec n, arma::mat U_inp);
        SEXP Gx_R(arma::mat U_x_inp, int x);
        SEXP Gstar_R(arma::vec n);
        SEXP Gstar_R(arma::vec n, arma::mat mu_inp);
        SEXP Gstarx_R(arma::mat mu_x_inp, int x);
        SEXP Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);

        empirical_R simul_R(int nbDraws);
};

class probit_R : public trame::probit
{
    public:
        empirical_R simul_R(int nbDraws);
};

class rsc_R : public trame::rsc
{
    public:
        SEXP G_R(arma::vec n);
        SEXP G_R(arma::vec n, arma::mat U_inp);
        SEXP Gx_R(arma::mat U_x_inp, int x);
        SEXP Gstar_R(arma::vec n);
        SEXP Gstar_R(arma::vec n, arma::mat mu_inp);
        SEXP Gstarx_R(arma::mat mu_x_inp, int x);
        SEXP Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);

        empirical_R simul_R(int nbDraws);
};

class rusc_R : public trame::rusc
{
    public:
        SEXP G_R(arma::vec n);
        SEXP G_R(arma::vec n, arma::mat U_inp);
        SEXP Gx_R(arma::mat U_x_inp, int x);
        SEXP Gstar_R(arma::vec n);
        SEXP Gstar_R(arma::vec n, arma::mat mu_inp);
        SEXP Gstarx_R(arma::mat mu_x_inp, int x);
        SEXP Gbar_R(arma::mat U_bar, arma::mat mu_bar, arma::vec n);

        empirical_R simul_R(int nbDraws);
};

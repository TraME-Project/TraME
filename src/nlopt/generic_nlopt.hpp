#ifndef _generic_nlopt_HPP
#define _generic_nlopt_HPP

#include <RcppArmadillo.h>
#include <nlopt.hpp>
#include "../aux/trame_structs.hpp"

bool generic_nlopt(int n_pars, std::vector<double>& io_val, double& opt_val, double* lb, double* ub,
                   double (*opt_objfn)(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data),
                   trame_nlopt_opt_data opt_data);
bool generic_nlopt(int n_pars, std::vector<double>& io_val, double& opt_val, double* lb, double* ub,
                   double (*opt_objfn)(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data),
                   double (*opt_constr)(const std::vector<double> &x_inp, std::vector<double> &grad, void *constr_data),
                   trame_nlopt_opt_data opt_data,
                   trame_nlopt_constr_data constr_data);

#endif

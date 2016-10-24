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
 * Generic NLopt
 *
 * Keith O'Hara
 * 08/08/2016
 */

#include "trame.hpp"

bool trame::generic_nlopt(int n_pars, std::vector<double>& io_val, double& opt_val, double* lb, double* ub,
                          double (*opt_objfn)(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data),
                          trame_nlopt_opt_data opt_data)
{
    bool success = false;

    nlopt::opt opt_trame(nlopt::LD_MMA, n_pars);

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

bool trame::generic_nlopt(int n_pars, std::vector<double>& io_val, double& opt_val, double* lb, double* ub,
                          double (*opt_objfn)(const std::vector<double> &x_inp, std::vector<double> &grad, void *opt_data),
                          double (*opt_constr)(const std::vector<double> &x_inp, std::vector<double> &grad, void *constr_data),
                          trame_nlopt_opt_data opt_data,
                          trame_nlopt_constr_data constr_data)
{
    bool success = false;

    nlopt::opt opt_trame(nlopt::LD_MMA, n_pars);

    if (lb) {
        opt_trame.set_lower_bounds(*lb);
    }
    if (ub) {
        opt_trame.set_upper_bounds(*ub);
    }

    opt_trame.set_min_objective(*opt_objfn, &opt_data);
    opt_trame.add_inequality_constraint(*opt_constr, &constr_data, 1e-8);
    
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

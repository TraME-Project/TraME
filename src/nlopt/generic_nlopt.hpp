/*
 * Generic NLopt
 *
 * Keith O'Hara
 * 05/08/2016
 *
 * clang -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include generic_lp.c -c -o generic_lp.o
 * gcc-mp-5 -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include generic_lp.c -c -o generic_lp.o
 */

#include <nlopt.hpp>
//#include "trame_opt_fns.hpp"

bool generic_nlopt(int n_pars, std::vector<double>& io_val, double& opt_val, double* lb, double* ub,
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
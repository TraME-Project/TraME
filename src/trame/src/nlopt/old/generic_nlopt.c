/*
 * Generic NLopt function
 *
 * Keith O'Hara
 * 05/08/2016
 *
 * clang -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include generic_lp.c -c -o generic_lp.o
 * gcc-mp-5 -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include generic_lp.c -c -o generic_lp.o
 */
 
#include <string.h>

#include "generic_nlopt.h"

//
// Dense setup; should be used with Armadillo memptr based passing 
int generic_nlopt_C(int n_pars, double* io_val, double* lb, double* ub,
                    double* c_data)
{
    bool success = false;
    nlopt_opt opt_trame; 

    opt = nlopt_create(NLOPT_LD_MMA, n); /* algorithm and dimensionality */

    if (lb) {
        nlopt_set_lower_bounds(opt_trame, lb);
    }

    if (ub) {
        nlopt_set_upper_bounds(opt_trame, ub);
    }

    if (c_data) {
        nlopt_add_inequality_constraint(opt, trame_gbar_cons_func, c_data, 1e-8);
    }

    nlopt_set_min_objective(opt_trame, trame_gbar_obj_func, NULL);
    nlopt_set_xtol_rel(opt_trame, 1e-7);

    double minf;
    int res_nlopt = nlopt_optimize(opt, io_val, &minf);

    /* Free model */
    
    nlopt_destroy(opt_trame);

    return success;
}
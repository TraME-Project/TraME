
#ifndef TRAMELIB_INCL
#define TRAMELIB_INCL

#include <cmath> // for abs value
#include <limits>

#include "misc/TRAME_OPTIONS.hpp"

namespace trame
{
    #include "aux/trame_structs.hpp"
    #include "aux/trame_aux.hpp"
    #include "aux/trame_stats.hpp"

    #include "aux/inv_pwa.hpp"
    #include "aux/zeroin.hpp"

    #include "lp/generic_lp.hpp"
    #include "nlopt/generic_nlopt.hpp"

    #include "optim/optim.hpp"

    // heterogeneity
    #include "arums/arums_empirical.hpp"
    #include "arums/arums_logit.hpp"
    #include "arums/arums_none.hpp"
    #include "arums/arums_probit.hpp"
    #include "arums/arums_rsc.hpp"
    #include "arums/arums_rusc.hpp"

    // markets
    #include "markets/mmf.hpp"
    #include "markets/transfers.hpp"
    #include "markets/mfe.hpp"
    #include "markets/dse.hpp"

    // solvers
    #include "solvers/aux_solvers.hpp"
    #include "solvers/ipfp.hpp"
    #include "solvers/oap_lp.hpp"
    #include "solvers/darum.hpp"
    #include "solvers/eap_nash.hpp"
    #include "solvers/jacobi.hpp"
    #include "solvers/cupids_lp.hpp"
    #include "solvers/max_welfare.hpp"

    // models
    #include "models/affinity.hpp"
    #include "models/model.hpp"
}

#endif

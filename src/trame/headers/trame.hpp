
#ifndef TRAMELIB_INCLUDES
#define TRAMELIB_INCLUDES

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

    // heterogeneity
    #include "arums/arums_empirical.hpp"
    #include "arums/arums_logit.hpp"
    #include "arums/arums_none.hpp"
    #include "arums/arums_probit.hpp"
    #include "arums/arums_rsc.hpp"
    #include "arums/arums_rusc.hpp"

    // markets
    #include "markets/mmf.hpp"
}

#endif

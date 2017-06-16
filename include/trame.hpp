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

#ifndef TRAMELIB_INCL
#define TRAMELIB_INCL

#include <cmath> // for abs value
#include <limits>

#include "misc/TRAME_OPTIONS.hpp"
#include "optim/optim.hpp"

namespace trame
{
    #include "ancillary/trame_aux.hpp"
    #include "ancillary/trame_stats.hpp"
    #include "ancillary/logit_transform.hpp"

    #include "ancillary/inv_pwa.hpp"
    #include "ancillary/zeroin.hpp"

    #include "lp/generic_lp.hpp"

    // heterogeneity
    namespace arums {
        #include "arums/arums_empirical.hpp"
        #include "arums/arums_logit.hpp"
        #include "arums/arums_none.hpp"
        #include "arums/arums_probit.hpp"
        #include "arums/arums_rsc.hpp"
        #include "arums/arums_rusc.hpp"
    }

    // matching functions
    namespace mmfs {
        #include "mmfs/mmfs_ces.hpp"
        #include "mmfs/mmfs_cd.hpp"
        #include "mmfs/mmfs_min.hpp"
        #include "mmfs/mmfs_geo.hpp"
    }

    // transfers
    namespace transfers {
        #include "transfers/transfers_etu.hpp"
        #include "transfers/transfers_ltu.hpp"
        #include "transfers/transfers_ntu.hpp"
        #include "transfers/transfers_tu.hpp"
    }

    // markets
    #include "markets/dse.hpp"
    #include "markets/mfe.hpp"

    // solvers
    #include "solvers/aux_solvers.hpp"
    
    #include "solvers/arc_newton.hpp"
    #include "solvers/cupids_lp.hpp"
    #include "solvers/darum.hpp"
    #include "solvers/eap_nash.hpp"
    #include "solvers/ipfp.hpp"
    #include "solvers/jacobi.hpp"
    #include "solvers/max_welfare.hpp"
    #include "solvers/nodal_newton.hpp"
    #include "solvers/oap_lp.hpp"

    #include "solvers/equil_solve.hpp"

    // models
    #include "models/affinity.hpp"
    #include "models/model.hpp"
}

#endif

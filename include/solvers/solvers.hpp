/*################################################################################
  ##
  ##   Copyright (C) 2015-2018 the TraME Team:
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

#ifndef TRAMELIB_SLVRS_INCL
#define TRAMELIB_SLVRS_INCL

namespace trame
{
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
}

#endif

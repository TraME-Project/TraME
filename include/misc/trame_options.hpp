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

#pragma once

#ifdef USE_RCPP_ARMADILLO
    #include <RcppArmadillo.h>
#else
    #ifndef ARMA_DONT_USE_WRAPPER
        #define ARMA_DONT_USE_WRAPPER
    #endif
    #include "armadillo"
#endif

#if defined(_OPENMP) && !defined(TRAME_DONT_USE_OPENMP)
    #undef TRAME_USE_OPENMP
    #define TRAME_USE_OPENMP
#endif

#if !defined(_OPENMP) && defined(TRAME_USE_OPENMP)
    #undef TRAME_USE_OPENMP

    #undef TRAME_DONE_USE_OPENMP
    #define TRAME_DONE_USE_OPENMP
#endif

#ifdef TRAME_USE_OPENMP
    // #include "omp.h" //  OpenMP
    #ifndef ARMA_USE_OPENMP
        #define ARMA_USE_OPENMP
    #endif
#endif

#ifdef TRAME_DONT_USE_OPENMP
    #ifdef TRAME_USE_OPENMP
        #undef TRAME_USE_OPENMP
    #endif

    #ifndef ARMA_DONT_USE_OPENMP
        #define ARMA_DONT_USE_OPENMP
    #endif
#endif

#ifndef TRAME_DEFAULT_SIM_DRAWS
    #define TRAME_DEFAULT_SIM_DRAWS 1000;
#endif

namespace trame
{
    using uint_t = unsigned int;
}

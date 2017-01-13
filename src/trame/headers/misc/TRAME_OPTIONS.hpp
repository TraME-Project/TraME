/*################################################################################
  ##
  ##   Copyright (C) 2015 - 2017 the TraME Team:
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

#pragma once

#ifndef TRAME_USE_GUROBI_C
    #define TRAME_USE_GUROBI_C
#endif

#ifdef TRAME_RCPP_ARMADILLO
    #include <RcppArmadillo.h>
#else
    #include "armadillo"
#endif

#ifndef TRAME_DEFAULT_SIM_DRAWS
    #define TRAME_DEFAULT_SIM_DRAWS 1000;
#endif
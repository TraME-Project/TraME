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
 * Slight modification of the zeroin.c file from Netlib
 *
 * Keith O'Hara
 * 05/10/2016
 */

#ifndef _zeroin_HPP
#define _zeroin_HPP

#include <cmath> // for abs value
#include <limits>
#include "../misc/TRAME_OPTIONS.hpp"
#include "trame_structs.hpp" // for abs value

double zeroin(double ax, double bx, double (*f)(double x, const trame_zeroin_data& opt_data), const trame_zeroin_data& zeroin_data, double* tol_inp, int* max_iter_inp);

#endif

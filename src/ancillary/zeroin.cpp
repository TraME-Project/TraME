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

/*
 * Slight modification of the zeroin.c file from Netlib
 *
 * Keith O'Hara
 * 05/10/2016
 *
 * This version:
 * 07/17/2017
 */

#include "trame.hpp"

double 
trame::zeroin(const double ax, const double bx, std::function<double (const double val_inp, void* opt_data)> zero_objfn, void* opt_data, const double* tol_inp, const int* max_iter_inp)
{
    const double tol = (tol_inp) ? *tol_inp : 1E-12;
    const int max_iter = (max_iter_inp) ? *max_iter_inp : 10000;
    //
    double fa = zero_objfn(ax,opt_data);
    double fb = zero_objfn(bx,opt_data);

    // check endpoints
    if (fa == 0.0) {
        return ax;
    }
    if (fb == 0.0) {
        return bx;
    }

    double a = ax, b = bx, c = ax;
    double fc = fa;
    
    //
    int iter = 0;
    
    const double eps_temp = std::numeric_limits<double>::epsilon();
    const double tol_act = 2*eps_temp*std::abs(b) + tol/2;

    double p, q, prev_step, new_step;
    double t1, cb, t2;
    
    new_step = (c - b)/2;
    
    while (std::abs(new_step) > tol_act && iter < max_iter) {
        iter++;
        prev_step = b - a;
            
        if ( std::abs(fc) < std::abs(fb) ) {
            a = b;  b = c;  c = a;
            fa=fb;  fb=fc;  fc=fa;
        }
        
        new_step = (c - b)/2;

        if ( std::abs(prev_step) >= tol_act && std::abs(fa) > std::abs(fb) ) {
            
            cb = c - b;

            if ( a==c ) {
                t1 = fb/fa;
                p = cb*t1;
                q = 1.0 - t1;
            } else {
                q = fa/fc;  t1 = fb/fc;  t2 = fb/fa;
                p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
                q = (q-1.0) * (t1-1.0) * (t2-1.0);
            }

            if ( p > 0.0 ) {
                q = -q;
            } else {
                p = -p;
            }

            if ( p < (0.75*cb*q-std::abs(tol_act*q)/2) && p < std::abs(prev_step*q/2) ) {
                new_step = p/q;
            }
        }

        if ( std::abs(new_step) < tol_act ) {
            if ( new_step > 0.0 ) {
                new_step = tol_act;
            } else {
                new_step = -tol_act;
            }
        }

        a = b;
        fa = fb;

        b += new_step;  
        fb = zero_objfn(b,opt_data);

        if ( (fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0) ) {
            c = a;
            fc = fa;
        }
    }
    //
    return b;
}


#ifndef _inv_pwa_R_H
#define _inv_pwa_R_H
#include <R.h>
#include <math.h>
#include "trame.hpp"

RcppExport SEXP inv_pwa_R(SEXP a_R, SEXP B_R, SEXP C_R, SEXP k_R);

#endif
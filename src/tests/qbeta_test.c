/*
 * test for qbeta.c from the Rnmath library
 *
 * cd ~/Desktop/SCM/GitHub/TraME/src/tests
 * clang -I/Builds/Rnmath/include qbeta_test.c -o qbeta.test -L/Builds/Rnmath/lib -lRmath
 * then add 
 * export DYLD_FALLBACK_LIBRARY_PATH=/Builds/Rnmath/lib
 * before running ./qbeta.test
 */
 
#include <stdio.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>

int main()
{
    double p = 0.2;
    double shape1 = 2.0;
    double shape2 = 3.0;
    int lwrtail = 1;
    int log_p = 0;
    //
    double ret_val = qbeta(p,shape1,shape2,lwrtail,log_p);
    printf("%lf", ret_val);
    //
    return 0;
}
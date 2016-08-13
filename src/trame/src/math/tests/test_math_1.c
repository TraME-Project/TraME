/*
 *
 * cd ~/Desktop/SCM/GitHub/TraME/src/math/tests
 * gcc-mp-5 -O2 test_math_1.c -c -o test_math_1.o
 * gfortran -O2 ../fn.f90  -c -o ../fn.o
 * gcc-mp-5 -o test_math_1.test ../fn.o test_math_1.o -lgfortran
 *
 * gfortran-mp-5 -O2 ../fn.f90  -c -o ../fn.o
 */

#include <stdio.h>

extern double r8_gamma_(double* x);

extern double r8_lbeta_(double* a, double* b);

int main()
{
	double x[] = {3.0};
    
	double res_gamma = r8_gamma_(x);
	
	printf("gamma fn: %.4f\n", res_gamma);
    
	//
	
	double a[] = {2.0};
	double b[] = {3.0};
	
	double res_lbeta = r8_lbeta_(a,b);
	
	printf("gamma fn: %.4f\n", res_lbeta);
	
	return 0;
}
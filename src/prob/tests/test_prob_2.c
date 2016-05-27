/*
 *
 * Test integrated beta quantile function 
 *
 * cd ~/Desktop/SCM/GitHub/TraME/src/prob/tests
 * gcc-mp-5 -O2 test_prob_2.c -c -o test_prob_2.o
 * gfortran-mp-5 -O2 ../prob.f90  -c -o ../prob.o
 * gfortran-mp-5 -O2 ../../math/quadpack_double.f90  -c -o ../../math/quadpack_double.o
 * gfortran-mp-5 -O2 ../aux.f90  -c -o ../aux.o
 * gcc-mp-5 -o test_prob_2.test ../prob.o ../../math/quadpack_double.o ../aux.o test_prob_2.o -lgfortran
 *
 *
 */

#include <stdio.h>

extern void beta_cdf_inv_int_(double* x, double* p, double* q, double* result);

int main()
{
	double result[1];
	
	double x[] = {0.7};
	double p[] = {2.0};
	double q[] = {3.0};
    
	beta_cdf_inv_int_(x,p,q,result);
	
	printf("Inv. Beta: %.4f\n", result[0]);
	
	return 0;
}
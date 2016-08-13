/*
 *
 * cd ~/Desktop/SCM/GitHub/TraME/src/prob/tests
 * gcc-mp-5 -O2 test_prob_1.c -c -o test_prob_1.o
 * gfortran -O2 ../prob.f90  -c -o ../prob.o
 * gcc-mp-5 -o test_prob_1.test ../prob.o test_prob_1.o -lgfortran
 *
 *
 */

#include <stdio.h>

extern void beta_cdf_inv_(double* cdf, double* p, double* q, double* x);
extern void beta_cdf_(double* x, double* a, double* b, double* cdf);

int main()
{
	double x[1];
	
	double cdf[] = {0.1808};
	double p[]   = {2.0};
	double q[]   = {3.0};
    
	beta_cdf_inv_(cdf,p,q,x);
	
	printf("Inv. Beta: %.4f\n", x[0]);
    
    //
    
    double cdf_2[1];
	
	double x_2[] = {0.2};
	double a[]   = {2.0};
	double b[]   = {3.0};
    
	beta_cdf_(x_2,a,b,cdf_2);
	
	printf("Beta: %.4f\n", cdf_2[0]);
	
	return 0;
}
/*
 * Test integrated beta quantile function 
 *
 * cd ~/Desktop/SCM/GitHub/TraME/src/prob/tests
 * g++-mp-5 -O2 -fopenmp test_prob_2.cpp -c -o test_prob_cpp_2.o
 * gfortran-mp-5 -O2 ../prob.f90  -c -o ../prob.o
 * gfortran-mp-5 -O2 ../../math/quadpack_double.f90  -c -o ../../math/quadpack_double.o
 * gfortran-mp-5 -O2 ../aux.f90  -c -o ../aux.o
 * g++-mp-5 -o test_prob_cpp_2.test ../prob.o ../../math/quadpack_double.o ../aux.o test_prob_cpp_2.o -lgfortran -fopenmp
 */

#include <stdio.h>
#include "omp.h"

extern "C" {
	void beta_cdf_inv_int_(double* x, double* p, double* q, double* result);
}

int main()
{    
	double result[1];
	
	double x[] = {0.7};
	double p[] = {2.0};
	double q[] = {3.0};
    
	beta_cdf_inv_int_(x,p,q,result);
	
	printf("Inv. Beta: %.4f\n", result[0]);
	
	//
    
    int i;
    int N_test = 100;
    double a = 0.01, b = 0.99;
    double seqby = (b - a)/(N_test-1);
    
    double* x_seq = new double[N_test];
    
    for (i=0; i<N_test; i++) {
        x_seq[i] = a + i*seqby;
    }
	
	//#pragma omp parallel for num_threads(2) private(result)
		for (i=0; i<N_test; i++) {
			beta_cdf_inv_int_(&x_seq[i],p,q,result);
            
            printf("Inv. Beta: %.4f\n", result[0]);
		}
	
    //
	return 0;
}
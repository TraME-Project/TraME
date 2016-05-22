/*
 * generic_LP test: Linear programming
 *
 * Keith O'Hara
 * 05/17/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/tests
 * clang++ -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include ktest_1.cpp -c -o ktest_1.o
 * clang++ -O2 -Wall -o cpptest.test ../generic_lp.o ktest_1.o -L/Library/gurobi650/mac64/lib -lgurobi65 -framework Accelerate
 */

#include "armadillo"

extern "C" {
#include "../generic_lp.h"
}

int main()
{
    double  obj[]   = {2, 4, 3};
    double  A[3][3] = {{3, 4, 2}, {2, 1, 2}, {1, 3, 2}};
    char    sense[] = {'<', '<', '<'};
    double  rhs[]   = {60, 40, 80};
    double  lb[]    = {0, 0, 0};
    double  sol_mat_1[3];
    double  sol_mat_2[3];
    double  dual_mat_1[3];
    double  dual_mat_2[3];
    int     modelSense = 1; // maximize
    int     solved;
    double  objval;

    /* Solve the model */

    solved = generic_LP_c(3, 3, obj, &A[0][0], modelSense, rhs, sense, NULL, lb,
                          NULL, &objval, sol_mat_1, sol_mat_2, dual_mat_1, dual_mat_2);
    //
    std::cout << sol_mat_1[0] << sol_mat_1[1] << sol_mat_1[2] << std::endl;
    //  
    return 0;
}
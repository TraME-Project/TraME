/*
 * test for calling C-version of generic_lp
 *
 * Keith O'Hara
 * 05/17/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/lp/tests
 * clang++ -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include c_with_cpp_test_1.cpp -c -o c_with_cpp_test_1.o
 * clang++ -O2 -Wall -o cpptest_1.test ../generic_lp.o c_with_cpp_test_1.o -L/Library/gurobi650/mac64/lib -lgurobi65 -framework Accelerate
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

    solved = generic_LP_C(3, 3, obj, &A[0][0], modelSense, rhs, sense, NULL, lb,
                          NULL, &objval, sol_mat_1, sol_mat_2, dual_mat_1, dual_mat_2);
    //
    std::cout << "x_1: " << sol_mat_1[0] << ". x_2: " << sol_mat_1[1] << ". x_3: " << sol_mat_1[2] << ".\n" << std::endl;
    //  
    return 0;
}
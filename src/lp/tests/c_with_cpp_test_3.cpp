/*
 * test for calling C-version of generic_lp
 *
 * Keith O'Hara
 * 05/17/2016
 * 
 * cd ~/Desktop/SCM/GitHub/TraME/src/lp/tests
 * clang++ -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include c_with_cpp_test_3.cpp -c -o c_with_cpp_test_3.o
 * clang++ -O2 -Wall -o cpptest_3.test ../generic_lp.o c_with_cpp_test_3.o -L/Library/gurobi650/mac64/lib -lgurobi65 -framework Accelerate
 */

#include "armadillo"

extern "C" {
#include "../generic_lp.h"
}

int main()
{
    double  obj[]   = {1, 1, 0};
    double  Q[3][3] = {{1, 1, 0}, {0, 1, 1}, {0, 0, 1}};
    double  A[2][3] = {{1, 2, 3}, {1, 1, 0}};
    char    sense[] = {'>', '>'};
    double  rhs[]   = {4, 1};
    double  lb[]    = {0, 0, 0};
    double  sol_mat_1[3];
    double  sol_mat_2[3];
    double  dual_mat_1[2];
    double  dual_mat_2[2];
    int     modelSense = 0; // minimize
    int     solved;
    double  objval;

    /* Solve the model */

    solved = generic_LP_C(2, 3, obj, &A[0][0], modelSense, rhs, sense, &Q[0][0], lb,
                          NULL, &objval, sol_mat_1, sol_mat_2, dual_mat_1, dual_mat_2);
    //
    std::cout << "x_1: " << sol_mat_1[0] << ". x_2: " << sol_mat_1[1] << ". x_3: " << sol_mat_1[2] << ".\n" << std::endl;
    //  
    return 0;
}
/* Copyright 2016, Gurobi Optimization, Inc. */

/* This example formulates and solves the following simple QP model:

     minimize    x + y + x^2 + x*y + y^2 + y*z + z^2
     subject to  x + 2 y + 3 z >= 4
                 x +   y       >= 1

   The example illustrates the use of dense matrices to store A and Q
   (and dense vectors for the other relevant data).  We don't recommend
   that you use dense matrices, but this example may be helpful if you
   already have your data in this format.
*/

// cd ~/Desktop/SCM/GitHub/TraME/src/lp/tests
// clang -O2 -Wall -I/opt/local/include -I/Library/gurobi650/mac64/include generic_lp_c_1.c -o generic_lp_c_1.test -L/Library/gurobi650/mac64/lib -lgurobi65 -framework Accelerate

#include "../generic_lp_c.h"

int main()
{
  double  c[]     = {1, 1, 0};
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

  solved = generic_LP_C(2, 3, c, &A[0][0], modelSense, rhs, sense, &Q[0][0], lb,
                        NULL, &objval, sol_mat_1, sol_mat_2, dual_mat_1, dual_mat_2);

  if (solved) {
    printf("Solved X: x=%.4f, y=%.4f, z=%.4f\n", sol_mat_1[0], sol_mat_1[1], sol_mat_1[2]);
    printf("Solved RC: x=%.4f, y=%.4f, z=%.4f\n", sol_mat_2[0], sol_mat_2[1], sol_mat_2[2]);
    printf("Dual: a=%.4f, b=%.4f\n", dual_mat_1[0], dual_mat_1[1]);
    printf("Dual: a=%.4f, b=%.4f\n", dual_mat_2[0], dual_mat_2[1]);
  }

  return 0;
}
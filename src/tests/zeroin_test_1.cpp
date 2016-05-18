//
// cd ~/Desktop/SCM/GitHub/TraME/src/tests
// clang++ -O2 -Wall -I/opt/local/include zeroin_test_1.cpp -o zeroin_test_1 -framework Accelerate
//

#include <iostream>
#include <cmath> // for abs-value

#include "../headers/zeroin.hpp"

struct my_data {
    double scalepar;
};

// linear function
double test_fn (double x, void *opt_data)
{
    double ret = x;
    return ret;
}

// cosine btwn [1,2]
double test_fn_2 (double x, void *opt_data)
{
    double ret = std::cos(x);
    return ret;
}

int main()
{
    double zero_tol = 0.00001;
    double max_iter = 500;
    
    //
    double lb_1 = -3;
    double ub_1 = 3;
    
    my_data root_data;
    root_data.scalepar = 1;
    
    double linear_zero = zeroin(lb_1, ub_1, test_fn, &root_data, zero_tol, max_iter);
    
    std::cout << "linear zero: " << linear_zero << std::endl;
    //

    double lb_2 = 1;
    double ub_2 = 2;

    double cos_zero = zeroin(lb_2, ub_2, test_fn_2, &root_data, zero_tol, max_iter);
    
    std::cout << "cosine [1,2] zero: " << cos_zero << std::endl;
    //
    return 0;
}
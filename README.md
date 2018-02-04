# TraME &nbsp; [![Build Status](https://travis-ci.org/TraME-Project/TraME.svg)](https://travis-ci.org/TraME-Project/TraME) [![Build Status](https://codecov.io/github/TraME-Project/TraME/coverage.svg?branch=master)](https://codecov.io/github/TraME-Project/TraME?branch=master)

<!-- ## Overview -->

TraME (**Tra**nsportation **M**ethods for **E**conometrics) is a C++ library for 
solving problems of equilibrium computation and estimation in consumer 
demand and matching frameworks via the Mass Transportation Approach.

## Installation and Testing

TraME can be installed via
```
./configure
make
make install
```

The last line will install TraME into `/usr/local`

There are several configuration options available:
* `-c` a coverage build (used with Codecov)
* `-d` a 'development' build with install names set to the build directory (as opposed to an install path)
* `-g` a debugging build (optimization flags set to: `-O0 -g`)
* `-l` specify the linear programming library to link against, either GLPK (`-l glpk`) or Gurobi (`-l gurobi`)
* `-m` specify the BLAS and Lapack libraries to link against; for example, `-m "-lopenblas"` or `-m "-framework Accelerate"`
* `-o` compiler optimization options; defaults to `-O3 -march=native -flto -DARMA_NO_DEBUG`
* `-p` enable OpenMP parallelization features

## Example

The following example constructs several matching function-based markets under different matching technologies, then solves for equilibrium using the iterated proportional fitting procedure (IPFP) algorithm.

``` cpp
#include "trame.hpp"

int main()
{
    int nbX = 18;      // number of x types
    int nbY = 5;       // number of y types
    double sigma = 1;  // scaling value

    arma::vec n = arma::ones(nbX,1); // number of agents of each type
    arma::vec m = arma::ones(nbY,1);

    // systematic utilities

    arma::mat alpha  = arma::randu(nbX,nbY);
    arma::mat gamma  = arma::randu(nbX,nbY);
    arma::mat lambda = 1 + arma::randu(nbX,nbY);
    arma::mat zeta   = 1 + arma::randu(nbX,nbY);

    arma::mat phi = alpha + gamma;

    arma::mat lambda_LTU = lambda/(lambda+zeta);
    arma::mat phi_LTU = (lambda%alpha + zeta%gamma) / (lambda+zeta);

    // build markets

    trame::mfe<trame::mmfs::geo> mfe_obj_TU(sigma,false);
    mfe_obj_TU.build(n,m,phi); // geometric matching function <=> perfectly transferable utility

    trame::mfe<trame::mmfs::cd> mfe_obj_LTU(sigma,false);
    mfe_obj_LTU.build(n,m,lambda_LTU,phi_LTU); // CD matching function <=> linearly transferable utility

    trame::mfe<trame::mmfs::min> mfe_obj_NTU(sigma,false);
    mfe_obj_NTU.build(n,m,alpha,gamma); // min matching function <=> non-transferable utility

    //

    arma::mat mu_TU, mu_LTU, mu_NTU;

    mfe_obj_TU.solve(mu_TU);
    arma::cout << "Solution of TU-logit problem using ipfp:\n" << mu_TU << arma::endl;

    trame::ipfp(mfe_obj_LTU,mu_LTU);
    arma::cout << "Solution of LTU-logit problem using ipfp:\n" << mu_LTU << arma::endl;

    trame::ipfp(mfe_obj_NTU,mu_NTU);
    arma::cout << "Solution of NTU-logit problem using ipfp:\n" << mu_NTU << arma::endl;
}
```

## Authors

Alfred Galichon and Keith O'Hara

## License

GPL (>= 2) 

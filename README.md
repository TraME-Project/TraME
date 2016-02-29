# TraME

Transportation Methods for Econometrics

### Overview

TraME (Transportation Methods for Econometrics) is an R package for for 
solving problems of equilibrium computation and estimation in consumer 
demand and matching frameworks via the Mass Transportation Approach.

### Status [![Build Status](https://travis-ci.org/TraME-Project/TraME.svg)](https://travis-ci.org/TraME-Project/TraME) [![Build status](https://ci.appveyor.com/api/projects/status/github/TraME-Project/TraME?branch=master)](https://ci.appveyor.com/project/alfredgalichon/trame/branch/master)

The package is under active development and should be considered as
`alpha stage' software.

### Installation and Testing

The quickest way to install TraME is via the devtools package.
```
install.packages("devtools")
library(devtools)
install_github("TraME-Project/TraME")
```
The TraME test routines are invoked as follows:
```
library(TraME)
library(gurobi)
tests_TraME()
```

Note that TraME requires compilation, and so an appropriate development environment is necessary to install the package.
* For Windows users, get Rtools.
* For Mac users, get Xcode and install gfortran 4.8.2.

### Authors

Alfred Galichon and the TraME team.

### License

GPL (>= 2)

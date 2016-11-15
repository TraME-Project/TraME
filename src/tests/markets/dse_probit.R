library(TraME)
rm(list=ls())
#library(gurobi)
nbX=5
nbY=3

n=rep(1,nbX)
m=rep(1,nbY)

phi =  matrix(runif(nbX*nbY),nrow=nbX)

rho = 0.5

probit_G <- new(probit_R)
probit_H <- new(probit_R)

probit_G$build(nbX,nbY,rho,TRUE)
probit_G$unifCorrelCovMatrices()

probit_H$build(nbY,nbX,rho,TRUE)
probit_H$unifCorrelCovMatrices()

dse_emp_obj_TU <- new(dse_empirical_R)
dse_emp_obj_TU$build_TU(n,m,phi,probit_G,probit_H,FALSE)

dse_emp_obj_TU$solve("cupidsLP")

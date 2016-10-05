library(TraME)
rm(list=ls())
#library(gurobi)

logit_obj <- new(logit_R)

U = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=TRUE)
mu = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=TRUE)

nbX = dim(U)[1]
nbY = dim(U)[2]
n = c(apply(mu,1,sum))

logit_obj$build(nbX,nbY)
logit_obj$U = U

sim_obj = logit_obj$simul(10000)

logit_obj$G(n)

resG = logit_obj$G(n,U)
resG
resGSim = sim_obj$G(n,U)
resGSim

resGstar = logit_obj$Gstar(n,resG$mu)
resGstar
resGstarSim = sim_obj$Gstar(n,resGSim$mu)
resGstarSim

mubar = matrix(2,2,3)

resGbar = logit_obj$Gbar(U,mubar,n)
resGbar
resGbarSim = sim_obj$Gbar(U,mubar,n)
resGbarSim
#
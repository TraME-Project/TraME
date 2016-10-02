library(TraME)
#library(gurobi)

probit_obj <- new(probit_R)

U = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=TRUE)
mu = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=TRUE)

nbX = dim(U)[1]
nbY = dim(U)[2]
n = c(apply(mu,1,sum))

rho = 0.5

probit_obj$build(nbX,nbY,rho,TRUE)
probit_obj$unifCorrelCovMatrices()

sim_obj = probit_obj$simul(10000)

resGSim = sim_obj$G(n,U)
resGSim

resGstarSim = sim_obj$Gstar(n,resGSim$mu)
resGstarSim

mubar = matrix(2,2,3)

resGbarSim = sim_obj$Gbar(U,mubar,n)
resGbarSim
#

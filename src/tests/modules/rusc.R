library(TraME)
#library(gurobi)

rusc_obj <- new(rusc_R)

U = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=TRUE)
mu = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=TRUE)

nbX = dim(U)[1]
nbY = dim(U)[2]
n = c(apply(mu,1,sum))

zeta = matrix(1,nbX,1) %*% matrix(c(0.1, 0.2, 0.3, 0),1,nbY+1)

rusc_obj$build(zeta,TRUE)
rusc_obj$U = U

sim_obj = rusc_obj$simul(10000)

rusc_obj$G(n)

resG = rusc_obj$G(n,U)
resG
resGSim = sim_obj$G(n,U)
resGSim

resGstar = rusc_obj$Gstar(n,resG$mu)
resGstar
resGstarSim = sim_obj$Gstar(n,resGSim$mu)
resGstarSim

mubar = matrix(2,2,3)

resGbar = rusc_obj$Gbar(U,mubar,n)
resGbar
resGbarSim = sim_obj$Gbar(U,mubar,n)
resGbarSim
#


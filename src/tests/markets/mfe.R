library(TraME)
rm(list=ls())
#library(gurobi)
#
nbX = 18
nbY = 5
sigma = 1

n = rep(1,nbX)
m = rep(1,nbY)

alpha = matrix(runif(nbX*nbY),nbX,nbY)
gamma = matrix(runif(nbX*nbY),nbX,nbY)

lambda = 1 + matrix(runif(nbX*nbY),nbX,nbY)
zeta = 1 + matrix(runif(nbX*nbY),nbX,nbY)

phi = alpha + gamma

lambda_LTU = lambda/(lambda + zeta)
phi_LTU = (lambda*alpha + zeta*gamma) / (lambda + zeta)
#
mfe_mmf_obj_TU <- new(mfe_mmf_R)
mfe_mmf_obj_TU$build_TU(n,m,phi,sigma,FALSE)

mfe_mmf_obj_LTU <- new(mfe_mmf_R)
mfe_mmf_obj_LTU$build_LTU(n,m,lambda_LTU,phi_LTU,sigma,FALSE)

mfe_mmf_obj_NTU <- new(mfe_mmf_R)
mfe_mmf_obj_NTU$build_NTU(n,m,alpha,gamma,sigma,FALSE)
#
mfe_mmf_obj_TU$solve()
mfe_mmf_obj_LTU$solve()
mfe_mmf_obj_NTU$solve()

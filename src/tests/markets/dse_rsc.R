library(TraME)
rm(list=ls())
#library(gurobi)
#
alpha = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=T)
gamma = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=T)
muhat = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=T)

n = c(1.2*apply(muhat,1,sum))
m = c(1.3*apply(t(muhat),1,sum))

nbX = length(n)
nbY = length(m)

lambda = 1 + alpha
zeta = 1 + gamma

phi = alpha + gamma

lambda_LTU = lambda/(lambda + zeta)
phi_LTU = (lambda*alpha + zeta*gamma) / (lambda + zeta)
#
dse_logit_obj_TU <- new(dse_logit_R)
dse_logit_obj_TU$build_TU(n,m,phi,FALSE)

dse_logit_obj_LTU <- new(dse_logit_R)
dse_logit_obj_LTU$build_LTU(n,m,lambda_LTU,phi_LTU,FALSE)

dse_logit_obj_NTU <- new(dse_logit_R)
dse_logit_obj_NTU$build_NTU(n,m,alpha,gamma,FALSE)
#
dse_logit_obj_TU$solve("jacobi")
#dse_logit_obj_LTU$solve()
dse_logit_obj_NTU$solve("darum")
#
arums_G = dse_logit_obj_NTU$get_arums_G()
arums_H = dse_logit_obj_NTU$get_arums_H()

arums_G$U = matrix(0,nbX,nbY)
dse_logit_obj_NTU$set_arums_G(arums_G)

arums_G2 = dse_logit_obj_NTU$get_arums_G()
arums_G2$U



nbX=5
nbY=3

n=rep(1,nbX)
m=rep(1,nbY)  

phi =  matrix(runif(nbX*nbY),nrow=nbX)

zetaG = matrix(1,nbX,1) %*% matrix(runif(nbY+1),1,nbY+1)
zetaH = matrix(1,nbY,1) %*% matrix(runif(nbX+1),1,nbX+1)

rscG = new(rsc_R)
rscH = new(rsc_R)

rscG$build_beta(zetaG,2,2)
rscH$build_beta(zetaH,2,2)

dse_rsc_obj_TU <- new(dse_rsc_R)
dse_rsc_obj_TU$build_TU(n,m,phi,rscG,rscH,FALSE)

dse_rsc_obj_TU$solve("maxWelfare")

m2 = build_market_TU_general(n,m,phi,rscG,rscH,FALSE)
m2Sim = build_market_TU_empirical(n,m,phi,rscG,rscH,nbDraws,seed)  

r2 = maxWelfare(m2,xFirst=T,notifications=T)
r2Sim = CupidsLP(m2Sim,xFirst=T,notifications=T)
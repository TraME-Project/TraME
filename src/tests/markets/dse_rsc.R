library(TraME)
rm(list=ls())
#library(gurobi)
nbX=5
nbY=3

n=rep(1,nbX)
m=rep(1,nbY)  

phi =  matrix(runif(nbX*nbY),nrow=nbX)

zetaG = matrix(1,nbX,1) %*% matrix(runif(nbY+1),1,nbY+1)
zetaH = matrix(1,nbY,1) %*% matrix(runif(nbX+1),1,nbX+1)
#
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
library(TraME)
#library(gurobi)

none_obj <- new(none_R)

U = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=TRUE)
mu = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=TRUE)

nbX = dim(U)[1]
nbY = dim(U)[2]
n = c(apply(mu,1,sum))

none_obj$build(nbX,nbY)
none_obj$U = U

none_obj$G(n)

resG = none_obj$G(n,U)
resG

resGstar = none_obj$Gstar(n,resG$mu)
resGstar

mubar = matrix(2,2,3)

resGbar = none_obj$Gbar(U,mubar,n)
resGbar
#
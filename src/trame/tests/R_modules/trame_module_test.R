library(TraME)
#library(gurobi)

logit_obj <- new(logit_R)

U = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=TRUE)
mu = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=TRUE)

nbX = dim(U)[1]
nbY = dim(U)[2]
n = c(apply(mu,1,sum))

logit_obj$build(nbX,nbY)
logit_obj$U = U

logit_obj$test_1(3,4)
logit_obj$test_add(3,4)
logit_obj$test_create(3,4)
logit_obj$test_mat_add(3,4)
logit_obj$test_2(3,4)
logit_obj$test_3(3,4)

logit_obj$G(n)

resG = logit_obj$G(n,U)

ans_star = logit_obj$Gstar(n,resG$mu)

mubar = matrix(2,2,3)

logit_obj$Gbar(U,mubar,n)

ans = logit_obj$simul(10)

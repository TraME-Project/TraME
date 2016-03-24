################################################################################
##
##   Copyright (C) 2015 - 2016 Alfred Galichon
##
##   This file is part of the R package TraME.
##
##   The R package TraME is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 2 of the License, or
##   (at your option) any later version.
##
##   The R package TraME is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with TraME. If not, see <http://www.gnu.org/licenses/>.
##
################################################################################

test_ipfp <- function(seed=777, nbX=18, nbY=5)
{
  noSingles = (nbX == nbY)
  set.seed(seed)
  tm = proc.time()
  #
  message('*===================   Start of test_ipfp   ===================*\n')
  #
  n=rep(1,nbX)
  m=rep(1,nbY)
  
  alpha = matrix(runif(nbX*nbY),nrow=nbX)
  gamma = matrix(runif(nbX*nbY),nrow=nbX)
  lambda = matrix(1+runif(nbX*nbY),nrow=nbX)
  zeta = matrix(1+runif(nbX*nbY),nrow=nbX)
  
  m1 = build_market_TU_logit(n,m,alpha+gamma,neededNorm = defaultNorm(noSingles))
  m2 = build_market_NTU_logit(n,m,alpha,gamma,neededNorm = defaultNorm(noSingles))
  m3 = build_market_LTU_logit(n,m,lambda/(lambda+zeta),(lambda*alpha+zeta*gamma)/(lambda+zeta),neededNorm = defaultNorm(noSingles))
  #
  r1 = ipfp(m1,xFirst=TRUE,notifications=TRUE)
  message("Solution of TU-logit problem using ipfp:")
  print(c(r1$mu))
  message("")
  #
  if (!noSingles) {
    r2 = ipfp(m2,xFirst=TRUE,notifications=TRUE)
    message("Solution of NTU-logit problem using ipfp:")
    print(c(r2$mu))
    message("")
  }
  #
  r3 = ipfp(m3,xFirst=TRUE,notifications=TRUE)
  message("Solution of LTU-logit problem using parallel ipfp:")
  print(c(r3$mu))
  #
  time <- proc.time() - tm
  message(paste0('\nEnd of test_ipfp. Time elapsed = ', round(time["elapsed"],5), 's.\n'))
  #
  ret <- c(r1$mu,r3$mu)
  return(ret)
}

test_newton <- function(seed=777, nbX=5, nbY=3, nbDraws=1e3)
{
    set.seed(seed)
    tm = proc.time()
    #
    message('*===================   Start of test_newton   ===================*\n')
    #
    n = rep(1,nbX)
    m = rep(1,nbY)
    
    phi =  matrix(runif(nbX*nbY),nrow=nbX)
    #
    m1 = build_market_TU_logit(n,m,phi)
    r1 = ipfp(m1,xFirst=TRUE,notifications=TRUE)
    r1bis = newton(m1,xFirst=TRUE,notifications=TRUE)
    #
    message("Solution of TU-logit:")
    message("mu using (i) IPFP and (ii) newton:")
    print(r1$mu)
    print(r1bis$mu)
    #
    message("U using (i) IPFP and (ii) newton:")
    print(c(r1$U)[1:min(5,nbX*nbY)])
    print(c(r1bis$U)[1:min(5,nbX*nbY)])
    #
    message("V using (i) IPFP and (ii) newton:")
    print(c(r1$V)[1:min(5,nbX*nbY)])
    print(c(r1bis$V)[1:min(5,nbX*nbY)])  
    #
    time = proc.time() - tm
    message(paste0('\nEnd of test_newton. Time elapsed = ', round(time["elapsed"],5), 's.\n')) 
    #
    ret <- c(r1$mu,r1bis$mu,r1$U,r1bis$U,r1$V,r1bis$V)
    return(ret)
}

test_maxWelfare = function(seed=777, nbX=5, nbY=3, nbDraws=1e3)
{
    set.seed(seed)
    tm = proc.time()
    #
    message('*===================   Start of test_maxWelfare   ===================*\n')
    #
    n=rep(1,nbX)
    m=rep(1,nbY)  
    
    phi =  matrix(runif(nbX*nbY),nrow=nbX)
    #
    m1 = build_market_TU_logit(n,m,phi)
    r1 = ipfp(m1,xFirst=TRUE,notifications=TRUE)
    r1bis = maxWelfare(m1,xFirst=TRUE,notifications=TRUE)
    #
    message("Solution of TU-logit:")
    #
    message("mu using (i) IPFP and (ii) maxWelfare:")
    print(r1$mu)
    print(r1bis$mu)
    #
    message("U using (i) IPFP and (ii) maxWelfare:")
    print(c(r1$U)[1:min(5,nbX*nbY)])
    print(c(r1bis$U)[1:min(5,nbX*nbY)])
    #
    message("V using (i) IPFP and (ii) maxWelfare:")
    print(c(r1$V)[1:min(5,nbX*nbY)])
    print(c(r1bis$V)[1:min(5,nbX*nbY)])
    message("")
    #
    zetaG = matrix(1,nbX,1) %*% matrix(runif(nbY+1),1,nbY+1)
    zetaH = matrix(1,nbY,1) %*% matrix(runif(nbX+1),1,nbX+1)
    
    rscG = build_RSCbeta(zetaG,2,2)
    rscH = build_RSCbeta(zetaH,2,2)
    
    m2 = build_market_TU_general(n,m,phi,rscG,rscH)
    m2Sim = build_market_TU_empirical(n,m,phi,rscG,rscH,nbDraws,seed)  
    
    r2 = maxWelfare(m2,xFirst=T,notifications=T)
    r2Sim = CupidsLP(m2Sim,xFirst=T,notifications=T)
    #
    message("Solution of TU-RSCbeta:")
    #
    message("val using (i) LP and (ii) maxWelfare:")
    print(r2Sim$val)
    print(r2$val)
    #
    message("mu using (i) LP and (ii) maxWelfare:")
    print(r2Sim$mu)
    print(r2$mu)
    #
    message("U using (i) LP and (ii) maxWelfare:")
    print(c(r2Sim$U)[1:min(5,nbX*nbY)])
    print(c(r2$U)[1:min(5,nbX*nbY)])
    #
    message("V using (i) LP and (ii) maxWelfare:")
    print(c(r2Sim$V)[1:min(5,nbX*nbY)])
    print(c(r2$V)[1:min(5,nbX*nbY)])  
    #
    time = proc.time() - tm
    message(paste0('\nEnd of test_maxWelfare. Time elapsed = ', round(time["elapsed"],5), 's.\n'))
    #
    ret <- c(r1$mu,r1bis$mu,r1$U,r1bis$U,r1$V,r1bis$V,r2Sim$val,r2$val,r2Sim$mu,r2$mu,r2Sim$U,r2$U,r2Sim$V,r2$V)
    return(ret)
}

test_jacobi <- function(nbDraws=1E3, seed=777, extensiveTesting = FALSE)
{
    set.seed(seed)
    ptm = proc.time()
    #
    message('*===================   Start of test_jacobi   ===================*\n')
    #
    alpha = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=T)
    gamma = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=T)
    muhat = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=T)
    
    n = c(1.2*apply(muhat,1,sum))
    m = c(1.3*apply(t(muhat),1,sum))
    #
    m1 = build_market_TU_logit(n,m,alpha+gamma)
    r1_jacobi = jacobi(m1)
    #
    message("Solution of TU-logit problem using Jacobi:")
    print(c(r1_jacobi$mu))
    message("")
    #
    m2 = build_market_NTU_logit(n,m,alpha,gamma)
    r2_jacobi = jacobi(m2)
    #
    message("Solution of NTU-logit problem using Jacobi:")
    print(c(r2_jacobi$mu))
    #
    if(extensiveTesting==TRUE){
        nbX = length(n)
        nbY = length(m)
        
        logitM = build_logits(nbX,nbY)
        logitW = build_logits(nbY,nbX)
        
        logitSimM = simul(logitM,nbDraws,seed)
        logitSimW = simul(logitW,nbDraws,seed)
        #
        m2Sim = build_market_NTU_general(n,m,alpha,gamma,logitSimM,logitSimW)
        r2Sim_jacobi = jacobi(m2Sim)
        #
        message("\nSolution of NTU-logitSim problem using Jacobi:")
        print(c(r2Sim_jacobi$mu))
    }
    #
    time = proc.time()-ptm
    message(paste0('\nEnd of test_jacobi. Time elapsed = ', round(time["elapsed"],5), 's.\n'))
    #
    ret <- c(r1_jacobi$mu,r2_jacobi$mu)
    return(ret) 
}

test_darum <- function(nbDraws=1E3,seed=777)
{
    set.seed(seed)
    ptm = proc.time()
    #
    message('*===================   Start of test_darum   ===================*\n')
    #
    alpha = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=T)
    gamma = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=T)
    muhat = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=T)
    
    n = c(1.2*apply(muhat,1,sum))
    m = c(1.3*apply(t(muhat),1,sum))
    #
    nbX = length(n)
    nbY = length(m)
    
    logitM = build_logits(nbX,nbY)
    logitW = build_logits(nbY,nbX)  
    
    logitSimM = simul(logitM,nbDraws,seed)
    logitSimW = simul(logitW,nbDraws,seed)  
    #
    m1 = build_market_NTU_logit(n,m,alpha,gamma)
    r1 = darum(m1,TRUE,TRUE)
    #
    message("Solution of NTU-logit problem using Jacobi:")  
    print(c(r1$mu))
    message("")
    #
    m1Sim = build_market_NTU_general(n,m,alpha,gamma,logitSimM,logitSimW)
    r1Sim = darum(m1Sim,TRUE,TRUE)
    #
    message("Solution of NTU-logitSim problem using Jacobi:")  
    print(c(r1Sim$mu))
    #
    time = proc.time() - ptm
    message(paste0('\nEnd of test_darum. Time elapsed = ', round(time["elapsed"],5), 's.\n'))
    #
    ret <- c(r1$mu,r1Sim$mu)
    return(ret)
}

test_cupidsLP <- function(nbX=5, nbY=3, nbDraws=1E3, seed=777)
{
    set.seed(seed)
    ptm = proc.time()
    #
    message('*===================   Start of test_cupidsLP   ===================*\n')
    #
    alpha = matrix(runif(nbX*nbY),nrow=nbX)
    gamma = matrix(runif(nbX*nbY),nrow=nbX)
    
    n = rep(1,nbX)
    m = rep(1, nbY)
    
    logitM = build_logits(nbX,nbY)
    logitW = build_logits(nbY,nbX)
    
    logitSimM = simul(logitM,nbDraws,seed)
    logitSimW = simul(logitW,nbDraws,seed)
    #
    m1Sim = build_market_TU_general(n,m,alpha+gamma,logitSimM,logitSimW)
    r1SimSmart = CupidsLP(m1Sim, TRUE, TRUE)
    #
    message("Solution of TU-logitSim problem using LP:")
    print(c(r1SimSmart$mu))
    message("")
    #
    probitMth = build_probit(unifCorrelCovMatrices(nbX,nbY,0.3))
    probitWth =  build_probit(unifCorrelCovMatrices(nbY,nbX,0.3))
    
    probitM = simul(probitMth,nbDraws,seed)
    probitW = simul(probitWth,nbDraws,seed)
    #
    mktTUProbit = build_market_TU_general(n,m,alpha+gamma,probitM,probitW)
    rTUProbit = CupidsLP(mktTUProbit, TRUE, TRUE)
    #
    message("Solution of TU-Probit problem using LP:")
    print(c(rTUProbit$mu))
    #
    time = proc.time() - ptm
    message(paste0('\nEnd of test_cupidsLP. Time elapsed = ', round(time["elapsed"],5), 's.\n')) 
    #
    ret <- c(r1SimSmart$mu,rTUProbit$mu)
    return(ret)
}

test_oapLP <- function(nbX=8,nbY=5,seed=777)
{
    set.seed(seed)
    ptm <- proc.time()
    #
    message('*===================   Start of test_plainOAP   ===================*\n')
    #
    n=rep(1,nbX)
    m=rep(1,nbY)
    
    alpha = matrix(runif(nbX*nbY),nrow=nbX)
    gamma = matrix(runif(nbX*nbY),nrow=nbX)
    #
    mkt = build_market_TU_none(n,m,alpha+gamma)
    res = oapLP(mkt,TRUE, TRUE)
    #
    #   print('mu:')
    #   print(res$mu)
    message('u:')
    print(res$u)
    message('v:')
    print(res$v)
    #
    time <- proc.time()-ptm
    message(paste0('\nEnd of test_oapLP. Time elapsed = ', round(time["elapsed"],5), 's.\n'))
    #
    ret <- c(res$u,res$v)
}  

test_eapNash <- function(nbX=8,nbY=5,seed=777,debugmode = FALSE)
{
    set.seed(seed)
    ptm <- proc.time()
    #
    message('*===================   Start of test_nashITU   ===================*\n')
    #
    n=rep(1,nbX)
    m=rep(1,nbY)
    
    alpha = matrix(runif(nbX*nbY),nrow=nbX)
    gamma = matrix(runif(nbX*nbY),nrow=nbX)
    lambda = matrix(1+runif(nbX*nbY),nrow=nbX)
    zeta = matrix(1+runif(nbX*nbY),nrow=nbX)
    #
    mkt = build_market_LTU_none(n,m,lambda/(lambda+zeta),(lambda*alpha+zeta*gamma)/(lambda+zeta) )
    #mkt = build_market_TU_none(n,m,alpha,gamma)
    #
    nash1 = eapNash(mkt,TRUE,TRUE)
    unash1 = nash1$u
    vnash1 = nash1$v
    #
    nash2 = eapNash(mkt,FALSE,TRUE)
    unash2 = nash2$u
    vnash2 = nash2$v
    #
    message("u[x] (upper and lower):")
    print(matrix(c(unash1,unash2),nrow=2,byrow=TRUE))
    message("v[y] (lower and upper):")
    print(matrix(c(vnash1,vnash2),nrow=2,byrow=TRUE)) 
    # print(c("Method","u1","u2","u3","v1","v2","v3","val"))
    # print(c("na1",round(c(unash1[1:3],vnash1[1:3]),2),round(sum(n*unash1)+sum(m*vnash1),4)))
    # print(c("na2",round(c(unash2[1:3],vnash2[1:3]),2),round(sum(n*unash2)+sum(m*vnash2),4)))
    # print("---------")
    #
    time <- proc.time() - ptm
    message(paste0('\nEnd of test_eapNash. Time elapsed = ', round(time["elapsed"],5), 's.\n'))
    #
    ret <- c(c(unash1,unash2),c(vnash1,vnash2))
    return(ret)
}  

tests_equilibrium = function(notifications=TRUE,nbDraws=1e3){
    ptm = proc.time()
    #
    res_darum  <- test_darum(nbDraws=nbDraws)
    res_ipfp   <- test_ipfp()
    res_newton <- test_newton(nbDraws=nbDraws)
    res_maxW   <- test_maxWelfare(nbDraws=nbDraws)
    res_jacobi <- test_jacobi(nbDraws=nbDraws)
    res_CLP    <- test_cupidsLP(nbDraws=nbDraws)
    res_oapLP  <- test_oapLP()
    res_nash   <- test_eapNash()
    # MD5 checksum
    res_all <- round(c(res_darum,res_ipfp,res_newton,res_maxW,res_jacobi,res_CLP,res_oapLP,res_nash),5)
    res_md5 <- digest(res_all,algo="md5")
    #
    time = proc.time() - ptm
    if(notifications){
        message(paste0('All tests of Equilibrium completed. Overall time elapsed = ', round(time["elapsed"],5), 's.'))
    }
    #
    return(res_md5)
}


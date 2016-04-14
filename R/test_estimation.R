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

test_loglikelihood <- function(seed=777, nbX=5, nbY=4, dX=3, dY=2)
{  
    set.seed(seed)
    tm = proc.time()
    #
    message('*===================   Start of test_loglikelihood   ===================*\n')
    #
    n = rep(1,nbX)
    m = rep(1,nbY)  
    
    xs = matrix(runif(nbX*dX),nrow=nbX)
    ys = matrix(runif(nbY*dY),nrow=nbY)
    #
    logitM = build_logits(nbX,nbY)
    logitW = build_logits(nbY,nbX)
    # logitSimM = simul(logitM,1,seed)
    # logitSimW = simul(logitW,1,seed)
    #
    A = matrix(1,nrow=dX,ncol=dY)
    phi = xs %*% A %*% t(ys) 
    #
    mktSim = build_market_TU_empirical(n,m,phi,logitM,logitW,1,seed)  
    muhat = CupidsLP(mktSim, T, F)$mu
    muhatx0  = n-apply(muhat,1,sum)
    muhat0y  = m-apply(muhat,2,sum)
    #
    affinitymodel = buildModel_affinity(array( kronecker(xs,ys) , dim = c(nbX,nbY,dX *dY)),n,m)
    theta0=initparam(affinitymodel)$param
    market = parametricMarket(affinitymodel,theta0)
    dtheta = diag(affinitymodel$nbParams)
    #dtheta = matrix(0.1,nrow=6,ncol=1)
    #
    tsf = proc.time()  
    #mudmu = dtheta_mu_default(affinitymodel,market,theta0,dtheta)
    mudmu = dtheta_mu_default(affinitymodel,market,theta0)
    message(paste0('Time elapsed - general: ', round((proc.time()-tsf)["elapsed"],5), 's.')) 
    #
    tsf = proc.time()  
    dmunum = dtheta_mu_numeric(affinitymodel,market,theta0,dtheta)
    message(paste0('Time elapsed - numerical: ', round((proc.time()-tsf)["elapsed"],5), 's.')) 
    #
    tsf = proc.time()    
    dmulogit=dtheta_mu_logit(affinitymodel,market,theta0,dtheta)
    message(paste0('Time elapsed - logit: ', round((proc.time()-tsf)["elapsed"],5), 's.')) 
    #
    mu = mudmu$mu 
    #
    LL = sum(2* muhat * log(mudmu$mu)) + sum(muhatx0 * log(mudmu$mux0s)) + sum(muhat0y * log(mudmu$mu0ys)) 
    gradLL = apply( c( t( t(2* muhat/matrix(mudmu$mu,nrow=nbX) - muhatx0 / mudmu$mux0s)  - muhat0y / mudmu$mu0ys ) )* mudmu$dmu   ,2,sum) 
    gradLLlogit = apply( c( t( t(2* muhat/matrix(dmulogit$mu,nrow=nbX) - muhatx0 / dmulogit$mux0s)  - muhat0y / dmulogit$mu0ys ) )* dmulogit$dmu   ,2,sum) 
    gradLLnum = apply( c( t( t(2* muhat/matrix(dmunum$mu,nrow=nbX) - muhatx0 / dmunum$mux0s)  - muhat0y / dmunum$mu0ys ) )* dmunum$dmu   ,2,sum)
    #
    print(gradLL)
    print(gradLLlogit)
    print(gradLLnum)
    #
    time = proc.time() - tm  
    message(paste0('\nEnd of test_loglikelihood. Time elapsed = ', round(time["elapsed"],5), 's.\n'))
    #
    ret <- c(LL,gradLL,gradLLlogit,gradLLnum)
    return(ret)
}

test_mle <- function(seed=777, nbX=80, nbY=72, noiseScale=0.1, dX=3, dY=3)
{  
    set.seed(seed)
    tm = proc.time()
    #
    message('*===================   Start of test_mle   ===================*\n')
    #
    n = rep(1,nbX)
    m = rep(1,nbY)  
    
    xs = matrix(runif(nbX*dX),nrow=nbX)
    ys = matrix(runif(nbY*dY),nrow=nbY)
    #
    logitM = build_logits(nbX,nbY,1)
    logitW = build_logits(nbY,nbX,1)
    
    logitSimM = simul(logitM,50,seed)
    logitSimW = simul(logitW,50,seed)
    #
    A = matrix(1,nrow=dX,ncol=dY)
    phi = xs %*% A %*% t(ys) 
    #
    mktLogit = build_market_TU_logit(n,m,phi)
    noise = matrix(1+ noiseScale*rnorm(nbX*nbY),nrow=nbX)
    muhat = ipfp(mktLogit, T, F)$mu * noise
    #
    affinitymodel = buildModel_affinity(array( kronecker(xs,ys) , dim = c(nbX,nbY,dX *dY)),n=n,m=m)
    thetahat = c(matrix(mle(affinitymodel,muhat, print_level=0)$thetahat,nrow = dX, byrow=T))
    #
    message("Estimator:")
    print(thetahat)
    #
    time = proc.time() - tm
    message(paste0('\nEnd of test_mle. Time elapsed = ', round(time["elapsed"],5), 's.\n')) 
    #
    ret <- c(thetahat)
    return(ret)
}

test_mme <- function(seed=777, nbX=80, nbY=72, noiseScale=0.1, dX=3, dY=3)
{  
    set.seed(seed)
    tm = proc.time()
    #
    message('*===================   Start of test_mme   ===================*\n')
    #
    n = rep(1,nbX)
    m = rep(1,nbY)  
    
    xs = matrix(runif(nbX*dX),nrow=nbX)
    ys = matrix(runif(nbY*dY),nrow=nbY)
    #
    logitM = build_logits(nbX,nbY,1)
    logitW = build_logits(nbY,nbX,1)
    
    logitSimM = simul(logitM,50,seed)
    logitSimW = simul(logitW,50,seed)
    #
    A = matrix(1,nrow=dX,ncol=dY)
    phi = xs %*% A %*% t(ys)
    #
    mktLogit = build_market_TU_logit(n,m,phi)
    noise = matrix(1 + noiseScale*rnorm(nbX*nbY),nrow=nbX)
    muhat = ipfp(mktLogit, T, F)$mu * noise
    #
    affinitymodel = buildModel_affinity(array( kronecker(xs,ys) , dim =c( nbX,nbY, dX *dY )) ,n,m)
    thetahat = c(matrix(mle(affinitymodel,muhat, print_level=0)$thetahat,nrow = dX, byrow=T))
    #
    message("Estimator:")
    print(thetahat)  
    #
    time = proc.time() - tm  
    message(paste0('\nEnd of test_mme. Time elapsed = ', round(time["elapsed"],5), 's.\n')) 
    #
    ret <- c(thetahat)
    return(ret)
}

###################################################################################
###################################################################################
###################################################################################
# mLogLikelihoodBenchmark = function(theta)
# {
#   theA = matrix(theta,dX,dY)
#   thephi = xs %*% theA %*% t(ys)
#   thephiBis = matrix(affinitymodel$kron %*% c(theta),nrow=affinitymodel$nbX)
#   themkt = build_market_TU_logit(n,m,thephi)
#   res = ipfp(themkt,notifications=F)
#   themu = res$mu
#   thegrad = matrix(0,dX,dY)
#   for (i in 1:dX)
#     for (j in 1:dY)
#     {thegrad[i,j]= sum( t( (themu - muhat) * xs[,i]) * ys[,j]) 
#     }
#   mux0s  = n-apply(themu,1,sum)
#   mu0ys  = m-apply(themu,2,sum) 
#   entrop = 2 * sum(themu *log(themu) ) + sum(mux0s *log(mux0s)  - n*log(n)) +  sum(mu0ys*log(mu0ys)- m*log(m))
#   theval = sum( (themu - muhat) *thephi  ) - entrop  
#   return(list(objective=theval,gradient=c(thegrad) ))
# }

tests_estimation = function(notifications=TRUE,nbDraws=1e3)
{
    ptm = proc.time()
    #
    res_LL  <- round(test_loglikelihood(),5)
    res_mle <- round(test_mle(),5)
    res_mme <- round(test_mme(),5)
    
    res_all <- c(res_LL,res_mle,res_mme)
    # MD5 checksum
    res_LL_md5 <- digest(res_LL,algo="md5")
    res_mle_md5 <- digest(res_mle,algo="md5")
    res_mme_md5 <- digest(res_mme,algo="md5")
    
    res_all_md5 <- digest(res_all,algo="md5")
    #
    time = proc.time() - ptm
    if(notifications){
        message(paste0('All tests of Estimation completed. Overall time elapsed = ', round(time["elapsed"],5), 's.'))
    }
    #
    ret <- list(res_all_md5=res_all_md5,res_LL_md5=res_LL_md5,res_mle_md5=res_mle_md5,res_mme_md5=res_mme_md5)
    #
    return(ret)
}

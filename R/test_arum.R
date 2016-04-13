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

test_Logit <- function(nbDraws=1E4, seed=777, outsideOption=TRUE)
{
    set.seed(seed)
    ptm = proc.time()
    #
    message('*===================   Start of testLogit   ===================*\n')
    #
    U = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=TRUE)
    mu = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=TRUE)
    
    nbX = dim(U)[1]
    nbY = dim(U)[2]
    n = c(apply(mu,1,sum))
    #
    logits = build_logits(nbX,nbY,outsideOption=outsideOption)
    logitsSim = simul(logits,nbDraws,seed)
    #
    resG = G(logits,U,n)
    resGstar = Gstar(logits,resG$mu,n)
    resGSim = G(logitsSim,U,n)
    resGstarSim = Gstar(logitsSim,resGSim$mu,n)
    #
    message("(i) U and \\nabla G*(\\nabla G(U)) in (ii) cf and (iii) simulated logit:")
    print(c(U))
    print(c(resGstar$U))
    print(c(resGstarSim$U))
    
    message("\nG(U) in (i) cf and (ii) simulated logit:")
    print(resG$val)
    print(resGSim$val)
    
    message("\nG*(mu) in (i) cf and (ii) simulated logit:")
    print(resGstar$val)
    print(resGstarSim$val)
    #
    if(outsideOption){
        mubar = matrix(2,2,3)
        resGbar = Gbar(logits,U,n,mubar)
        resGbarS = Gbar(logitsSim,U,n,mubar)
        
        message("\nGbar(mu) in (i) cf and (ii) simulated logit:")
        print(resGbar$val)
        print(resGbarS$val)
    }
    #
    time = proc.time() - ptm
    message(paste0('\nEnd of testLogit. Time elapsed = ', time["elapsed"], 's.\n'))
    #
    ret <- c(resGstar$U,resGstarSim$U,resG$val,resGSim$val,resGstar$val,resGstarSim$val,resGbar$val,resGbarS$val)
    return(ret)
}

test_Probit <- function(nbDraws=1E4, seed=777, outsideOption=TRUE)
{
    set.seed(seed)
    ptm = proc.time()
    #
    message('*===================   Start of testProbit   ===================*\n')
    #
    U = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=T)
    mu = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=T)
    
    nbX = dim(U)[1]
    nbY = dim(U)[2]
    n = c(apply(mu,1,sum))
    
    rho = 0.5
    #
    Covar = unifCorrelCovMatrices(nbX,nbY,rho,outsideOption=outsideOption)
    probits = build_probit(Covar,outsideOption=outsideOption)
    probitsSim = simul(probits,nbDraws,seed)
    #
    resGSim = G(probitsSim,U,n)
    resGstarSim = Gstar(probitsSim,resGSim$mu,n)
    #
    message("(i) U and \\nabla G*(\\nabla G(U)) in simulated probit:")
    print(c(U))
    print(c(resGstarSim$U))
    #
    time = proc.time() - ptm
    message(paste0('\nEnd of testProbit. Time elapsed = ', time["elapsed"], 's.\n')) 
    #
    ret <- c(resGstarSim$U)
    return(ret)
}

test_RUSC <- function(nbDraws=1E4,seed=777)
{
    set.seed(seed)
    ptm = proc.time()
    #
    message('*===================   Start of test_RUSC   ===================*\n')
    #
    U = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=T)
    mu = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=T)
    
    nbX = dim(U)[1]
    nbY = dim(U)[2]
    n = c(apply(mu,1,sum)) + c(1,1)
    
    zeta = matrix(1,nbX,1) %*% matrix(c(0.1, 0.2, 0.3, 0),1,nbY+1)
    #
    RUSCs = build_RUSC(zeta)
    RUSCsSim = simul(RUSCs,nbDraws,seed)  
    #
    r1 = G(RUSCs,U,n)
    r1Sim = G(RUSCsSim,U,n)
    r2 = Gstar(RUSCs,r1$mu,n)
    r2Sim = Gstar(RUSCsSim,r1Sim$mu,n)
    #
    message("G(U) in (i) cf and (ii) simulated RUSC:")
    print(c(r1$val))
    print(c(r1Sim$val))  
    #
    message("\n\\nabla G(U) in (i) cf and (ii) simulated RUSC:")
    print(c(r1$mu))
    print(c(r1Sim$mu))  
    #
    message("\n(i) U and \\nabla G*(\\nabla G(U)) in (ii) cf and (iii) simulated RUSC:")
    message("(Note: in RUSC, (ii) should be approx equal to (iii) but not to (i).)")
    print(c(U))
    print(c(r2$U))
    print(c(r2Sim$U))
    #
    r3 = Gstar(RUSCs,mu,n)
    r3Sim = Gstar(RUSCsSim,mu,n)
    #
    message("\n\\nabla G*(mu) in (i) closed form and (ii) simulated RUSC:")
    print(c(r3$U))
    print(c(r3Sim$U))
    message("\nG*(mu) in (i) closed form and (ii) simulated RUSC:")
    print(c(r3$val))
    print(c(r3Sim$val))
    #
    r4 = G(RUSCs,r3$U, n)
    r4Sim = G(RUSCsSim,r3Sim$U,n)
    message("\n\\nabla G \\nabla G*(mu) in (i) closed form and (ii) simulated RUSC:")
    print(c(r4$mu))
    print(c(r4Sim$mu))
    #
    mubar = matrix(2,2,3)
    r5 = Gbar(RUSCs,U,n,mubar)
    r5Sim = Gbar(RUSCsSim,U,n,mubar)
    #
    message("\nGbar(U,mubar) in (i) cf and (ii) simulated RUSC:")
    print(r5$val)
    print(r5Sim$val)
    message("\n\\nabla Gbar(U,mubar) in (i) cf and (ii) simulated RUSC:")
    print(c(r5$mu))
    print(c(r5Sim$mu))
    #
    time = proc.time()-ptm
    message(paste0('\nEnd of test_RUSC. Time elapsed = ', time["elapsed"], 's.\n')) 
    #
    ret <- c(r1$val,r1Sim$val,r1$mu,r1Sim$mu,r2$U,r2Sim$U,r3$U,r3Sim$U,r3$val,r3Sim$val,r4$mu,r4Sim$mu,r5$val,r5Sim$val,r5$mu,r5Sim$mu)
    return(ret)
}

test_RSC <- function(nbDraws=1E4,seed=777)
{
    set.seed(seed)
    ptm = proc.time()
    #
    message('*===================   Start of test_RSC   ===================*\n')
    #
    U = matrix(c(1.6, 3.2, 1.1, 2.9, 1.0, 3.1),nrow=2,byrow=TRUE)
    mu = matrix(c(1, 3, 1, 2, 1, 3), nrow=2, byrow=TRUE)
    
    nbX = dim(U)[1]
    nbY = dim(U)[2]
    n = c(apply(mu,1,sum)) + c(1,1)
    
    zeta = matrix(1,nbX,1) %*% matrix(c(0.1, 0.2, 0.3, 0),1,nbY+1)
    #
    RSCs = build_RSCbeta(zeta,2,2)
    #RSCs = build_RSCnorm(zeta)
    RSCsSim = simul(RSCs,nbDraws,seed)  
    #
    r1 = G(RSCs,U,n)
    r1Sim = G(RSCsSim,U,n)
    r2 = Gstar(RSCs,r1$mu,n)
    r2Sim = Gstar(RSCsSim,r1Sim$mu,n)
    #
    message("G(U) in (i) cf and (ii) simulated RSC:")
    print(c(r1$val))
    print(c(r1Sim$val))  
    #
    message("\n\\nabla G(U) in (i) cf and (ii) simulated RSC:")
    print(c(r1$mu))
    print(c(r1Sim$mu))  
    #
    message("\n(i) U and \\nabla G*(\\nabla G(U)) in (ii) cf and (iii) simulated RSC:")
    message("(Note: in RSC, (ii) should be approx equal to (iii) but not to (i).)")
    print(c(U))
    print(c(r2$U))
    print(c(r2Sim$U))
    #
    r3 = Gstar (RSCs,mu, n)
    r3Sim = Gstar (RSCsSim,mu,n)
    #
    message("\n\\nabla G*(mu) in (i) closed form and (ii) simulated RSC:")
    print(c(r3$U))
    print(c(r3Sim$U))
    message("\nG*(mu) in (i) closed form and (ii) simulated RSC:")
    print(c(r3$val))
    print(c(r3Sim$val))
    #
    r4 = G (RSCs,r3$U, n)
    r4Sim = G (RSCsSim,r3Sim$U,n)
    message("\n\\nabla G \\nabla G*(mu) in (i) closed form and (ii) simulated RSC:")
    print(c(r4$mu))
    print(c(r4Sim$mu))
    #
    mubar = matrix(2,2,3)
    r5 = Gbar(RSCs,U,n,mubar)
    r5Sim = Gbar(RSCsSim,U,n,mubar)
    #
    message("\nGbar(U,mubar) in (i) cf and (ii) simulated RSC:")
    print(r5$val)
    print(r5Sim$val)
    message("\n\\nabla Gbar(U,mubar) in (i) cf and (ii) simulated RSC:")
    print(c(r5$mu))
    print(c(r5Sim$mu))
    #
    hess = D2Gstar.RSC(RSCs,mu,n)
    thef = function(themu) (Gstar(RSCs,themu,n)$val)
    hessNum = hessian(thef,mu)
    message("\nD^2G^* (i) in cf and (ii) using numerical hessian:")
    print(hess)
    print(round(hessNum,6))
    #
    time = proc.time()-ptm
    message(paste0('\nEnd of test_RSC. Time elapsed = ', time["elapsed"], 's.\n')) 
    #
    ret <- c(r1$val,r1Sim$val,r1$mu,r1Sim$mu,r2$U,r2Sim$U,r3$U,r3Sim$U,r3$val,r3Sim$val,r4$mu,r4Sim$mu,r5$val,r5Sim$val,r5$mu,r5Sim$mu)
    return(ret)
}

tests_arum <- function(notifications=TRUE,nbDraws=1e4)
{
    ptm = proc.time()
    #
    res_logit  <- round(test_Logit(nbDraws=nbDraws),5)
    res_probit <- round(test_Probit(nbDraws=nbDraws),5)
    res_RUSC   <- round(test_RUSC(nbDraws=nbDraws),5)
    res_RSC    <- round(test_RSC(nbDraws=nbDraws),5)
    
    res_all <- c(res_logit,res_probit,res_RUSC,res_RSC)
    # MD5 checksum
    res_logit_md5  <- digest(res_logit,algo="md5")
    res_probit_md5 <- digest(res_probit,algo="md5")
    res_RUSC_md5   <- digest(res_RUSC,algo="md5")
    res_RSC_md5    <- digest(res_RSC,algo="md5")
    
    res_all_md5 <- digest(res_all,algo="md5")
    #
    time = proc.time() - ptm
    #
    if(notifications){
        message(paste0('All tests of arums completed. Overall time elapsed = ', time["elapsed"], 's.'))
    }
    #
    ret <- list(res_all_md5=res_all_md5,res_logit_md5=res_logit_md5,res_probit_md5=res_probit_md5,
                res_RUSC_md5=res_RUSC_md5,res_RSC_md5=res_RSC_md5)
    #
    return(ret)
}

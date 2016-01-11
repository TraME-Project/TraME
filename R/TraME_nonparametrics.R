#################################################
###                  TraME:                   ###
###   Transportation Methods in Econometrics  ###
#################################################
#################################################
###      (c) 2015 The TraME Project           ###
###    by The TraME core development team     ###
###      coordinated by Alfred Galichon       ###
###          sources available at:            ###
###       www.alfredgalichon.com/TraME        ###
#################################################
#################################################
#################################################
############  Library estimation:   #############
#####      Matching surplus estimation     ######
#################################################
#################################################


########################################
#' @title Nonparametric Estimation (TU empirical case)
#' 
#' @description Description
#' @param Parameters
#' @export
#' @author Alfred Galichon
#' @examples Examples
NonparametricEsimtationTUEmpirical=function(n,m,hetG,hetH,kron,muhat,theta0,xtol_rel  = 1e-4, maxeval=1e5,print_level=0)
{
  if (print_level>0) {print(paste0("LP optimization used."))}
  nbX = length (n)
  nbY = length (m)
  nbParams = nbX*nbY
  
  
  res1 = build_disaggregate_epsilon(n,nbX,nbY,hetG)
  res2 = build_disaggregate_epsilon(m,nbY,nbX,hetH)
  
  epsilon_iy = res1$epsilon_iy
  epsilon0_i = c(res1$epsilon0_i)
  I_ix = res1$I_ix
  eta_xj = t(res2$epsilon_iy)
  eta0_j = c(res2$epsilon0_i)  
  I_yj = t(res2$I_ix)
  
  ni = c(I_ix %*% n)/res1$nbDraws
  mj = c( m %*% I_yj)/res2$nbDraws
  
  nbI = length(ni)
  nbJ = length(mj)

  A_11 = kronecker(matrix(1,nbY,1),sparseMatrix(1:nbI,1:nbI,x=1))
  A_12 = sparseMatrix(i=NULL,j=NULL,dims=c(nbI*nbY,nbJ),x=0)
  A_13 = kronecker(sparseMatrix(1:nbY,1:nbY,x=-1),I_ix)
  A_14 = sparseMatrix(i=NULL,j=NULL,dims=c(nbI*nbY,nbParams),x=0)
  A_21 = sparseMatrix(i=NULL,j=NULL,dims=c(nbX*nbJ,nbI),x=0)
  A_22 = kronecker(sparseMatrix(1:nbJ,1:nbJ,x=1),matrix(1,nbX,1))
  A_23 = kronecker(t(I_yj),sparseMatrix(1:nbX,1:nbX,x=1))
  A_24 = - kronecker(t(I_yj),sparseMatrix(1:nbX,1:nbX,x=1))
  A_1  = cbind(A_11,A_12,A_13, A_14)
  A_2  = cbind(A_21,A_22,A_23, A_24)
  A    = rbind(A_1,A_2)
  nbconstr = dim(A)[1]
  nbvar = dim(A)[2]
  
  
  lb  =  c(epsilon0_i,t(eta0_j), rep(-Inf,nbX*nbY+nbParams))
  rhs = c(epsilon_iy, eta_xj)
  obj = c(ni,mj,rep(0,nbX*nbY),c(-muhat))
  
  
  gurobiModel = list(A=A,obj=obj,modelsense="min",rhs=rhs,sense=rep(">=",nbconstr),lb=lb)
  result = gurobi ( gurobiModel, params=list(OutputFlag=0) ) 
  
  if (result$status=="OPTIMAL") {
    U = matrix(result$x[(nbI+nbJ+1):(nbI+nbJ+nbX*nbY)],nrow=nbX)
    
    muiy = matrix(result$pi[1:(nbI*nbY)],nrow=nbI)
    phihat = matrix(result$x[(nbI+nbJ+nbX*nbY+1):(nbI+nbJ+nbX*nbY+nbParams)], nbX,nbY)
    V = phihat - U
    mu = t(I_ix) %*% muiy 
    val = result$objval 
  }
  else {
    print("optimization problem with Gurobi")
    return(result)
    }
  
  return(list(phihat = phihat,
              U = U,
              V=V,
              val=val))
  
}

########################################
#' @title Nonparametric Surplus Estimation
#' 
#' @description Description
#' @param Parameters
#' @export
#' @author Alfred Galichon
#' @examples Examples
npe = function(model,muhat, print_level=0)
{
  if (print_level>0) {print(paste0("Moment Matching Estimation of ",class(model)," model."))}
  
  market= parametricMarket(model,initparam(model)$param)
  
  if (class(market$transfers) != "TU") {stop("Nonparametric estimation currently only applies to TU models.")}

  if ((class(market$hetG)=="empirical") & (class(market$hetH)=="empirical"))
  {outcome = NonparametricEsimtationTUEmpirical(market$n,market$m,market$hetG,market$hetH,kron,muhat,theta0,print_level=print_level)}
  else
    {stop("Nonparametric estimation currently defined in empirical case only.")}
  return(outcome)
}

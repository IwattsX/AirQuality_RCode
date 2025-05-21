rm(list = ls())

library(mvtnorm)
library(foreach)
library(doMC)

options("scipen"=10) # show all digits

1# data simulation
ClipSimulation = function(ci, theta, rho, K, Ts, DesignX, seed=NULL){
  
  if(length(ci) != K-2){stop("Number of cut points and categories NOT match!!!")}
  if(length(theta) != NCOL(DesignX)){stop("Number of Design Matrix columns and coefficients NOT match!!!")}
  
  len_par = length(ci) + length(theta) + length(rho)
  
  if(!is.null(seed)){set.seed(seed)}
  mst = DesignX%*%theta#mean
  et = arima.sim(n=Ts, list(ar=rho), sd=sqrt(1-rho^2)) # sd argument is for WN arima
  Z = et + mst
  
  ci = c(0, ci) # c1 = 0
  clip = function(X=NULL){return(length(which(X>ci))+1)}
  X_hour = sapply(X=as.vector(Z), FUN = clip)
  
  X_hour_wide = matrix(0, nrow=Ts, ncol=K)
  for (h in 1:Ts) {X_hour_wide[h,X_hour[h]] = 1}# cator
  
  res = list(Z, X_hour, X_hour_wide)
  names(res) = c("Z", "X_hour", "X_hour_wide")
  
  return(res)
}

################################################################################
############################ Simulation parameter ##############################
################################################################################
Ts = 500 #sample size

# Model: STAT
K = 3
ciT = 0.8615
alpha0T = 0.43075
thetaT = c(alpha0T)
names(thetaT) = "alpha0"
rhoT = 0.2

DesignXT = matrix(1, nrow=Ts, ncol=1) #intercept only 
colnames(DesignXT) = "Intercept"

################################################################################
############################# Simulation Data ##################################
################################################################################
resSim = ClipSimulation(ci=ciT, theta=thetaT, rho=rhoT, K=K, Ts=Ts, DesignX=DesignXT, seed=NULL)
Z = resSim$Z
X_hour = resSim$X_hour#is a cat
X_hour_wide = resSim$X_hour_wide
plot(Z)

# when sample size is large, empirical variance converges to 1/(1-rhoT^2)
1/(1-rhoT^2)
var(Z)


################################################################################
########################## CLS Objective function ##############################
################################################################################
#
AopCLS = function(par, Xl, DesignX){
  
  Ts = length(Xl)
  K = nlevels(factor(Xl))
  
  num_par = length(par)
  
  seg = c(-Inf, 0, par[1:(K-2)], Inf)#
  theta = par[(K-1):(num_par-1)]#
  rho = par[num_par]
  
  # mean vector
  mst = DesignX%*%theta
  CondEX = rep(0, Ts)#
  
  # t = 1 (using marginal expectation)
  tmpprob = pnorm(seg[2:(K+1)]-mst[1]) - pnorm(seg[1:K]-mst[1]) #marginal probablity
  CondEX[1] = sum(tmpprob*(1:K))#summation
  #equation 5
  
  # t > 2 (using conditional expectation) !!! only lag 1 now !!!
  corr_nom = rho^abs(matrix(1:2 - 1, nrow = 2, ncol = 2, byrow = TRUE) - (1:2 - 1))
  diag(corr_nom) = 1
  for(t in 2:Ts){
    prevX = Xl[t-1]
    denom = pnorm(seg[prevX+1]-mst[t-1]) - pnorm(seg[prevX]-mst[t-1])
    mu_nom = mst[t:(t-1)]
    for(k in 1:K){
      a_nom = seg[c(k, prevX)]
      b_nom = seg[c(k, prevX)+1]
      nom = pmvnorm(lower=a_nom, upper=b_nom, mean=mu_nom, corr=corr_nom)
      CondEX[t] = CondEX[t] + k*nom/denom
    }
  }
  
  Q = sum((Xl - CondEX)^2)#equation 4
  
  # cat("\n par =", par)
  return(Q)
}

parT = c(ciT, thetaT, rhoT)
par = parT
Xl = X_hour
DesignX = DesignXT #intercept



parT = c(ciT, thetaT, rhoT)
tim1 = Sys.time()
AopCLS(parT, X_hour, DesignXT)
tim2 = Sys.time()

tim2 - tim1

# parameter estimation
cout = summary(as.factor(X_hour))
ci_initial = as.vector(qnorm(cumsum(cout/sum(cout))))[1:(K-1)]
phi_initial = pacf(X_hour, plot = F)$acf[1]
par_initial = c(ci_initial[2] - ci_initial[1], -ci_initial[1], phi_initial)#parT is true

constrLSE <- matrix(0, nrow=K+1, ncol=length(par_initial))
constrLSE[1,1] = 1
constrLSE[2,1:2] = c(1, -1)
constrLSE[3:4,length(par_initial)] = c(1, -1)
constrLSE_ci = c(0, 0, -1, -1)

OptimFirst <- constrOptim(par_initial, f=AopCLS,#constrain optimization
                          method = "Nelder-Mead", ui=constrLSE,
                          ci=constrLSE_ci, hessian=F,
                          control= list(reltol=1e-07, maxit=1000),
                          Xl=X_hour, DesignX=DesignXT)
OptimFirst




nSim = 100 #500 tries

paramRes = matrix(NA, nrow=nSim, ncol=3)

for(i in 1:nSim){
  
  resSim = ClipSimulation(ci=ciT, theta=thetaT, rho=rhoT, K=K, Ts=Ts, DesignX=DesignXT, seed=1234+i)
  Z = resSim$Z
  X_hour = resSim$X_hour
  X_hour_wide = resSim$X_hour_wide
  
  # parameter estimation
  cout = summary(as.factor(X_hour))
  ci_initial = as.vector(qnorm(cumsum(cout/sum(cout))))[1:(K-1)]
  phi_initial = pacf(X_hour, plot = F)$acf[1]
  par_initial = c(ci_initial[2] - ci_initial[1], -ci_initial[1], phi_initial)
  OptimFirst <- constrOptim(par_initial, f=AopCLS,
                            method = "Nelder-Mead", ui=constrLSE,
                            ci=constrLSE_ci, hessian=F,
                            control= list(reltol=1e-07, maxit=1000),
                            Xl=X_hour, DesignX=DesignXT)
  paramRes[i,] = OptimFirst$par
  cat("\n", i)
}

colMeans(paramRes)
# [1] 1.0494025 0.5200773 0.2225020


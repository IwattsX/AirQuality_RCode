rm(list = ls())

options("scipen"=10) # show all digits
install.packages("doMC")
install.pacakges("covr")
install.packages("Rcpp")
install.packages("RcppArmadillo")
install.packages("readxl")
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doMC)


sourceCpp("CLSTools.cpp")

##### data





library(readxl)



y1t <- read_excel("/home/iwatts/AirQuality_RCode/sortedPhillyPADS.xlsx")
y1t
colnames(y1t)

dates = y1t[ 1:10, "D"]
dates
ClipSimulation = function(ci, theta, rho, K, Ts, DesignX, seed=NULL){
  
  if(length(ci) != K-2){stop("Number of cut points and categories NOT match!!!")}
  if(length(theta) != NCOL(DesignX)){stop("Number of Design Matrix columns and coefficients NOT match!!!")}
  
  len_par = length(ci) + length(theta) + length(rho)
  
  if(!is.null(seed)){set.seed(seed)}
  mst = DesignX%*%theta
  et = arima.sim(n=Ts, list(ar=rho), sd=sqrt(1-rho^2)) # sd argument is for WN
  Z = et + mst
  
  ci2 = c(0, ci) # c1 = 0
  clip = function(X=NULL){return(length(which(X>ci2))+1)}
  X_hour = sapply(X=as.vector(Z), FUN = clip)
  
  X_hour_wide = matrix(0, nrow=Ts, ncol=K)
  for (h in 1:Ts) {X_hour_wide[h,X_hour[h]] = 1}
  
  res = list(Z, X_hour, X_hour_wide)
  names(res) = c("Z", "X_hour", "X_hour_wide")
  
  return(res)
}
### Conditional least square for univariate form by Innovation method objective function


AopCLS = function(par, Xl, DesignX){
  
  Ts = length(Xl)
  K = nlevels(factor(Xl))
  num_par = length(par)
  
  ci = c(0, par[1:(K-2)])
  seg = c(-Inf, ci, Inf)
  theta = par[(K-1):(num_par-1)]
  rho = par[num_par]
  
  ## mean vector
  mst = DesignX%*%theta
  seg.m.1 = matrix(rep(seg[2:(K+1)], each=Ts), nrow=Ts , ncol=K)
  seg.m.1.minusmu = seg.m.1 - matrix(rep(mst, times=K), nrow=Ts, ncol=K)
  seg.m.2 = matrix(rep(seg[1:K], each=Ts), nrow=Ts , ncol=K)
  seg.m.2.minusmu = seg.m.2 - matrix(rep(mst, times=K), nrow=Ts, ncol=K)
  EX = (pnorm(seg.m.1.minusmu) - pnorm(seg.m.2.minusmu)) %*% (1:K)
  
  rm(seg.m.1, seg.m.1.minusmu, seg.m.2, seg.m.2.minusmu)
  
  ##### Innovation Algorithm
  InnovRes = UniInnovRcpp(EXL=EX, mst=mst, seg=seg, phi_est=rho, K=K, numCoef=1500)
  HQ = InnovRes$HQ
  V = InnovRes$V
  mylag = InnovRes$lag
  rm(InnovRes)
  ## One Step Ahead Prediction
  MyPred = ClipPred(EXL=EX, X_hour=Xl, mylag=mylag, HQ=HQ)
  ## conditional expectation
  CondEX = MyPred$Prediction
  
  Q = sum((Xl - CondEX)^2)
  
  # cat("\n par =", par)
  # cat("\n Q =", Q)
  # cat("\n ======== \n")
  
  return(Q)
}



### Single simulation investigation

#```{r code4}
myseed = 12345 + 1

##### Main scripts Start from here
Ts = 500

K = 3
ciT = 0.8615
thetaT = 0.43075
names(thetaT) = c("beta0")
rhoT = 0.5
parT = c(ciT, thetaT, rhoT)

DesignXT = matrix(1, nrow=Ts, ncol=1)
colnames(DesignXT) = "Intercept"

# data simulation
resSim = ClipSimulation(ci=ciT, theta=thetaT, rho=rhoT, K=K, Ts=Ts, DesignX=DesignXT, seed=myseed)
Z = resSim$Z
X_hour = resSim$X_hour
X_hour_wide = resSim$X_hour_wide
#DesignXT = resSim$DesignX

#------------------- Parameter Estimation   ------------------------#
## initial value calculation
cout = summary(as.factor(X_hour))
ci_initial = as.vector(qnorm(cumsum(cout/sum(cout))))[1:(K-1)]
p <- which(!is.finite(ci_initial))
if(length(p)!=0){
  if(p==1){
    ci_initial[p] <- -6
  }else{
    ci_initial[p] <-  6
  }
}
phi_initial = pacf(X_hour, plot = F)$acf[1]
par_initial = c(ci_initial[2:length(ci_initial)] - ci_initial[1], -ci_initial[1], 
                rep(0, NCOL(DesignXT)-1), phi_initial)
## parameter constraints
constrLSE = matrix(0, nrow=K, ncol=length(par_initial))
constrLSE[1,1] = 1
# for(ii in 1:(K-3)){constrLSE[ii+1,c(ii, ii+1)] <- c(-1,1)} # ci
constrLSE[(K-1):K,length(par_initial)] = c(1, -1)
constrLSE_ci = c(rep(0,K-2), -1, -1)

## optimization comparison
OptimLInnov = constrOptim(par_initial, f=AopCLS,
                          method = "Nelder-Mead", ui=constrLSE,
                          ci=constrLSE_ci, hessian=F,
                          control= list(reltol=1e-4),
                          Xl=X_hour, DesignX=DesignXT)
OptimLInnov

OptimLInnov = constrOptim(par_initial, f=AopCLS,
                          method = "Nelder-Mead", ui=constrLSE,
                          ci=constrLSE_ci, hessian=F,
                          control= list(reltol=1e-5),
                          Xl=X_hour, DesignX=DesignXT)
OptimLInnov

OptimLInnov = constrOptim(par_initial, f=AopCLS,
                          method = "Nelder-Mead", ui=constrLSE,
                          ci=constrLSE_ci, hessian=F,
                          control= list(reltol=1e-10),
                          Xl=X_hour, DesignX=DesignXT)
OptimLInnov


#Multiple simulations

nSim = 200
resW = resL = matrix(0, nrow=nSim, ncol=3)
for(iii in 1:nSim){
  
  # data simulation
  resSim = ClipSimulation(ci=ciT, theta=thetaT, rho=rhoT, K=K, Ts=Ts, DesignX=DesignXT, seed=myseed+iii)
  Z = resSim$Z
  X_hour = resSim$X_hour
  X_hour_wide = resSim$X_hour_wide
  
  #------------------- Parameter Estimation   ------------------------#
  myseed = 12345 + iii
  
  ##### Main scripts Start from here
  Ts = 500
  
  K = 3
  ciT = 0.8615
  thetaT = 0.43075
  names(thetaT) = c("beta0")
  rhoT = 0.5
  parT = c(ciT, thetaT, rhoT)
  
  DesignXT = matrix(1, nrow=Ts, ncol=1)
  colnames(DesignXT) = "Intercept"
  
  # data simulation
  resSim = ClipSimulation(ci=ciT, theta=thetaT, rho=rhoT, K=K, Ts=Ts, DesignX=DesignXT, seed=myseed)
  Z = resSim$Z
  X_hour = resSim$X_hour
  X_hour_wide = resSim$X_hour_wide
  # DesignXT = resSim$DesignX
  
  #------------------- Parameter Estimation   ------------------------#
  ## initial value calculation
  cout = summary(as.factor(X_hour))
  ci_initial = as.vector(qnorm(cumsum(cout/sum(cout))))[1:(K-1)]
  p <- which(!is.finite(ci_initial))
  if(length(p)!=0){
    if(p==1){
      ci_initial[p] <- -6
    }else{
      ci_initial[p] <-  6
    }
  }
  phi_initial = pacf(X_hour, plot = F)$acf[1]
  par_initial = c(ci_initial[2:length(ci_initial)] - ci_initial[1], -ci_initial[1], 
                  rep(0, NCOL(DesignXT)-1), phi_initial)
  ## parameter constraints
  constrLSE = matrix(0, nrow=K, ncol=length(par_initial))
  constrLSE[1,1] = 1
  # for(ii in 1:(K-3)){constrLSE[ii+1,c(ii, ii+1)] <- c(-1,1)} # ci
  constrLSE[(K-1):K,length(par_initial)] = c(1, -1)
  constrLSE_ci = c(rep(0,K-2), -1, -1)
  
  ## optimization comparison
  OptimLInnov = constrOptim(par_initial, f=AopCLS,
                            method = "Nelder-Mead", ui=constrLSE,
                            ci=constrLSE_ci, hessian=F,
                            control= list(reltol=1e-04),
                            Xl=X_hour, DesignX=DesignXT)
  resL[iii,] = OptimLInnov$par
}

#Summary of the comparison

mySummary = matrix(NA, nrow=4, ncol=3)
colnames(mySummary) = c("c2", "Intercept", "rho")
rownames(mySummary) = c("TRUE", 
                        "Est", "Bias", "sd")
mySummary[1,] = parT
mySummary[2,] = colMeans(resL)
mySummary[3,] = colMeans(resL) - parT
mySummary[4,] = apply(resL, 2, sd)
mySummary

#Histogram of the comparison

myXLowB = apply(resL, 2, min) - 0.1
myXUppB = apply(resL, 2, max) + 0.1
par(mfrow=c(1,3))
hist(resL[,1], xlim = c(myXLowB[1], myXUppB[1]))
hist(resL[,2], xlim = c(myXLowB[2], myXUppB[2]))
hist(resL[,3], xlim = c(myXLowB[3], myXUppB[3]))

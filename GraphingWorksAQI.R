rm(list = ls())

options("scipen"=10) # show all digits
install.packages("doMC")

library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doMC)
library(readxl)
library(ggplot2)
sourceCpp("CLSTools.cpp")

##### Load Actual Data #####

y1t <- read_excel("C:\\Users\\Ju Wang\\Downloads\\aqi file merger\\sortedCBSAdata\\sortedPhillyPADSB.xlsx", col_names = FALSE)
#y1t <-read_excel("C:\\Users\\Ju Wang\\Downloads\\aqi file merger\\sortedCBSAdata\\sortedPhillyPADS.xlsx", header = TRUE)
colnames(y1t) <- c("date", "AQI")
X1 <- ts(y1t$AQI,frequency = 365)
Y1 <-ts(y1t$date, frequency = 365)


X2<- X1(start=1,end = 100)

ts.plot(X1,
        main = "Time Series Plot of AQI",  # Title
        xlab = "Starting 1/1/2021 and Ending 10/31/2024 ",                    # X-axis label
        ylab = "AQI",
        type ="p")                   # Y-axis label     
plot(as.ts(X1))

#top2_indices <- order(X1, decreasing = TRUE)[1:2]
#print(top2_indices)
# 888 and 889 is June 7 2023 and June 8 2023 because of Canada wildfires

# Load or define the design matrix based on actual data Error in data used
DesignXT = as.matrix(X1)


##### Conditional least square objective function (unchanged)
AopCLS = function(par, Xl, DesignX){
  
  Ts = length(X1)
  K = nlevels(factor(X1))
  num_par = length(par)
  
  Ts = 1400
  K = 4
  num_par = 1
  
  ci = c(0, par[1:2])
  seg = c(-Inf, ci, Inf)
  theta = par[3:(num_par-1)]
  rho = par[num_par]
  
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

  
  ## mean vector
  mst = DesignX %*% theta
  seg.m.1 = matrix(rep(seg[3:(K+1)], each=Ts), nrow=Ts , ncol=K)
  seg.m.1.minusmu = seg.m.1 - matrix(rep(mst, times=K), nrow=Ts, ncol=K)
  seg.m.2 = matrix(rep(seg[2:K], each=Ts), nrow=Ts , ncol=K)
  seg.m.2.minusmu = seg.m.2 - matrix(rep(mst, times=K), nrow=Ts, ncol=K)
  seg.m.3 = matrix(rep(seg[1:K], each=Ts), nrow=Ts , ncol=K)
  seg.m.3.minusmu = seg.m.3 - matrix(rep(mst, times=K), nrow=Ts, ncol=K)
  EX = (pnorm(seg.m.1.minusmu) - pnorm(seg.m.2.minusmu)-pnorm(seq.m.3.minusms)) %*% (1:K)
  rm(seg.m.1, seg.m.1.minusmu, seg.m.2, seg.m.2.minusmu,seg.m.3,seg)
  
  ##### Innovation Algorithm
  InnovRes = UniInnovRcpp(EXL=EX, mst=mst, seg=seg, phi_est = rho, K=K, numCoef=1500)
  HQ = InnovRes$HQ
  V = InnovRes$V
  mylag = InnovRes$lag
  rm(InnovRes)
  
  ## One Step Ahead Prediction
  MyPred = ClipPred(EXL=EX, X_hour=Xl, mylag=mylag, HQ=HQ)
  CondEX = MyPred$Prediction
  
  Q = sum((Xl - CondEX)^2)
  
  return(Q)


##### Parameter Estimation with Actual Data #####
cout = summary(as.factor(X1))
ci_initial = as.vector(qnorm(cumsum(cout/sum(cout))))[1:(K-1)]
p <- which(!is.finite(ci_initial))
if(length(p) != 0){
  if(p == 1){
    ci_initial[p] <- -6
  }else{
    ci_initial[p] <- 6
  }
}
phi_initial = pacf(X1, plot = FALSE)$acf[1]
par_initial = c(ci_initial[2:length(ci_initial)] - ci_initial[1], -ci_initial[1], 
                rep(0, NCOL(DesignXT)-1), phi_initial)

## Parameter constraints
constrLSE = matrix(0, nrow=K, ncol=length(par_initial))
constrLSE[1,1] = 1
for(ii in 1:(K-3)){constrLSE[ii+1,c(ii, ii+1)] <- c(-1,1)} # ci
constrLSE[(K-1):K,length(par_initial)] = c(1, -1)
constrLSE_ci = c(rep(0,K-2), -1, -1)

## Optimization
OptimFirst = constrOptim(par_initial, f=AopCLS,
                         method = "Nelder-Mead", ui=constrLSE,
                         ci=constrLSE_ci, hessian=FALSE,
                         control= list(reltol=1e-04, maxit=1000),
                         X1=X_hour, DesignX=DesignXT)

res = OptimFirst$par

## Save results
save(Z, X_hour, X_hour_wide, DesignXT, ciT, thetaT, rhoT, res, 
     file = "actual_data_results.RData")

print("Parameter estimation completed successfully!")

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






## initial value calculation
K=6
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
for(ii in 1:(K-3)){constrLSE[ii+1,c(ii, ii+1)] <- c(-1,1)} # ci
constrLSE[(K-1):K,length(par_initial)] = c(1, -1)
constrLSE_ci = c(rep(0,K-2), -1, -1)
## optimization
OptimFirst = constrOptim(par_initial, f=AopCLS,
                         method = "Nelder-Mead", ui=constrLSE,
                         ci=constrLSE_ci, hessian=F,
                         control= list(reltol=1e-04, maxit=1000),
                         Xl=X_hour, DesignX=DesignXT)
Param[subid,] = OptimFirst$par
RES$count[subid]  = OptimFirst$counts[1] # Number of function evaluation
RES$convg[subid]  = OptimFirst$convergence

tim.end <- Sys.time()
timeused <- difftime(tim.end, tim.start, units="secs")
cat("\n No.job", job, "task", subid,
    "completed within", timeused, "at seed", myseed, "!")

#################################################################
#---------------------- Result Return --------------------------#
#################################################################
RES$seed[subid]  = myseed
RES$tim[subid]   = timeused

RES = cbind(RES, Param)

return(RES)



useMultiCore <- function(nTasks=NULL, nCore=NULL, iii=NULL,
                         Ts=NULL,
                         K=NULL,
                         ciT=NULL,
                         thetaT=NULL,
                         rhoT=NULL)
{
  
  cat("Multicores working, please wait ... \n")
  registerDoMC(cores = nCore)
  tim.start = Sys.time()
  FinalRes <- foreach(i=1:nCore, .combine = "rbind") %dopar%
    one.job(job=i, nTasks=nTasks, nCore=nCore, iii=iii,
            Ts=Ts, K=K, ciT=ciT, thetaT=thetaT, rhoT=rhoT)
  tim.end = Sys.time()
  cat("Done.\n")
  cat("\n\n nTasks =", nTasks, "\t nCore =", nCore,
      "\t Aveg. time =", mean(FinalRes$tim),
      "\t Total Time =", difftime(tim.end, tim.start, units="hours"),
      "hours \n\n")
  return(FinalRes)
}



rm(list = ls())

options("scipen"=10) # show all digits
install.packages("doMC")

library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doMC)
sourceCpp("CLSTools.cpp")



##### data simulation
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

##### Conditional least square objective function
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
  InnovRes = UniInnovRcpp(EXL=EX, mst=mst, seg=seg, phi_est = rho, K=K, numCoef=1500)
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

##### One job for each core
one.job <- function(job=NULL, nTasks=NULL, nCore=NULL, iii=NULL,
                    Ts=NULL,
                    K=NULL,
                    ciT=NULL,
                    thetaT=NULL,
                    rhoT=NULL)
{
  n.ciT = length(ciT)
  n.thetaT = length(thetaT)
  n.rhoT = length(rhoT)
  n.param = n.ciT + n.thetaT + n.rhoT
  
  nSubtasks = round(nTasks/nCore)
  RES = data.frame(job=job, task=1:nSubtasks, tim=rep(NA, nSubtasks),
                   seed=rep(NA, nSubtasks), count=rep(NA, nSubtasks),
                   convg=rep(NA, nSubtasks))
  Param = matrix(NA, nrow=nSubtasks, ncol=n.param)
  colnames(Param) = c(paste0("ci", 2:(K-1)), names(thetaT), paste0("rho",1:n.rhoT))
  
  for(subid in 1:nSubtasks){
    
    tim.start <- Sys.time()
    #################################################################
    #----------------------- Calculation ---------------------------#
    #################################################################
    myseed = 0 + job*nSubtasks + subid + iii*100000
    # myseed <- 00000 # Complete Random
    
    Intercept = rep(1, Ts)
    x1 = rnorm(n=Ts, mean=-1, sd=1)
    x2 = rnorm(n=Ts, mean=-0.25, sd=sqrt(0.0324))
    DesignXT = cbind(Intercept, x1, x2)
    
    # data simulation
    resSim = ClipSimulation(ci=ciT, theta=thetaT, rho=rhoT, K=K, Ts=Ts, DesignX=DesignXT, seed=myseed)
    Z = resSim$Z
    X_hour = resSim$X_hour
    X_hour_wide = resSim$X_hour_wide
    
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
  }
  
  RES = cbind(RES, Param)
  
  return(RES)
}


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


##### Main scripts Start from here
Ts = 2000

# Model: Varion and Vidoni (2006)
K = 7
ciT = c(1.2, 2.2, 3.1, 4.1, 5.3)
thetaT = c(2.9, -0.6,  9)
names(thetaT) = c("beta0", "beta1", "beta2")
rhoT = 0.5

### Warning: !!!!!!!!!!!
# check your computer available cores before run code below!!!!!!!!!!!!!
detectCores()
### Warning: !!!!!!!!!!!
# You have to set "nCores < number of Cores available in your computer" !!!!!
### Warning: !!!!!!!!!!!
# be careful of your computer memory !!!!!!!!!!!

nCore  = 8    # number of cores you for computation
nTasks = 500   # number of simulations

res = useMultiCore(nTasks=nTasks, nCore=nCore, iii=1, Ts=Ts, 
                   K=K, ciT=ciT, thetaT=thetaT, rhoT=rhoT)
res

save(Ts, ciT, thetaT, rhoT, nCore, nTasks, res, 
     file = paste0("res/CLSE_InnoEst_VarinVidoni2006_", rhoT, ".RData"))
   

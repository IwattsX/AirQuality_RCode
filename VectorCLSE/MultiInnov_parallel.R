rm(list = ls())

options("scipen"=10) # show all digits

library(mvtnorm)
library(foreach)
library(doMC)

multikappa = function(i, j, msti, mstj, EXwi, EXwj, seg, rho){
  
  K = length(seg)-1
  res = matrix(NA, nrow=K-1, ncol=K-1)
  
  # different time
  lag = abs(i-j)
  corr_nom = (rho^lag)^abs(matrix(1:2 - 1, nrow = 2, ncol = 2, byrow = TRUE) - (1:2 - 1))
  for(k in 1:(K-1)){
    for(kp in 1:(K-1)){
      l_nom = seg[c(k, kp)]
      u_nom = seg[c(k, kp)+1]
      meanvec = c(msti, mstj)
      res[k,kp] = pmvnorm(lower=l_nom, upper=u_nom, mean=meanvec, corr=corr_nom)
    }
  }
  res = res - EXwi %*% t(EXwj)
  
  return(res)
}

multiInnov = function(mst, EXW, seg, rho, num_coef=NULL, mytol=1e-06){
  
  Ts = length(mst)
  K = length(seg)-1
  if(is.null(num_coef)){num_coef = Ts}
  
  V = vector(mode="list", length=Ts)
  HQ = array(NA, c(Ts, num_coef, K-1, K-1)) # set dimension to num_coef to save memory for allocation
  
  # t=0 with index 1
  # t=1 with index 2
  
  #--------------------- t=0 ---------------------#
  V[[1]] = multikappa(i=1, j=1, msti=mst[1], mstj=mst[1], EXwi=EXW[1,], EXwj=EXW[1,], seg, rho)
  
  #--------------------- t=1 ---------------------#
  HQ[1,1,,] = multikappa(i=2, j=1, msti=mst[2], mstj=mst[1], EXwi=EXW[2,], EXwj=EXW[1,], seg, rho) %*% solve(V[[1]]) # k=0
  tmp = HQ[1,1,,] %*% V[[1]] %*% t(HQ[1,1,,])
  V[[2]] = multikappa(i=2, j=2, msti=mst[2], mstj=mst[2], EXwi=EXW[2,], EXwj=EXW[2,], seg, rho) - tmp
  
  #--------------------- t>=2 --------------------#
  mylag = 0
  tmpcheck = matrix(100, nrow=K-1, ncol=K-1)
  for(t in 2:(Ts-1)){
    if(any(tmpcheck > mytol) & mylag==0){
      # k = 0
      HQ[t,t,,] = multikappa(i=t+1, j=1, msti=mst[t+1], mstj=mst[1], EXwi=EXW[t+1,], EXwj=EXW[1,], seg, rho) %*% solve(V[[1]])
      # k > 0
      for(k in 1:(t-1)){
        tmp = matrix(0, nrow=K-1, ncol=K-1)
        for(j in 0:(k-1)){
          tmp = tmp + HQ[t,t-j,,] %*% V[[j+1]] %*% t(HQ[k,k-j,,])
        }
        HQ[t,t-k,,] = (multikappa(i=t+1, j=k+1, msti=mst[t+1], mstj=mst[k+1], EXwi=EXW[t+1,], EXwj=EXW[k+1,], seg, rho) - tmp) %*% solve(V[[k+1]])
      }
      # Prediction MSE
      tmp = matrix(0, nrow=K-1, ncol=K-1)
      for(j in 0:(t-1)){
        tmp = tmp + HQ[t,t-j,,] %*% V[[j+1]] %*% t(HQ[t,t-j,,])
      }
      V[[t+1]] = multikappa(i=t+1, j=t+1, msti=mst[t+1], mstj=mst[t+1], EXwi=EXW[t+1,], EXwj=EXW[t+1,], seg, rho) - tmp
      
      tmpcheck = abs(HQ[t-1,t-1,,])
      # set mylag = t-1 if the innovation coefficient small enough
      if(all(abs(HQ[t,t,,]) < mytol)){
        # set lags to save computation time
        mylag = t;
        # erase the rest lags to save memory since they are all 0
        HQ = HQ[,1:mylag,,]
      }
    }else{
      # Only to cover the effective coefficients (lag), the rest are less than 1e-05
      for(k in (t-mylag):(t-1)){
        tmp = matrix(0, nrow=K-1, ncol=K-1)
        for(j in max(k-mylag,0):(k-1)){
          if(t-j<mylag){
            tmp = tmp + HQ[t,t-j,,] %*% V[[j+1]] %*% t(HQ[k,k-j,,])
          }
        }
        HQ[t,t-k,,] = (multikappa(i=t+1, j=k+1, msti=mst[t+1], mstj=mst[k+1], EXwi=EXW[t+1,], EXwj=EXW[k+1,], seg, rho) - tmp) %*% solve(V[[k+1]])
      }
      # Prediction MSE
      tmp = matrix(0, nrow=K-1, ncol=K-1)
      for(j in (t-mylag):(t-1)){
        if(t-j<mylag){
          tmp = tmp + HQ[t,t-j,,] %*% V[[j+1]] %*% t(HQ[t,t-j,,])
        }
      }
      V[[t+1]] = multikappa(i=t+1, j=t+1, msti=mst[t+1], mstj=mst[t+1], EXwi=EXW[t+1,], EXwj=EXW[t+1,], seg, rho) - tmp
      tmpcheck = abs(HQ[t-1,mylag,,])
    }
  }
  
  RES = list(HQ=HQ, V=V, mylag=mylag)
  return(RES)
}

multiInnovALL = function(mst, EXW, seg, rho){
  
  Ts = length(mst)
  K = length(seg)-1
  
  V = vector(mode="list", length=Ts)
  HQ = array(NA, c(Ts, Ts, K-1, K-1))
  
  # t=0 with index 1
  # t=1 with index 2
  
  #--------------------- t=0 ---------------------#
  V[[1]] = multikappa(i=1, j=1, msti=mst[1], mstj=mst[1], EXwi=EXW[1,], EXwj=EXW[1,], seg, rho)
  
  #--------------------- t=1 ---------------------#
  HQ[1,1,,] = multikappa(i=2, j=1, msti=mst[2], mstj=mst[1], EXwi=EXW[2,], EXwj=EXW[1,], seg, rho) %*% solve(V[[1]]) # k=0
  tmp = HQ[1,1,,] %*% V[[1]] %*% t(HQ[1,1,,])
  V[[2]] = multikappa(i=2, j=2, msti=mst[2], mstj=mst[2], EXwi=EXW[2,], EXwj=EXW[2,], seg, rho) - tmp
  
  #--------------------- t>=2 --------------------#
  mylag = 0
  tmpcheck = matrix(100, nrow=K-1, ncol=K-1)
  for(t in 2:(Ts-1)){
    # k = 0
    HQ[t,t,,] = multikappa(i=t+1, j=1, msti=mst[t+1], mstj=mst[1], EXwi=EXW[t+1,], EXwj=EXW[1,], seg, rho) %*% solve(V[[1]])
    # k > 0
    for(k in 1:(t-1)){
      tmp = matrix(0, nrow=K-1, ncol=K-1)
      for(j in 0:(k-1)){
        tmp = tmp + HQ[t,t-j,,] %*% V[[j+1]] %*% t(HQ[k,k-j,,])
      }
      HQ[t,t-k,,] = (multikappa(i=t+1, j=k+1, msti=mst[t+1], mstj=mst[k+1], EXwi=EXW[t+1,], EXwj=EXW[k+1,], seg, rho) - tmp) %*% solve(V[[k+1]])
    }
    # Prediction MSE
    tmp = matrix(0, nrow=K-1, ncol=K-1)
    for(j in 0:(t-1)){
      tmp = tmp + HQ[t,t-j,,] %*% V[[j+1]] %*% t(HQ[t,t-j,,])
    }
    V[[t+1]] = multikappa(i=t+1, j=t+1, msti=mst[t+1], mstj=mst[t+1], EXwi=EXW[t+1,], EXwj=EXW[t+1,], seg, rho) - tmp
  }
  
  RES = list(HQ=HQ, V=V, mylag=mylag)
  return(RES)
}

MultiOneStepPred = function(EXW, Xw, mylag, HQ){
  
  Err = Xw - EXW
  
  PredErr = matrix(0, nrow=NROW(Err), ncol=NCOL(Err))
  for(t in 1:(Ts-1)){
    
    if(mylag > 0){
      predlag = mylag
    }else{
      predlag = t
    }
    
    for(j in 1:min(t,predlag)){
      PredErr[t+1,] = PredErr[t+1,] + HQ[t,j,,]%*%(Err[t+1-j,]-PredErr[t+1-j,])
    }
    
  }
  
  ErrErr = Err - PredErr
  PredX = PredErr + EXW
  
  return(list(ErrErr=ErrErr, PredX=PredX))
}

CLSWInnov = function(par, Xl, Xw, DesignX, UseAllLag=FALSE, num_coef=NULL){
  
  Ts = length(Xl)
  K = nlevels(factor(Xl))
  num_par = length(par)
  
  ci = c(0, par[1:(K-2)])
  seg = c(-Inf, ci, Inf)
  theta = par[(K-1):(num_par-1)]
  rho = par[num_par]
  
  ## mean vector
  mst = DesignX%*%theta
  
  ## Conditional expectation
  tmp = t(matrix(rep(seg, Ts), nrow=length(seg), ncol=Ts)) - matrix(rep(mst, length(seg)), nrow=Ts)
  EXW = pnorm(tmp[,2:(K+1)]) - pnorm(tmp[,1:K])
  EXW = EXW[,1:(K-1)]
  
  # Innovation function starts from here
  if(UseAllLag){
    ## all lags version
    InnovRes = multiInnovALL(mst=mst, EXW=EXW, seg=seg, rho=rho)
    HQ = InnovRes$HQ
    V = InnovRes$V
    mylag=InnovRes$mylag
  }else{
    ## necessary lags version
    InnovRes = multiInnov(mst=mst, EXW=EXW, seg=seg, rho=rho, num_coef=NULL, mytol=1e-06)
    HQ = InnovRes$HQ
    V = InnovRes$V
    mylag=InnovRes$mylag
  }
  
  Predw = MultiOneStepPred(EXW, Xw[,1:(K-1)], mylag, HQ)
  Predw = cbind(Predw$PredX, 1-rowSums(Predw$PredX))
  
  res = sum((Xw - Predw)^2)
  
  # cat("\n=====================")
  # cat("\n par = ", par)
  # cat("\n res = ", res)
  
  return(res)
}

AopCLSW = function(par, Xl, Xw, DesignX, lag=1){
  
  Ts = length(Xl)
  K = NCOL(Xw)
  
  num_par = length(par)
  
  seg = c(-Inf, 0, par[1:(K-2)], Inf)
  theta = par[(K-1):(num_par-1)]
  rho = par[num_par]
  
  # mean vector
  mst = DesignX%*%theta
  CondEX = matrix(0, nrow=Ts, ncol=K)
  
  # t = 1 (using marginal expectation)
  tmpprob = pnorm(seg[2:(K+1)]-mst[1]) - pnorm(seg[1:K]-mst[1])
  CondEX[1,] = tmpprob
  
  # t = 2 (using conditonal expectation with lag=1)
  corr_nom = rho^abs(matrix(1:2 - 1, nrow = 2, ncol = 2, byrow = TRUE) - (1:2 - 1))
  diag(corr_nom) = 1
  prevX = Xl[1]
  denom = pnorm(seg[prevX+1]-mst[1]) - pnorm(seg[prevX]-mst[1])
  mu_nom = mst[2:1]
  for(k in 1:K){
    a_nom = seg[c(k, prevX)]
    b_nom = seg[c(k, prevX)+1]
    nom = pmvnorm(lower=a_nom, upper=b_nom, mean=mu_nom, corr=corr_nom)
    CondEX[2,k] = nom/denom
  }
  
  # t > 2 (using conditional expectation with lag)
  for(t in 3:Ts){
    if(t<=lag){
      # 2 < t <= lag
      corr_nom = rho^abs(matrix(1:t-1, nrow=t, ncol=t, byrow=TRUE) - (1:t - 1))
      diag(corr_nom) = 1
      corr_denom = rho^abs(matrix(1:(t-1)-1, nrow=t-1, ncol=t-1, byrow = TRUE) - (1:(t-1) - 1))
      diag(corr_denom) = 1
      
      prevID = 1:(t-1)
      prevX = Xl[prevID]
      a_denom = seg[prevX]
      b_denom = seg[prevX+1]
      mu_denom = mst[prevID]
      denom = pmvnorm(lower=a_denom, upper=b_denom, mean=mu_denom, corr=corr_denom)
      
      currID = 1:t
      mu_nom = mst[currID]
      
      for(k in 1:K){
        a_nom = seg[c(k, prevX)]
        b_nom = seg[c(k, prevX)+1]
        nom = pmvnorm(lower=a_nom, upper=b_nom, mean=mu_nom, corr=corr_nom)
        CondEX[t,k] = nom/denom
      }
    }else{
      # t > lag (using conditional expectation)
      corr_nom = rho^abs(matrix(1:(lag+1) - 1, nrow = (lag+1), ncol = (lag+1), byrow = TRUE) - (1:(lag+1) - 1))
      diag(corr_nom) = 1
      corr_denom = rho^abs(matrix(1:lag - 1, nrow = lag, ncol = lag, byrow = TRUE) - (1:lag - 1))
      diag(corr_denom) = 1
      prevID = (t-1):(t-lag)
      prevX = Xl[prevID]
      a_denom = seg[prevX]
      b_denom = seg[prevX+1]
      mu_denom = mst[prevID]
      if(lag<2){
        denom = pnorm(b_denom-mu_denom) - pnorm(a_denom-mu_denom)
      }else{
        denom = pmvnorm(lower=a_denom, upper=b_denom, mean=mu_denom, corr=corr_denom)  
      }
      
      currID = t:(t-lag)
      mu_nom = mst[currID]
      for(k in 1:K){
        a_nom = seg[c(k, prevX)]
        b_nom = seg[c(k, prevX)+1]
        nom = pmvnorm(lower=a_nom, upper=b_nom, mean=mu_nom, corr=corr_nom)
        CondEX[t,k] = nom/denom
      }
    }
  }
  
  Q = sum((Xw - CondEX)^2)
  
  # cat("\n ===============")
  # cat("\n par =", par)
  # cat("\n Q =", Q)
  return(Q)
}

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
  
  res = list(Z, X_hour, X_hour_wide, DesignX)
  names(res) = c("Z", "X_hour", "X_hour_wide", "DesignX")
  
  return(res)
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
    
    DesignXT = matrix(1, nrow=Ts, ncol=1)
    colnames(DesignXT) = "Intercept"
    
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
    # for(ii in 1:(K-3)){constrLSE[ii+1,c(ii, ii+1)] <- c(-1,1)} # ci
    constrLSE[(K-1):K,length(par_initial)] = c(1, -1)
    constrLSE_ci = c(rep(0,K-2), -1, -1)
    ## optimization
    OptimWInnov = constrOptim(par_initial, f=CLSWInnov,
                              method = "Nelder-Mead", ui=constrLSE,
                              ci=constrLSE_ci, hessian=F,
                              Xl=X_hour, Xw=X_hour_wide, DesignX=DesignXT)
    Param[subid,] = OptimWInnov$par
    RES$count[subid]  = OptimWInnov$counts[1] # Number of function evaluation
    RES$convg[subid]  = OptimWInnov$convergence
    
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
Ts = 500

# Model
K = 3
ciT = 0.8615
thetaT = 0.43075
names(thetaT) = c("beta0")

### Warning: !!!!!!!!!!!
# check your computer available cores before run code below!!!!!!!!!!!!!
detectCores()
### Warning: !!!!!!!!!!!
# You have to set "nCores < number of Cores available in your computer" !!!!!
### Warning: !!!!!!!!!!!
# be careful of your computer memory !!!!!!!!!!!

nCore  = 20     # number of cores you for computation
nTasks = 1000   # number of simulations

rhoT = 0.5

res = useMultiCore(nTasks=nTasks, nCore=nCore, iii=1, Ts=Ts, 
                   K=K, ciT=ciT, thetaT=thetaT, rhoT=rhoT)
res

save(Ts, ciT, thetaT, rhoT, nCore, nTasks, res, 
     file = paste0("MultiInnov_rho", rhoT, "_iii1.RData"))

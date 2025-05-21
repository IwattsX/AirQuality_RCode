### Set Environment ###
options(scipen = 10) # show more digits
# install.packages("doMC") # Only need to install once

### Load Libraries ###
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(foreach)
library(doMC)
library(readxl)
library(ggplot2)

### Load C++ Functions ###
sourceCpp("CLSTools.cpp")

### Load and Prepare Data ###
y1t <- read_excel("C:/Users/Ju Wang/Downloads/aqi file merger/sortedCBSAdata/sortedPhillyPADSB.xlsx", col_names = FALSE)
colnames(y1t) <- c("date", "AQI")
X1 <- ts(y1t$AQI, frequency = 365)
Y1 <- ts(y1t$date, frequency = 365)

# Plot AQI Time Series
ts.plot(X1, main = "Time Series Plot of AQI", xlab = "Time (2021-2024)", ylab = "AQI", type = "p")
plot(as.ts(X1))

# Design Matrix
DesignXT <- as.matrix(X1)

##### Conditional Least Squares Objective Function #####
AopCLS <- function(par, Xl, DesignX) {
  Ts <- length(Xl)
  K <- nlevels(factor(Xl))
  
  # Parameter decomposition
  ci <- c(0, par[1:(K-1)])
  seg <- c(-Inf, ci, Inf)
  theta <- par[K:(length(par)-1)]
  rho <- par[length(par)]
  
  # Mean + innovations
  mst <- DesignX %*% theta
  et <- arima.sim(n = Ts, list(ar = rho), sd = sqrt(1 - rho^2))
  Z <- mst + et
  
  # Mapping to categories
  ci2 <- c(0, ci)
  clip <- function(x) sum(x > ci2) + 1
  X_hour <- sapply(Z, clip)
  
  X_hour_wide <- matrix(0, nrow = Ts, ncol = K)
  for (h in 1:Ts) X_hour_wide[h, X_hour[h]] <- 1
  
  # Innovation Algorithm
  mst <- DesignX %*% theta
  seg_matrices <- function(seg, mst, K, Ts) {
    sapply(1:K, function(k) pnorm(seg[k+1] - mst) - pnorm(seg[k] - mst))
  }
  EX <- seg_matrices(seg, mst, K, Ts) %*% (1:K)
  
  InnovRes <- UniInnovRcpp(EXL = EX, mst = mst, seg = seg, phi_est = rho, K = K, numCoef = 1500)
  HQ <- InnovRes$HQ
  V <- InnovRes$V
  mylag <- InnovRes$lag
  
  # Prediction and Loss
  MyPred <- ClipPred(EXL = EX, X_hour = Xl, mylag = mylag, HQ = HQ)
  CondEX <- MyPred$Prediction
  
  Q <- sum((Xl - CondEX)^2)
  return(Q)
}

##### Estimation Preparation #####
initialize_parameters <- function(X, K, DesignX) {
  cout <- summary(as.factor(X))
  ci_initial <- qnorm(cumsum(cout/sum(cout)))[1:(K-1)]
  
  p <- which(!is.finite(ci_initial))
  if (length(p) > 0) {
    ci_initial[p] <- ifelse(p == 1, -6, 6)
  }
  
  phi_initial <- pacf(X, plot = FALSE)$acf[1]
  par_initial <- c(ci_initial[2:length(ci_initial)] - ci_initial[1], -ci_initial[1], rep(0, NCOL(DesignX)-1), phi_initial)
  
  # Constraints for optimization
  constrLSE <- matrix(0, nrow = K, ncol = length(par_initial))
  constrLSE[1,1] <- 1
  for (ii in 1:(K-3)) constrLSE[ii+1, c(ii, ii+1)] <- c(-1,1)
  constrLSE[(K-1):K, length(par_initial)] <- c(1, -1)
  constrLSE_ci <- c(rep(0, K-2), -1, -1)
  
  list(par_initial = par_initial, constrLSE = constrLSE, constrLSE_ci = constrLSE_ci)
}

##### Run Optimization #####
run_optimization <- function(X, DesignX, K) {
  init <- initialize_parameters(X, K, DesignX)
  
  OptimResult <- constrOptim(
    par = init$par_initial, 
    f = AopCLS, 
    method = "Nelder-Mead",
    ui = init$constrLSE, 
    ci = init$constrLSE_ci, 
    control = list(reltol = 1e-5, maxit = 1000),
    Xl = X, 
    DesignX = DesignX
  )
  
  return(OptimResult)
}

##### Parameter Estimation #####
K <- 4  # Categories (AQI levels)
OptimResult <- run_optimization(X1, DesignXT, K)

# Save the result
save(OptimResult, file = "actual_data_results.RData")
print("Parameter estimation completed successfully!")

##### Parallel Estimation Function #####
one.job <- function(job, nTasks, nCore, iii, Ts, K, ciT, thetaT, rhoT) {
  RES <- data.frame(count = NA, convg = NA, seed = NA, tim = NA)
  Param <- matrix(NA, nrow = 1, ncol = length(thetaT)+length(ciT)+1)
  
  tim.start <- Sys.time()
  
  # optimization
  OptimResult <- run_optimization(X1, DesignXT, K)
  Param[1,] <- OptimResult$par
  RES$count <- OptimResult$counts[1]
  RES$convg <- OptimResult$convergence
  
  tim.end <- Sys.time()
  RES$tim <- difftime(tim.end, tim.start, units="secs")
  
  RES$seed <- iii  # assuming `iii` is random seed
  
  return(cbind(RES, Param))
}

##### Multicore Setup #####
useMultiCore <- function(nTasks, nCore, iii, Ts, K, ciT, thetaT, rhoT) {
  cat("Multicores working, please wait ...\n")
  registerDoMC(cores = nCore)
  tim.start <- Sys.time()
  
  FinalRes <- foreach(i = 1:nCore, .combine = "rbind") %dopar% 
    one.job(i, nTasks, nCore, iii, Ts, K, ciT, thetaT, rhoT)
  
  tim.end <- Sys.time()
  cat("Done.\n")
  cat("\n\nnTasks =", nTasks, "\tnCore =", nCore, 
      "\tAverage Time =", mean(FinalRes$tim),
      "\tTotal Time =", difftime(tim.end, tim.start, units = "hours"), "hours\n\n")
  
  return(FinalRes)
}


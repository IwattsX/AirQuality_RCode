library(readxl)

data_raw <- read_excel("C:\\Users\\Ju Wang\\Downloads\\aqi file merger\\sortedCBSAdata\\sortedPhillyPADSB.xlsx", col_names = FALSE)
colnames(data_raw) <- c("date", "AQI")

AQI_ts <- ts(data_raw$AQI, frequency=365)
date_ts <- ts(data_raw$date, frequency=365)

# For example, make DesignXT
Intercept = rep(1, length(AQI_ts))
DesignXT = as.matrix(Intercept) # simple design matrix with just intercept for now
one.job <- function(job=NULL, nTasks=NULL, nCore=NULL, iii=NULL,
                    Ts=NULL,
                    K=NULL,
                    ciT=NULL,
                    thetaT=NULL,
                    rhoT=NULL,
                    X_real=NULL,    # add this argument for real X data
                    DesignX_real=NULL) { # add this too
  
  nSubtasks = round(nTasks/nCore)
  RES = data.frame(job=job, task=1:nSubtasks, tim=rep(NA, nSubtasks),
                   seed=rep(NA, nSubtasks), count=rep(NA, nSubtasks),
                   convg=rep(NA, nSubtasks))
  
  Param = matrix(NA, nrow=nSubtasks, ncol=length(ciT)+length(thetaT)+length(rhoT))
  colnames(Param) = c(paste0("ci", 2:(K-1)), names(thetaT), paste0("rho",1))
  
  for(subid in 1:nSubtasks){
    
    tim.start <- Sys.time()
    
    ## Use real data instead of simulation
    X_hour <- as.integer(factor(X_real))
    
    ## initial values
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
    par_initial = c(ci_initial[2:length(ci_initial)] - ci_initial[1], -ci_initial[1], rep(0, ncol(DesignX_real)-1), phi_initial)
    
    constrLSE = matrix(0, nrow=K, ncol=length(par_initial))
    constrLSE[1,1] = 1
    for(ii in 1:(K-3)){constrLSE[ii+1,c(ii, ii+1)] <- c(-1,1)}
    constrLSE[(K-1):K,length(par_initial)] = c(1, -1)
    constrLSE_ci = c(rep(0,K-2), -1, -1)
    
    ## optimization
    OptimFirst = constrOptim(par_initial, f=AopCLS,
                             method = "Nelder-Mead", ui=constrLSE,
                             ci=constrLSE_ci, hessian=F,
                             control= list(reltol=1e-04, maxit=1000),
                             Xl=X_hour, DesignX=DesignX_real)
    
    Param[subid,] = OptimFirst$par
    RES$count[subid]  = OptimFirst$counts[1]
    RES$convg[subid]  = OptimFirst$convergence
    
    tim.end <- Sys.time()
    timeused <- difftime(tim.end, tim.start, units="secs")
    
    cat("\n No.job", job, "task", subid,
        "completed within", timeused)
    
    RES$tim[subid]   = timeused
  }
  
  RES = cbind(RES, Param)
  return(RES)
}

K <- nlevels(factor(xl))
ci_initial <- seq(0.5, K-1.5, by=1) # example initial cutpoints
theta_initial <- rep(0, ncol(DesignX)) # starting theta
rho_initial <- 0.5

par_initial <- c(ci_initial, theta_initial, rho_initial)

# Constraint matrices
constrLSE <- matrix(0, nrow=K, ncol=length(par_initial))
constrLSE[1,1] <- 1
for(ii in 1:(K-3)){ constrLSE[ii+1, c(ii, ii+1)] <- c(-1,1) }
constrLSE[(K-1):K, length(par_initial)] <- c(1, -1)
constrLSE_ci <- c(rep(0,K-2), -1, -1)

# Optimization
result <- constrOptim(par=par_initial, f=AopCLS,
                      method="Nelder-Mead", 
                      ui=constrLSE, ci=constrLSE_ci,
                      hessian=FALSE,
                      Xl=Xl, DesignX=DesignX)

result$par

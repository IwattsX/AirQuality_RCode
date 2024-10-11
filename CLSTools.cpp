// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(mvtnorm)]]
#include <mvtnormAPI.h>

using namespace Rcpp;
using namespace arma;

int len_par, len_mean;
int Ts, Th, Td, K, t, i, j, k, p, n;

// [[Rcpp::export]]
double pmvnorm_cpp(arma::vec& lb, arma::vec& ub, arma::vec& mu, arma::vec& lowertrivec, double abseps = 1e-3){
  
  int n = lb.n_elem;
  int nu = 0;
  int maxpts = 25000;     // default in mvtnorm: 25000
  double releps = 0;      // default in mvtnorm: 0
  int rnd = 1;            // Get/PutRNGstate
  
  double* lb_ = lb.memptr();                        // lower bound;
  double* ub_ = ub.memptr();                        // upper bound;
  double* correlationMatrix = lowertrivec.memptr(); // array of correlation coefficients;
  int* infin = new int[n];                          // Integer, array of integration limits flag;
  double* mu_ = mu.memptr();                        // array of non-centrality parameters;
  
  // if INFIN(I) < 0, Ith limits are (-infinity, infinity);
  // if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
  // if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
  // if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)];
  
  for (int i = 0; i < n; ++i) {
    if(lb(i) == R_NegInf){
      infin[i] = 0;
    }else if (ub(i) == R_PosInf){
      infin[i] = 1;
    }else{
      infin[i] = 2;
    }
  }
  
  // return values
  double error;
  double value;
  int inform;
  
  // Fortran function details check:
  // https://github.com/cran/mvtnorm/blob/master/src/mvt.f
  mvtnorm_C_mvtdst(&n, &nu, lb_, ub_,
                   infin, correlationMatrix, mu_,
                   &maxpts, &abseps, &releps,
                   &error, &value, &inform, &rnd);
                   delete[] (infin);
                   
                   return value;
}

// [[Rcpp::export]]
double KappLongRcpp(int ti, int tj, double mui, double muj, double EXi, double EXj, int K, vec seg, double phi){
  
  double res=0.0, prob=0.0;
  vec mymu(2), lb(2), ub(2), correle(1);
  mymu = {mui, muj};
  correle(0) = pow(phi, abs(ti-tj));
  
  for(i=0;i<K;i++){
    for(j=0;j<K;j++){
      lb = {seg(i), seg(j)};
      ub = {seg(i+1), seg(j+1)};
      prob = pmvnorm_cpp(lb, ub, mymu, correle, 1e-5);
      res += (i+1)*(j+1)*prob;
    }
  }
  return(res - EXi*EXj);
}


// [[Rcpp::export]]
List UniInnovRcpp(vec EXL, vec mst, vec seg, double phi_est, int K, int numCoef){
  
  Ts = EXL.size();
  int nn;
  
  double tmp;
  vec V(Ts, fill::zeros);
  mat HQ(Ts, numCoef, fill::zeros);
  
  //---------- n = 0 ----------//
  V(0) = KappLongRcpp(0, 0, mst(0), mst(0), EXL(0), EXL(0), K, seg, phi_est);
  
  //---------- n = 1 ----------//
  HQ(0,0) = KappLongRcpp(1, 0, mst(1), mst(0), EXL(1), EXL(0), K, seg, phi_est)/V(0);
  tmp = V(0)*pow(HQ(0,0), 2.0);
  V(1) = KappLongRcpp(1, 1, mst(1), mst(1), EXL(1), EXL(1), K, seg, phi_est) - tmp;
  
  //---------- n >= 2 ----------// n=nn+1
  int mylag=0;
  double tmpcheck=100;
  for(nn=1;nn<Ts-1;nn++){
    //for(nn=1;nn<100;nn++){
    if((tmpcheck>0.00001) & (mylag==0)){
      // k=0
      HQ(nn,nn) = KappLongRcpp(nn+1, 0, mst(nn+1), mst(0), EXL(nn+1), EXL(0), K, seg, phi_est)/V(0);
      // k>0
      for(k=1;k<nn+1;k++){
        tmp = 0.0;
        for(j=0;j<k;j++){
          tmp += HQ(k-1,k-1-j)*HQ(nn,nn-j)*V(j);
        }
        HQ(nn,nn-k) = (KappLongRcpp(nn+1+1, k+1, mst(nn+1), mst(k), EXL(nn+1), EXL(k), K, seg, phi_est)-tmp)/V(k);
      }
      // Prediction MSE //
      tmp = 0.0;
      for(j=0;j<nn+1;j++){
        tmp += HQ(nn,nn-j)*HQ(nn,nn-j)*V(j);
      }
      V(nn+1) = KappLongRcpp(nn+1, nn+1, mst(nn+1), mst(nn+1), EXL(nn+1), EXL(nn+1), K, seg, phi_est) - tmp;
      
      tmpcheck = fabs(HQ(nn-1,nn-1));
      // set the lag=n-1 if the innovation coefficient small enough //
      if(fabs(HQ(nn,nn)) < 0.00001){
        // set lags to save computation time
        mylag = nn;
        // erase therest lags to save memory since they are all 0
        HQ.shed_cols(mylag+1, numCoef-1);
      }
    }else{
      // Only to cover the effective coefficients (lag), the rest are less than 1e-05
      for(k=nn+1-mylag;k<nn+1;k++){
        tmp = 0.0;
        for(j=std::max(k-mylag,0);j<k;j++){
          if(nn-j<mylag+1){
            tmp += HQ(k-1,k-1-j)*HQ(nn,nn-j)*V(j);
          }
        }
        HQ(nn,nn-k) = (KappLongRcpp(nn+1+1, k+1, mst(nn+1), mst(k), EXL(nn+1), EXL(k), K, seg, phi_est)-tmp)/V(k);
      }
      // Prediction MSE //
      tmp = 0.0;
      for(j=nn+1-mylag;j<nn+1;j++){
        if(nn-j<mylag+1){
          tmp += HQ(nn,nn-j)*HQ(nn,nn-j)*V(j);
        }
      }
      V(nn+1) = KappLongRcpp(nn+1, nn+1, mst(nn+1), mst(nn+1), EXL(nn+1), EXL(nn+1), K, seg, phi_est) - tmp;
      
      tmpcheck = fabs(HQ(nn-1,mylag));
    }
  }
  
  List res;
  res["HQ"]=HQ;
  res["V"]=V;
  res["lag"]=mylag;
  
  return(res);
}


// [[Rcpp::export]]
List ClipPred(vec EXL, vec X_hour, int mylag, mat HQ){
  
  Ts = EXL.size();
  
  int predlag;
  vec Err = X_hour - EXL;
  vec PredX(Ts, fill::zeros);
  vec PredErr(Ts, fill::zeros), ErrErr(Ts, fill::zeros);
  
  for(n=0;n<Ts-1;n++){
    if(mylag>0){
      predlag = mylag;
    }else{
      predlag = n;
    }
    for(j=0;j<std::min(n+1,predlag);j++){
      PredErr(n+1) += HQ(n,j)*(Err(n+1-j-1)-PredErr(n+1-j-1));
    }
  }
  ErrErr = Err - PredErr;
  PredX = PredErr + EXL;
  
  //---------- Return Result ----------//
  List res;
  res["Innovation"] = ErrErr;
  res["Prediction"] = PredX;
  return(res);
}
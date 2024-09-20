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
double LSE_cpp1(arma::vec par, arma::mat X, arma::mat DesignX){
  
  Th = X.n_rows;
  K = X.n_cols;
  len_par = par.size();
  
  //--------------------------- Declar identifiers ----------------------//
  double h;
  double phi, beta0, beta1, beta2, mst;
  arma::vec ci(K-1, fill::zeros);
  arma::vec seg(K+1, fill::zeros);
  arma::vec norm_quan;
  arma::vec tmp(Th, fill::zeros);
  arma::mat EX(Th, K, fill::zeros);
  
  //-------------------------- Parameter Extraction ---------------------//
  ci(0)   = 0;
  ci.subvec(1,K-2) = par.subvec(0,K-3);
  seg(0)  = R_NegInf;
  seg.subvec(1, K-1) = ci;
  seg(K)  = R_PosInf;
  beta0   = par(K-2);
  beta1   = par(K-1);
  beta2   = par(K);
  phi     = par(K+1);
  
  for(t=0;t<Th;t++){
    tmp(t) = pow(phi, t);
  }
  double tmp1, tmp2, tmp3;
  for(h=0;h<Th;h++){
    tmp1 = 0;
    tmp2 = 0;
    tmp3 = 0;
    for(t=0;t<h;t++){
      tmp1 += tmp(t);
      tmp2 += tmp(t)*DesignX(h-t,1);
      tmp3 += tmp(t)*DesignX(h-t,2);
    }
    mst = beta0*tmp1 + beta1*tmp2 + beta2*tmp3;
    // norm_quan = normcdf((seg-mst)/pow(1-phi*phi, 0.5));
    norm_quan = normcdf((seg-mst)/(1/pow(1-phi*phi, 0.5)));
    EX.row(h) = trans(norm_quan.subvec(1,K) - norm_quan.subvec(0,K-1));
  }
  
  //-------------------------    Results Output   ----------------------//
  // Rcout << "par = " << par << std::endl;
  double LSE = accu(pow(X-EX, 2.0));
  return(LSE);
}



// [[Rcpp::export]]
double LSE_cpp2(arma::vec par, arma::vec X, int K, arma::mat DesignX){
  
  Th = X.size();
  len_par = par.size();
  
  //--------------------------- Declar identifiers ----------------------//
  double h;
  double phi, beta0, beta1, beta2;
  arma::vec ci(K-1, fill::zeros);
  arma::vec seg(K+1, fill::zeros);
  arma::vec norm_quan(K+1, fill::zeros);
  arma::vec denorm_quan(K+1, fill::zeros);
  arma::vec mst(Th, fill::zeros);
  arma::vec tmp(Th, fill::zeros);
  arma::vec CLEX(Th, fill::zeros);
  
  //-------------------------- Parameter Extraction ---------------------//
  ci(0)   = 0;
  ci.subvec(1,K-2) = par.subvec(0,K-3);
  seg(0)  = R_NegInf;
  seg.subvec(1, K-1) = ci;
  seg(K)  = R_PosInf;
  beta0   = par(K-2);
  beta1   = par(K-1);
  beta2   = par(K);
  phi     = par(K+1);
  
  for(t=0;t<Th;t++){
    tmp(t) = pow(phi, t);
  }
  double tmp1, tmp2, tmp3;
  for(h=0;h<Th;h++){
    tmp1 = 0;
    tmp2 = 0;
    tmp3 = 0;
    for(t=0;t<h;t++){
      tmp1 += tmp(t);
      tmp2 += tmp(t)*DesignX(h-t,1);
      tmp3 += tmp(t)*DesignX(h-t,2);
    }
    mst(h) = beta0*tmp1 + beta1*tmp2 + beta2*tmp3;
  }
  
  double tmpphi = 1/pow(1-phi*phi, 0.5);
  
  // t=1
  norm_quan = normcdf((seg-mst(0))/tmpphi);
  for(k=1;k<(K+1);k++){
    CLEX(0) = CLEX(0) + (k+1)*(norm_quan(k) - norm_quan(k-1));
  }
  
  // t>1
  arma::mat tmpsigma = { {1/(1-phi*phi), phi/(1-phi*phi)}, {phi/(1-phi*phi), 1/(1-phi*phi)} };
  Environment pkg = Environment::namespace_env("mvtnorm");
  Function f = pkg["pmvnorm"];
  
  double kp;
  arma::vec denorm;
  NumericVector lb(2), ub(2), mymu(2);
  
  // t>1
  for(t=1;t<Th;t++){
    denorm_quan = normcdf((seg-mst(t-1))/tmpphi);
    denorm = denorm_quan.subvec(1,K) - denorm_quan.subvec(0,K-1);
    // Rcout << "denorm_quan = " << denorm_quan << std::endl;
    kp = X(t-1)-1;
    for(k=0;k<K;k++){
      lb = {seg(kp), seg(k)};
      ub = {seg(kp+1), seg(k+1)};
      mymu = {mst(t-1), mst(t)};
      CLEX(t) = CLEX(t) + (k+1)*as<double>(f(lb, ub, mymu, R_NilValue, tmpsigma))/denorm(kp);
    }
  }
  
  //-------------------------    Results Output   ----------------------//
  // Rcout << "par = " << trans(par) << std::endl;
  double CLSE = accu(pow(X-CLEX, 2.0));
  return(CLSE);
}


// [[Rcpp::export]]
double LSE_cpp3(arma::vec par, arma::vec Xl, arma::mat X, arma::mat DesignX){
  
  Th = X.n_rows;
  K = X.n_cols;
  len_par = par.size();
  
  //--------------------------- Declar identifiers ----------------------//
  double h;
  double phi, beta0, beta1, beta2;
  arma::vec ci(K-1, fill::zeros);
  arma::vec seg(K+1, fill::zeros);
  arma::vec norm_quan(K+1, fill::zeros);
  arma::vec denorm_quan(K+1, fill::zeros);
  arma::vec mst(Th, fill::zeros);
  arma::vec tmp(Th, fill::zeros);
  arma::mat CLEX(Th, K, fill::zeros);
  
  //-------------------------- Parameter Extraction ---------------------//
  ci(0)   = 0;
  ci.subvec(1,K-2) = par.subvec(0,K-3);
  seg(0)  = R_NegInf;
  seg.subvec(1, K-1) = ci;
  seg(K)  = R_PosInf;
  beta0   = par(K-2);
  beta1   = par(K-1);
  beta2   = par(K);
  phi     = par(K+1);
  
  for(t=0;t<Th;t++){
    tmp(t) = pow(phi, t);
  }
  double tmp1, tmp2, tmp3;
  for(h=0;h<Th;h++){
    tmp1 = 0;
    tmp2 = 0;
    tmp3 = 0;
    for(t=0;t<h;t++){
      tmp1 += tmp(t);
      tmp2 += tmp(t)*DesignX(h-t,1);
      tmp3 += tmp(t)*DesignX(h-t,2);
    }
    mst(h) = beta0*tmp1 + beta1*tmp2 + beta2*tmp3;
  }
  
  double tmpphi = 1/pow(1-phi*phi, 0.5);
  
  // t=1
  norm_quan = normcdf((seg-mst(0))/tmpphi);
  CLEX.row(0) = trans(norm_quan.subvec(1,K) - norm_quan.subvec(0,K-1));
  
  // t>1
  arma::mat tmpsigma = { {1/(1-phi*phi), phi/(1-phi*phi)}, {phi/(1-phi*phi), 1/(1-phi*phi)} };
  Environment pkg = Environment::namespace_env("mvtnorm");
  Function f = pkg["pmvnorm"];
  
  double kp;
  arma::vec denorm;
  NumericVector lb(2), ub(2), mymu(2);
  // t>1
  for(t=1;t<Th;t++){
    kp = Xl(t-1)-1;
    denorm_quan = normcdf((seg-mst(t-1))/tmpphi);
    denorm = denorm_quan.subvec(1,K) - denorm_quan.subvec(0,K-1);
    // Rcout << "denorm(kp) = " << denorm(kp) << std::endl;
    for(k=0;k<K;k++){
      lb = {seg(kp), seg(k)};
      ub = {seg(kp+1), seg(k+1)};
      mymu = {mst(t-1), mst(t)};
      // Rcout << "kp = "<< kp << std::endl;
      // Rcout << "lb = "<< lb << std::endl;
      // Rcout << "ub = "<< ub << std::endl;
      // Rcout << "===="<< pmvnorm_cpp(lb, ub, mymu, correle, 1e-3) << std::endl;
      // CLEX(t,k) = pmvnorm_cpp(lb, ub, mymu, correle, 1e-3)/denorm(kp);
      CLEX(t,k) = as<double>(f(lb, ub, mymu, R_NilValue, tmpsigma))/denorm(kp);
    }
  }
  
  //-------------------------    Results Output   ----------------------//
  // Rcout << "par = " << trans(par) << std::endl;
  double CLSE = accu(pow(X-CLEX, 2.0));
  return(CLSE);
}
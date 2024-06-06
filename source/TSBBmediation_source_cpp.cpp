// #define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>
// #include <vector>
// #include <cmath>
// #include <Rmath.h>
//#include "utilFunctions.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace Rcpp;

template <typename T>
inline bool approx_equal_cpp(const T& lhs, const T& rhs, double tol = 1e-10) {
  return arma::approx_equal(lhs, rhs, "absdiff", tol);
}

umat unique_rows(const umat& M) {
  int n_temp = M.n_rows;
  uvec ind_temp = arma::zeros<arma::uvec>(n_temp);
  
  for (int i = 0; i < n_temp; i++) {
    for (int j = i + 1; j < n_temp; j++) {
      if (approx_equal_cpp(M.row(i), M.row(j))) { 
        ind_temp(j) = 1; 
        break; 
      }
    }
  }
  
  ind_temp = find(ind_temp == 0);
  umat unique_M = M.rows(ind_temp);
  
  uvec unique_M1 = unique_M.col(0);
  uvec unique_M2 = unique_M.col(1);
  
  uvec sorted_unique_M1 = sort(unique_M1);
  uvec sorted_unique_M2 = sort(unique_M2);
  
  uvec M1_temp = unique(sorted_unique_M1);
  int n_temp1 = M1_temp.size();
  int count = 0;
  for (int i = 0; i < n_temp1; i++) {
    ind_temp = find(unique_M1 == M1_temp(i));
    int n_temp2 = ind_temp.size();
    uvec unique_M2_temp = sort(unique_M2(ind_temp));
    
    for (int j = 0; j < n_temp2; j++) {
      sorted_unique_M2(count) = unique_M2_temp(j);
      count++;
    }
  }
  
  umat sorted_unique_M = join_rows(sorted_unique_M1, sorted_unique_M2);
  
  return sorted_unique_M;
}

int mod(int const& a, int const& b){
  return a - floor(a/b)*b;
}

// [[Rcpp::export]]
mat inv_cpp(mat const& MAT) {
  int p_temp = MAT.n_cols;
  int rank_temp = rank(MAT);
  
  mat MATinv;
  if(p_temp > rank_temp){
    MATinv = pinv(MAT);
  } else if(MAT.is_sympd()){
    MATinv = inv_sympd(MAT);
  } else {
    MATinv = inv(MAT);
  }
  
  return MATinv;
}

// [[Rcpp::export]]
double det_cpp(mat const& MAT, bool const& logt) {
  // int p_temp = MAT.n_cols;
  // int rank_temp = rank(MAT);
  // 
  // double detMAT;
  // if(p_temp > rank_temp){
  //   detMAT = 0;
  // } else if(MAT.is_sympd()){
  //   detMAT = log_det_sympd(MAT);
  // } else {
  //   detMAT = real(log_det(MAT));
  // }
  
  double detMAT = real(log_det(MAT));
  
  if(logt == false) {detMAT = exp(detMAT);}
  
  return detMAT;
}

// [[Rcpp::export]]
int rank_cpp(mat const& MAT) {
  return rank(MAT);
}

// [[Rcpp::export]]
vec invlogit_cpp(vec const& x_temp) { 
  return((1/(1+exp(-x_temp))));
}

// [[Rcpp::export]]
vec count_cpp(vec const& VEC, int const& p) {
  vec count(p);
  for(int q = 0; q < p; q++) {
    uvec index = find(VEC == q);
    count(q) = index.n_elem; 
  }
  
  return(count);
}

// [[Rcpp::export]]
double mvbeta_cpp(vec const& VEC, bool const& logt) {
  double mvbeta = sum(lgamma(VEC))-lgamma(sum(VEC));
  if(logt == false) {mvbeta = exp(mvbeta);}
  
  return(mvbeta);
}

// [[Rcpp::export]]
vec factorial_cpp(vec const& N, bool const& logt) {
  int p = N.size();
  
  vec factorial = zeros(p);
  for(int q = 0; q < p; q++) {
    int N_temp = N(q);
    if(N_temp > 0){
      for(int i = 1; i <= N_temp; i++) {
        factorial(q) += log(i);
      }
    }
  }
  
  if(logt == false) {factorial = exp(factorial);}
  
  return(factorial);
}

int rmultinom_cpp(vec const& prob) {
  int p = prob.size();
  
  vec cumprob = cumsum(prob/sum(prob));
  
  int rmultinom = 1;
  double rnd = randu(1)[0];
  for(int q = 0; q < p; q++) {
    if(rnd>cumprob(q)) {
      rmultinom = rmultinom+1;
    }
  }
  
  return(rmultinom);
}

int dmultinom_cpp(vec const& X, vec const& prob, bool const& logt) {
  vec N(1);
  N(0) = X.size();
  vec log_prob = log(prob);
  log_prob.replace(-datum::inf,0);
  
  double dmultinom = factorial_cpp(N,true)[0] - sum(factorial_cpp(X,true) + X % log_prob);
  if(logt == false) {dmultinom = exp(dmultinom);}
  
  return(dmultinom);
}

// [[Rcpp::export]]
mat rdirichlet_cpp(int const& N, vec const& alpha) {
  int p = alpha.size();
  
  mat rdirichlet(N,p);
  for(int q = 0; q < p; q++) {
    if (alpha(q) == 0) {
      rdirichlet.col(q) = zeros(1);
    } else {
      rdirichlet.col(q) = randg(N,distr_param(alpha(q),1.0));
    }
  }
  vec sum_rdirichlet = sum(rdirichlet,1);
  rdirichlet.each_col() /= sum_rdirichlet;
  
  return rdirichlet.t();
}

// [[Rcpp::export]]
double ddirichlet_cpp(vec const& X, vec const& alpha, bool const& logt) {
  int p = alpha.size();
  vec ones_p = ones(p);
  
  double ddirichlet = lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - ones_p) % log(X));
  if(logt == false) {ddirichlet = exp(ddirichlet);}
  
  return ddirichlet;
}

// [[Rcpp::export]]
vec rscainvchisq_cpp(int const& N, double const& a_tau2, double const& b_tau2) {
  vec ones_N = ones(N);
  vec ab_tau_N = a_tau2 * b_tau2 * ones_N;
  
  vec rtau = ab_tau_N/chi2rnd(a_tau2,N);
  return rtau;
}

// [[Rcpp::export]]
double dscainvchisq_cpp(double const& tau, double const& a_tau2, double const& b_tau2, bool const& logt) {
  double ab_tau = a_tau2 * b_tau2;
  double a_tau_half = a_tau2/2;
  double ab_tau_half = ab_tau/2;
  
  double dtau = a_tau_half * log(ab_tau_half) - ab_tau_half/tau - 
    lgamma(a_tau_half) - (1+a_tau_half) * log(tau);
  if(logt == false) {dtau = exp(dtau);}
  
  return dtau;
}

// [[Rcpp::export]]
mat rmvn_cpp(int const& N, vec const& MU, mat const& SIG) {
  return mvnrnd(MU, SIG, N);
}

// [[Rcpp::export]]
double dmvn_cpp(vec const& X_temp, vec const& MU, mat const& SIG, bool const& logt) {
  int p = SIG.n_cols;
  mat Chol_Sig = chol(SIG);
  vec std_temp = solve(Chol_Sig, X_temp - MU); //check
  double std2_temp = sum(std_temp % std_temp);
  
  double dmvn = - (sum(log(Chol_Sig.diag())) + p*log(2*M_PI)/2 + std2_temp/2);
  if(logt == false) {dmvn = exp(dmvn);}
  
  return dmvn;
}

/* FOR CONTINUOUS OUTCOME */
List update_reg_con_cpp(double const& Y_temp, mat const& matX_temp, int const& p_temp,
                        double const& a_var, double const& b_var,
                        double const& a_diff, double const& b_diff,
                        vec const& a_coef, mat const& B_coef,
                        mat const& Binv_coef, mat const& aBinv_coef) {
  
  // Set initial value
  double diff_curr = rscainvchisq_cpp(1, a_diff, b_diff)[0];
  vec coef_curr = mvnrnd(a_coef, diff_curr*B_coef);
  
  // Update variance
  vec var_diff_temp = Y_temp - (matX_temp * coef_curr);
  double var_sum2_temp = dot(var_diff_temp,var_diff_temp);
  
  double a_var_new = (a_var + 1);
  double b_var_new = (a_var * b_var + var_sum2_temp)/a_var_new;
  double var_prop = rscainvchisq_cpp(1, a_var_new, b_var_new)[0];
  
  // Update diffusion
  double diff_sum2_temp = as_scalar((coef_curr - a_coef).t() * Binv_coef * (coef_curr - a_coef));
  
  double a_diff_new = (a_diff + p_temp);
  double b_diff_new = (a_diff * b_diff + diff_sum2_temp)/a_diff_new;
  double diff_prop = rscainvchisq_cpp(1, a_diff_new, b_diff_new)[0];
  
  // Update coefficient
  double constant = diff_prop/var_prop;
  mat Binv_coef_new = constant * (matX_temp.t() * matX_temp) + Binv_coef;
  mat B_coef_new = inv_cpp(Binv_coef_new);
  vec a_coef_new = B_coef_new * (constant * (matX_temp.t() * Y_temp) + aBinv_coef);
  vec coef_prop = mvnrnd(a_coef_new, diff_prop*B_coef_new);
  
  return List::create(_["var_par"]  = var_prop,
                      _["coef_par"] = coef_prop);
}

/* FOR BINARY OUTCOME */
vec update_reg_bin_cpp(int const& Y_temp, vec const& matX_temp, vec const& coef_curr,
                       vec const& a_coef, mat const& B_coef) {
  vec coef_prop = mvnrnd(coef_curr, B_coef);
  
  //careful this is a sum in original
  double loglik_prop = invlogit_cpp(matX_temp*coef_prop)[0] + 
    dmvn_cpp(coef_prop,a_coef,B_coef,true);
  double loglik_curr = invlogit_cpp(matX_temp*coef_curr)[0] + 
    dmvn_cpp(coef_curr,a_coef,B_coef,true);
  
  double ratio = exp(loglik_prop-loglik_curr);
  if(ratio < randu(1)[0]){ coef_prop = coef_curr; }
  
  return coef_prop;
}


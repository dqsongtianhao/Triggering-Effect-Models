#include <RcppArmadillo.h>
using namespace arma;  
// [[Rcpp::depends(RcppArmadillo)]]


//////// data generation
/*** R
library(dplyr)
# cumulative self-triggering hazard, without terms of covariates (time-invariant during a period)
# requires a vector of alpha, logbeta, tlag with the same length m; integral over t1 to t2; lambda0 is baseline
Lambda <- function(t1, t2, alpha, beta, tlag, lambda0) integrate(Vectorize(function(t) exp(sum(alpha * exp(-beta * (t - tlag)), na.rm = T)) * lambda0(t)), t1, t2)$value
# inverse of cumulative self-triggering hazard
# returns a value t2 such that integral of Lambda over t1 to t2 equals x (not exceeding tmax)
Lambda_inv <- function(t1, x, alpha, beta, tlag, lambda0, tmax) {
  if(Lambda(t1, tmax, alpha, beta, tlag, lambda0) < x) return(tmax)
  uniroot(function(t2) Lambda(t1, t2, alpha, beta, tlag, lambda0) - x, c(t1, tmax))$root
}

gen0 <- function(
    id, gamma, alpha, beta,
    lambda0, start_func, censor_func, p_drop, tmax) {
  x = rbinom(length(gamma), 1, 0.5)
  xr = sum(x * gamma)
  t = start_func()
  c = censor_func()
  if(t >= c) return(NULL)
  tlag = rep(NA, length(alpha))
  n = 0
  data = c(id, t, 0, 1, 0, x)
  
  while(1) {
    t0 = t
    t = Lambda_inv(t0, -log(runif(1)) * exp(-xr), alpha, beta, tlag, lambda0, tmax)
    if(t == t0) next
    if(c <= t) {
      data = rbind(data, c(id, c, 0, 0, 0, x))
      break
    }
    else {
      n = n + 1
      if(runif(1) < p_drop) {
        data = rbind(data, c(id, t, 1, 0, n, x))
        break
      }
      else {
        data = rbind(data, c(id, t, 1, 1, n, x))
        tlag = c(t, tlag[-length(alpha)])
      }
    }
  }
  colnames(data) = c("id", "time", "delta", "inout", "n", paste0("X", 1:length(gamma)))
  rownames(data) = NULL
  data.frame(data)
}
gen <- function(
    n, p, m, gamma, alpha, beta,
    lambda0 = function(t) 1,
    start_func = function() 0,
    censor_func = function() runif(1, 1, 4),
    p_drop = 0.1,
    tmax = 4) {
  if(length(gamma) == 1) gamma = rep(gamma, p)
  if(length(alpha) == 1) alpha = rep(alpha, m)
  if(length(beta) == 1) beta = rep(beta, m)
  bind_rows(lapply(1:n, function(id) gen0(id, gamma, alpha, beta, lambda0, start_func, censor_func, p_drop, tmax)))
}
*/

/*** R
pl <- function(par, cov.list, m, data) {
  p <- length(cov.list)
  gamma <- par[1:p]
  alpha <- par[(p+1):(p+m)]
  beta <- par[(p+m+1):(p+m*2)]
  data <- data %>%
    mutate(id = id %>% factor() %>% as.integer() - 1)
  X <- data %>%
    select(all_of(cov.list)) %>%
    as.matrix()
  with(data, pl_c(gamma, alpha, beta, id, time, delta, inout, X))
}
*/
// [[Rcpp::export]]
double pl_c(vec gamma, vec alpha, vec beta,
            vec id, vec time, vec delta, vec inout, mat X) {
  int n = max(id) + 1;
  int p = gamma.n_elem;
  int m = alpha.n_elem;
  // Sort all observations by time.
  uvec idx = stable_sort_index(time);
  
  double pl = 0;
  double t = 0;
  mat Z = zeros(n, p); // fixed-effect covariate matrix
  mat Z0 = Z;
  mat H = zeros(n, m); // history event time
  mat H0 = H;
  mat HI = H; // history event type
  mat HI0 = H;
  vec Y = zeros(n); // at-risk status
  vec Y0 = Y;
  for(uword i : idx) {
    if(time(i) != t) {
      Z = Z0;
      H = H0;
      HI = HI0;
      Y = Y0;
      t = time(i);
    }
    if(delta(i) == 1) {
      double sum_exp = 0;
      for(int j = 0; j < n; j++) {
        if(Y(j) == 1) {
          sum_exp += exp(as_scalar(Z.row(j) * gamma + HI.row(j) % alpha.t() * exp(-beta % (time(i) - H.row(j).t()))));
        }
      }
      pl += as_scalar(Z.row(id(i)) * gamma + HI.row(id(i)) % alpha.t() * exp(-beta % (time(i) - H.row(id(i)).t()))) - log(sum_exp);
      if(m > 1) {
        H0(id(i), span(1, m - 1)) = H0(id(i), span(0, m - 2));
        HI0(id(i), span(1, m - 1)) = HI0(id(i), span(0, m - 2));
      }
      H0(id(i), 0) = time(i);
      HI0(id(i), 0) = 1;
    }
    Y0(id(i)) = inout(i);
    Z0.row(id(i)) = X.row(i);
  }
  return pl;
}

/*** R
gr <- function(par, cov.list, m, data) {
  p <- length(cov.list)
  gamma <- par[1:p]
  alpha <- par[(p+1):(p+m)]
  beta <- par[(p+m+1):(p+m*2)]
  data <- data %>%
    mutate(id = id %>% factor() %>% as.integer() - 1)
  X <- data %>%
    select(all_of(cov.list)) %>%
    as.matrix()
  with(data, gr_c(gamma, alpha, beta, id, time, delta, inout, X))
}
*/
// [[Rcpp::export]]
vec gr_c(vec gamma, vec alpha, vec beta,
         vec id, vec time, vec delta, vec inout, mat X) {
  int n = max(id) + 1;
  int p = gamma.n_elem;
  int m = alpha.n_elem;
  // Sort all observations by time.
  uvec idx = stable_sort_index(time);
  
  double t = 0;
  mat Z = zeros(n, p);
  mat Z0 = Z;
  mat H = zeros(n, m);
  mat H0 = H;
  mat HI = H;
  mat HI0 = H;
  vec Y = zeros(n);
  vec Y0 = Y;
  
  vec gr = zeros(p + m + m);
  for(uword i : idx) {
    if(time(i) != t) {
      Z = Z0;
      H = H0;
      HI = HI0;
      Y = Y0;
      t = time(i);
    }
    if(delta(i) == 1) {
      double sum_exp = 0;
      vec sum_exp_Z = zeros(p);
      vec sum_exp_beta = zeros(m);
      vec sum_exp_alpha_beta = zeros(m);
      
      for(int j = 0; j < n; j++) {
        if(Y(j) == 1) {
          double agg0 = exp(as_scalar(Z.row(j) * gamma + HI.row(j) % alpha.t() * exp(-beta % (time(i) - H.row(j).t()))));
          sum_exp += agg0;
          sum_exp_Z += agg0 * Z.row(j).t();
          sum_exp_beta += agg0 * HI.row(j).t() % exp(-beta % (time(i) - H.row(j).t()));
          sum_exp_alpha_beta -= agg0 * HI.row(j).t() % alpha % (time(i) - H.row(j).t()) % exp(-beta % (time(i) - H.row(j).t()));
        }
      }
      gr(span(0, p - 1)) += Z.row(id(i)).t() - sum_exp_Z / sum_exp;
      gr(span(p, p + m - 1)) += HI.row(id(i)).t() % exp(-beta % (time(i) - H.row(id(i)).t())) - sum_exp_beta / sum_exp;
      gr(span(p + m, p + m + m - 1)) -= HI.row(id(i)).t() % alpha % (time(i) - H.row(id(i)).t()) % exp(-beta % (time(i) - H.row(id(i)).t())) + sum_exp_alpha_beta / sum_exp;
      
      if(m > 1) {
        H0(id(i), span(1, m - 1)) = H0(id(i), span(0, m - 2));
        HI0(id(i), span(1, m - 1)) = HI0(id(i), span(0, m - 2));
      }
      H0(id(i), 0) = time(i);
      HI0(id(i), 0) = 1;
    }
    Y0(id(i)) = inout(i);
    Z0.row(id(i)) = X.row(i);
  }
  return gr;
}

/*** R
hs <- function(par, cov.list, m, data) {
  p <- length(cov.list)
  gamma <- par[1:p]
  alpha <- par[(p+1):(p+m)]
  beta <- par[(p+m+1):(p+m*2)]
  data <- data %>%
    mutate(id = id %>% factor() %>% as.integer() - 1)
  X <- data %>%
    select(all_of(cov.list)) %>%
    as.matrix()
  with(data, hs_c(gamma, alpha, beta, id, time, delta, inout, X))
}
*/
// [[Rcpp::export]]
mat hs_c(vec gamma, vec alpha, vec beta,
         vec id, vec time, vec delta, vec inout, mat X) {
  int n = max(id) + 1;
  int p = gamma.n_elem;
  int m = alpha.n_elem;
  // Sort all observations by time.
  uvec idx = stable_sort_index(time);
  
  double t = 0;
  mat Z = zeros(n, p);
  mat Z0 = Z;
  mat H = zeros(n, m);
  mat H0 = H;
  mat HI = H;
  mat HI0 = H;
  vec Y = zeros(n);
  vec Y0 = Y;
  
  mat hs = zeros(p + m + m, p + m + m);
  for(uword i : idx) {
    if(time(i) != t) {
      Z = Z0;
      H = H0;
      HI = HI0;
      Y = Y0;
      t = time(i);
    }
    if(delta(i) == 1) {
      // define phi as gamma*z + u + alpha * exp(-beta * (t - tlag))
      double sum_exp = 0;
      vec sum_exp_Z = zeros(p);
      mat sum_exp_ZZ = zeros(p, p);
      mat sum_exp_beta2 = zeros(m, m);
      mat sum_exp_alpha_beta2 = zeros(m, m);
      vec sum_exp_beta = zeros(m);
      vec sum_exp_alpha_beta = zeros(m);
      mat sum_exp_Z_beta = zeros(p, m);
      mat sum_exp_Z_alpha_beta = zeros(p, m);
      mat sum_exp_alpha_beta_beta = zeros(m, m);
      
      for(int j = 0; j < n; j++) {
        if(Y(j) == 1) {
          double agg0 = exp(as_scalar(Z.row(j) * gamma + HI.row(j) % alpha.t() * exp(-beta % (time(i) - H.row(j).t()))));
          sum_exp += agg0;                                                      // exp(phi)
          sum_exp_Z += agg0 * Z.row(j).t();                                     // d[exp(phi)]/d[Z]
          sum_exp_ZZ += agg0 * Z.row(j).t() * Z.row(j);                         // d^2[exp(phi)]/d[Z]d[Z']
          
          vec exp_beta = agg0 * HI.row(j).t() % exp(-beta % (time(i) - H.row(j).t()));
          vec exp_alpha_beta = -alpha % (time(i) - H.row(j).t()) % exp_beta;
          sum_exp_beta2 += exp_beta * exp_beta.t() / agg0;                      // d^2[exp(phi)]/d[alpha]d[alpha']
          sum_exp_alpha_beta2 += diagmat(-(time(i) - H.row(j).t()) % exp_alpha_beta) +
            exp_alpha_beta * exp_alpha_beta.t() / agg0;                         // d^2[exp(phi)]/d[beta]d[beta']
          sum_exp_beta += exp_beta;                                             // d[exp(phi)]/d[alpha]
          sum_exp_alpha_beta += exp_alpha_beta;                                 // d[exp(phi)]/d[beta]
          
          mat exp_Z_beta = Z.row(j).t() * exp_beta.t();
          mat exp_Z_alpha_beta = Z.row(j).t() * exp_alpha_beta.t();
          sum_exp_Z_beta += exp_Z_beta;                                         // d^2[exp(phi)]/d[Z]d[alpha]
          sum_exp_Z_alpha_beta += exp_Z_alpha_beta;                             // d^2[exp(phi)]/d[Z]d[beta]
          
          sum_exp_alpha_beta_beta += diagmat(exp_alpha_beta / alpha) + exp_beta * exp_alpha_beta.t() / agg0; // d^2[exp(phi)]/d[alpha]d[beta]
        }
      }
      hs(span(0, p - 1), span(0, p - 1)) += (sum_exp_Z * sum_exp_Z.t() / sum_exp - sum_exp_ZZ) / sum_exp;
      hs(span(0, p - 1), span(p, p + m - 1)) += (sum_exp_Z * sum_exp_beta.t() / sum_exp - sum_exp_Z_beta) / sum_exp;
      hs(span(0, p - 1), span(p + m, p + m + m - 1)) += (sum_exp_Z * sum_exp_alpha_beta.t() / sum_exp - sum_exp_Z_alpha_beta) / sum_exp;
      hs(span(p, p + m - 1), span(p, p + m - 1)) += (sum_exp_beta * sum_exp_beta.t() / sum_exp - sum_exp_beta2) / sum_exp;
      hs(span(p, p + m - 1), span(p + m, p + m + m - 1)) += -diagmat(HI.row(id(i)).t() % (time(i) - H.row(id(i)).t()) % exp(-beta % (time(i) - H.row(id(i)).t()))) +
        (sum_exp_beta * sum_exp_alpha_beta.t() / sum_exp - sum_exp_alpha_beta_beta) / sum_exp;
      hs(span(p + m, p + m + m - 1), span(p + m, p + m + m - 1)) += diagmat(alpha % HI.row(id(i)).t() % (time(i) - H.row(id(i)).t()) % (time(i) - H.row(id(i)).t()) % exp(-beta % (time(i) - H.row(id(i)).t()))) +
        (sum_exp_alpha_beta * sum_exp_alpha_beta.t() / sum_exp - sum_exp_alpha_beta2) / sum_exp;
      hs = symmatu(hs);
      
      if(m > 1) {
        H0(id(i), span(1, m - 1)) = H0(id(i), span(0, m - 2));
        HI0(id(i), span(1, m - 1)) = HI0(id(i), span(0, m - 2));
      }
      H0(id(i), 0) = time(i);
      HI0(id(i), 0) = 1;
    }
    Y0(id(i)) = inout(i);
    Z0.row(id(i)) = X.row(i);
  }
  return hs;
}

// Cumulative baseline hazard
/*** R
cH <- function(par, cov.list, m, data) {
  p <- length(cov.list)
  gamma <- par[1:p]
  alpha <- par[(p+1):(p+m)]
  beta <- par[(p+m+1):(p+m*2)]
  data <- data %>%
    mutate(id = id %>% factor() %>% as.integer() - 1)
  X <- data %>%
    select(all_of(cov.list)) %>%
    as.matrix()
  with(data, cH_c(gamma, alpha, beta, id, time, delta, inout, X))
}
*/
// [[Rcpp::export]]
mat cH_c(vec gamma, vec alpha, vec beta,
            vec id, vec time, vec delta, vec inout, mat X) {
  int n = max(id) + 1;
  int p = gamma.n_elem;
  int m = alpha.n_elem;
  // Sort all observations by time.
  uvec idx = stable_sort_index(time);
  
  double t = 0;
  mat Z = zeros(n, p); // fixed-effect covariate matrix
  mat Z0 = Z;
  mat H = zeros(n, m); // history event time
  mat H0 = H;
  mat HI = H; // history event type
  mat HI0 = H;
  vec Y = zeros(n); // at-risk status
  vec Y0 = Y;
  
  mat cH = zeros(2, sum(delta));
  double cH0 = 0;
  int pos = 0;
  for(uword i : idx) {
    if(time(i) != t) {
      Z = Z0;
      H = H0;
      HI = HI0;
      Y = Y0;
      t = time(i);
    }
    if(delta(i) == 1) {
      double sum_exp = 0;
      for(int j = 0; j < n; j++) {
        if(Y(j) == 1) {
          sum_exp += exp(as_scalar(Z.row(j) * gamma + HI.row(j) % alpha.t() * exp(-beta % (time(i) - H.row(j).t()))));
        }
      }
      cH0 += 1 / sum_exp;
      cH(0, pos) = time(i);
      cH(1, pos) = cH0;
      if(m > 1) {
        H0(id(i), span(1, m - 1)) = H0(id(i), span(0, m - 2));
        HI0(id(i), span(1, m - 1)) = HI0(id(i), span(0, m - 2));
      }
      H0(id(i), 0) = time(i);
      HI0(id(i), 0) = 1;
      pos += 1;
    }
    Y0(id(i)) = inout(i);
    Z0.row(id(i)) = X.row(i);
  }
  return cH;
}
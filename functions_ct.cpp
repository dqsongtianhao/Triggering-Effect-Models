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

event <- function(start, type, tlag, elag, x, gamma_mat, alpha_mat, beta_mat, m_vec, lambda0_list, tmax) {
  m = m_vec[type]
  if(length(elag) > m) {
    tlag = tlag[1:m]
    elag = elag[1:m]
  }
  Lambda_inv(start,
             -log(runif(1)) * exp(-sum(x * gamma_mat[, type])),
             alpha_mat[elag, type], beta_mat[elag, type], tlag, lambda0_list[[type]], tmax)
}
gen0 <- function(
    id,
    m_vec, # a vector of m's: with length k
    gamma_mat, # a matrix of p*k gamma's, gamma[, j] are covariates for type j
    alpha_mat, # a matrix of k*k alpha's, alpha[i, j] is effect from ith type to jth type
    beta_mat, # a matrix of k*k beta's, beta[i, j] is effect from ith type to jth type
    lambda0_list, # a list of functions of baseline hazard, with length k
    start_func, censor_func, p_drop, tmax) {
  p = dim(gamma_mat)[1]
  k = dim(gamma_mat)[2]
  x = rbinom(p, 1, 0.5)
  t = start_func()
  c = censor_func()
  tlag = c()
  elag = c()
  data = c(id, t - 1, 0, 1, x)
  
  while(1) {
    t0 = t
    t_events = sapply(1:k, event,
                      start = t, tlag = tlag, elag = elag, x = x,
                      gamma_mat = gamma_mat, alpha_mat = alpha_mat, beta_mat = beta_mat,
                      m_vec = m_vec, lambda0_list = lambda0_list, tmax = tmax)
    t = min(t_events)
    e = which.min(t_events)
    if(t == t0) next
    if(c <= t) {
      data = rbind(data, c(id, c, 0, 0, x))
      break
    }
    else {
      if(runif(1) < p_drop) {
        data = rbind(data, c(id, t, e, 0, x))
        break
      }
      else {
        data = rbind(data, c(id, t, e, 1, x))
        tlag = c(t, tlag)
        elag = c(e, elag)
      }
    }
  }
  colnames(data) = c("id", "time", "delta", "inout", paste0("X", 1:p))
  rownames(data) = NULL
  data.frame(data)
}
gen <- function(
    n, p, k, m_vec, gamma_mat, alpha_mat, beta_mat,
    lambda0_list = list(function(t) 1),
    start_func = function() 0,
    censor_func = function() runif(1, 1, 4),
    p_drop = 0.1,
    tmax = 4) {
  if(length(m_vec) == 1) m_vec = rep(m_vec, k)
  if(length(gamma_mat) == 1) gamma_mat = matrix(gamma_mat, p, k)
  if(length(alpha_mat) == 1) alpha_mat = matrix(alpha_mat, k, k)
  if(length(beta_mat) == 1) beta_mat = matrix(beta_mat, k, k)
  if(length(lambda0_list) == 1) lambda0_list = rep(lambda0_list, k)
  bind_rows(lapply(1:n, function(id) gen0(id, m_vec, gamma_mat, alpha_mat, beta_mat, lambda0_list, start_func, censor_func, p_drop, tmax)))
}
# find effective sample size
find_ss <- function(k, k0, m, data, total = T) {
  res <- rep(0, m)
  id0 <- 0
  for(i in 1:nrow(data)) {
    if(data$id[i] != id0) {
      id0 <- data$id[i]
      v <- c(rep(NA, m), data$delta[i])
    }
    else {
      v <- c(v[2:(m+1)], data$delta[i])
      if(v[m+1] == k) {
        for(i in 1:m) res[m+1-i] = res[m+1-i] + (!is.na(v[i]) & v[i] == k0)
      }
    }
  }
  if(total) return(sum(res))
  else return(res)
}
*/

/*** R
pl <- function(e, par, cov.list, m, k, data) {
  p <- length(cov.list)
  gamma <- par[1:p]
  alpha <- par[(p+1):(p+k)]
  beta <- par[(p+k+1):(p+k*2)]
  data <- data %>%
    mutate(id = id %>% factor() %>% as.integer() - 1)
  X <- data %>%
    select(all_of(cov.list)) %>%
    as.matrix()
  with(data, pl_c(e, m, gamma, alpha, beta, id, time, delta, inout, X))
}
*/
// [[Rcpp::export]]
double pl_c(int e, int m, vec gamma, vec alpha, vec beta, vec id, vec time, vec delta, vec inout, mat X)  {
  int n = max(id) + 1;
  int p = gamma.n_elem;
  // Sort all observations by time.
  uvec idx = stable_sort_index(time);
  
  double pl = 0;
  double t = 0;
  mat Z = zeros(n, p);
  mat Z0 = Z;
  mat H = zeros(n, m);
  mat H0 = H;
  mat alpham = H;
  mat alpham0 = H;
  mat betam = H;
  mat betam0 = H;
  vec Y = zeros(n);
  vec Y0 = Y;
  for(uword i : idx) {
    if(time(i) != t) {
      Y = Y0;
      Z = Z0;
      H = H0;
      alpham = alpham0;
      betam = betam0;
      t = time(i);
    }
    if(delta(i) > 0) {
      if(delta(i) == e) {
        double sum_exp = 0;
        for(int j = 0; j < n; j++) {
          if(Y(j) == 1) {
            sum_exp += exp(as_scalar(Z.row(j) * gamma + alpham.row(j) * exp(-betam.row(j).t() % (time(i) - H.row(j).t()))));
          }
        }
        //std::cout << "id=" << id(i) << ", time=" << time(i) << ", Y=" << Y.t() << endl;
        //std::cout << "e=" << e << " t=" << time(i) << " i=" << i << " sum_exp=" << sum_exp << std::endl;
        pl += as_scalar(Z.row(id(i)) * gamma + alpham.row(id(i)) * exp(-betam.row(id(i)).t() % (time(i) - H.row(id(i)).t()))) - log(sum_exp);
        //std::cout << "pl=" << pl << std::endl;
        //if(sum_exp == 0) {
        //  std::cout << "H=" << H << std::endl;
        //  std::cout << "Z=" << Z << std::endl;
        //}
      }
      
      if(m > 1) {
        H0(id(i), span(1, m - 1)) = H0(id(i), span(0, m - 2));
        alpham0(id(i), span(1, m - 1)) = alpham0(id(i), span(0, m - 2));
        betam0(id(i), span(1, m - 1)) = betam0(id(i), span(0, m - 2));
      }
      H0(id(i), 0) = time(i);
      alpham0(id(i), 0) = alpha(delta(i) - 1);
      betam0(id(i), 0) = beta(delta(i) - 1);
    }
    Y0(id(i)) = inout(i);
    Z0.row(id(i)) = X.row(i);
  }
  return pl;
}

/*** R
gr <- function(e, par, cov.list, m, k, data) {
  p = length(cov.list)
  gamma = par[1:p]
  alpha = par[(p+1):(p+k)]
  beta = par[(p+k+1):(p+k*2)]
  data <- data %>%
    mutate(id = id %>% factor() %>% as.integer() - 1)
  X <- data %>%
    select(all_of(cov.list)) %>%
    as.matrix()
  with(data, gr_c(e, m, gamma, alpha, beta, id, time, delta, inout, X))
}
*/
// [[Rcpp::export]]
vec gr_c(int e, int m, vec gamma, vec alpha, vec beta,
             vec id, vec time, vec delta, vec inout, mat X) {
  int n = max(id) + 1;
  int p = gamma.n_elem;
  int k = alpha.n_elem;
  // Sort all observations by time.
  uvec idx = stable_sort_index(time);
  
  double t = 0;
  mat Z = zeros(n, p);
  mat Z0 = Z;
  mat H = zeros(n, m);
  mat H0 = H;
  mat E = H;
  mat E0 = H;
  mat alpham = H;
  mat alpham0 = H;
  mat betam = H;
  mat betam0 = H;
  vec Y = zeros(n);
  vec Y0 = Y;
  
  vec gr = zeros(p + k + k);
  for(uword i : idx) {
    if(time(i) != t) {
      Y = Y0;
      Z = Z0;
      H = H0;
      E = E0;
      alpham = alpham0;
      betam = betam0;
      t = time(i);
    }
    if(delta(i) > 0) {
      if(delta(i) == e) {
        double sum_exp = 0;
        vec sum_exp_Z = zeros(p);
        vec sum_exp_beta = zeros(k);
        vec sum_exp_alpha_beta = zeros(k);
        
        for(int j = 0; j < n; j++) {
          if(Y(j) == 1) {
            double agg = exp(as_scalar(Z.row(j) * gamma + alpham.row(j) * exp(-betam.row(j).t() % (time(i) - H.row(j).t()))));
            sum_exp += agg;
            sum_exp_Z += agg * Z.row(j).t();
            for(int k0 = 0; k0 < k; k0++) {
              sum_exp_beta(k0) += agg * as_scalar((E.row(j) == k0 + 1) * exp(-betam.row(j).t() % (time(i) - H.row(j).t())));
              sum_exp_alpha_beta(k0) -= agg * as_scalar((E.row(j) == k0 + 1) * (alpham.row(j).t() % (time(i) - H.row(j).t()) % exp(-betam.row(j).t() % (time(i) - H.row(j).t()))));
            }
          }
        }
        gr(span(0, p - 1)) += Z.row(id(i)).t() - sum_exp_Z / sum_exp;
        for(int k0 = 0; k0 < k; k0++) {
          gr(p + k0) += as_scalar((E.row(id(i)) == k0 + 1) * exp(-betam.row(id(i)).t() % (time(i) - H.row(id(i)).t())));
          gr(p + k + k0) -= as_scalar((E.row(id(i)) == k0 + 1) * (alpham.row(id(i)).t() % (time(i) - H.row(id(i)).t()) % exp(-betam.row(id(i)).t() % (time(i) - H.row(id(i)).t()))));
        }
        gr(span(p, p + k - 1)) -= sum_exp_beta / sum_exp;
        gr(span(p + k, p + k + k - 1)) -= sum_exp_alpha_beta / sum_exp;
      }
      
      if(m > 1) {
        H0(id(i), span(1, m - 1)) = H0(id(i), span(0, m - 2));
        E0(id(i), span(1, m - 1)) = E0(id(i), span(0, m - 2));
        alpham0(id(i), span(1, m - 1)) = alpham0(id(i), span(0, m - 2));
        betam0(id(i), span(1, m - 1)) = betam0(id(i), span(0, m - 2));
      }
      H0(id(i), 0) = time(i);
      E0(id(i), 0) = delta(i);
      alpham0(id(i), 0) = alpha(delta(i) - 1);
      betam0(id(i), 0) = beta(delta(i) - 1);
    }
    Y0(id(i)) = inout(i);
    Z0.row(id(i)) = X.row(i);
  }
  return gr;
}

/*** R
hs <- function(e, par, cov.list, m, k, data) {
  p = length(cov.list)
  gamma = par[1:p]
  alpha = par[(p+1):(p+k)]
  beta = par[(p+k+1):(p+k*2)]
  data <- data %>%
    mutate(id = id %>% factor() %>% as.integer() - 1)
  X <- data %>%
    select(all_of(cov.list)) %>%
    as.matrix()
  with(data, hs_c(e, m, gamma, alpha, beta, id, time, delta, inout, X))
}
*/
// [[Rcpp::export]]
mat hs_c(int e, int m, vec gamma, vec alpha, vec beta,
             vec id, vec time, vec delta, vec inout, mat X) {
  int n = max(id) + 1;
  int p = gamma.n_elem;
  int k = alpha.n_elem;
  // Sort all observations by time.
  uvec idx = stable_sort_index(time);
  
  double t = 0;
  mat Z = zeros(n, p);
  mat Z0 = Z;
  mat H = zeros(n, m);
  mat H0 = H;
  mat E = H;
  mat E0 = H;
  mat alpham = H;
  mat alpham0 = H;
  mat betam = H;
  mat betam0 = H;
  vec Y = zeros(n);
  vec Y0 = Y;
  
  mat hs = zeros(p + k + k, p + k + k);
  for(uword i : idx) {
    if(time(i) != t) {
      Y = Y0;
      Z = Z0;
      H = H0;
      E = E0;
      alpham = alpham0;
      betam = betam0;
      t = time(i);
    }
    if(delta(i) > 0) {
      // define phi as gamma*z + alpha * exp(-beta * (t - tlag))
      if(delta(i) == e) {
        double sum_exp = 0;
        vec sum_exp_Z = zeros(p);
        mat sum_exp_ZZ = zeros(p, p);
        mat sum_exp_beta2 = zeros(k, k);
        mat sum_exp_alpha_beta2 = zeros(k, k);
        vec sum_exp_beta = zeros(k);
        vec sum_exp_alpha_beta = zeros(k);
        mat sum_exp_Z_beta = zeros(p, k);
        mat sum_exp_Z_alpha_beta = zeros(p, k);
        mat sum_exp_alpha_beta_beta = zeros(k, k);
        
        for(int j = 0; j < n; j++) {
          if(Y(j) == 1) {
            double agg = exp(as_scalar(Z.row(j) * gamma + alpham.row(j) * exp(-betam.row(j).t() % (time(i) - H.row(j).t()))));
            sum_exp += agg;                                                       // exp(phi)
            sum_exp_Z += agg * Z.row(j).t();                                      // d[exp(phi)]/d[Z]
            sum_exp_ZZ += agg * Z.row(j).t() * Z.row(j);                          // d^2[exp(phi)]/d[Z]d[Z']
            
            vec exp_beta = zeros(k);
            vec exp_alpha_beta = zeros(k);
            vec exp_alpha_beta2 = zeros(k);
            for(int k0 = 0; k0 < k; k0++) {
              exp_beta(k0) = agg * as_scalar((E.row(j) == k0 + 1) * exp(-betam.row(j).t() % (time(i) - H.row(j).t())));
              exp_alpha_beta(k0) = -agg * as_scalar((E.row(j) == k0 + 1) * (alpham.row(j).t() % (time(i) - H.row(j).t()) % exp(-betam.row(j).t() % (time(i) - H.row(j).t()))));
              exp_alpha_beta2(k0) = agg * as_scalar((E.row(j) == k0 + 1) * (alpham.row(j).t() % square(time(i) - H.row(j).t()) % exp(-betam.row(j).t() % (time(i) - H.row(j).t()))));
            }
            sum_exp_beta2 += exp_beta * exp_beta.t() / agg;                       // d^2[exp(phi)]/d[alpha]d[alpha']
            sum_exp_alpha_beta2 += diagmat(exp_alpha_beta2) +
              exp_alpha_beta * exp_alpha_beta.t() / agg;                          // d^2[exp(phi)]/d[beta]d[beta']
            sum_exp_beta += exp_beta;                                             // d[exp(phi)]/d[alpha]
            sum_exp_alpha_beta += exp_alpha_beta;                                 // d[exp(phi)]/d[beta]
            
            mat exp_Z_beta = Z.row(j).t() * exp_beta.t();
            mat exp_Z_alpha_beta = Z.row(j).t() * exp_alpha_beta.t();
            sum_exp_Z_beta += exp_Z_beta;                                         // d^2[exp(phi)]/d[Z]d[alpha]
            sum_exp_Z_alpha_beta += exp_Z_alpha_beta;                             // d^2[exp(phi)]/d[Z]d[beta]
            
            sum_exp_alpha_beta_beta += diagmat(exp_alpha_beta / alpha) + exp_beta * exp_alpha_beta.t() / agg; // d^2[exp(phi)]/d[alpha]d[beta]
          }
        }
        hs(span(0, p - 1), span(0, p - 1)) += (sum_exp_Z * sum_exp_Z.t() / sum_exp - sum_exp_ZZ) / sum_exp;
        hs(span(0, p - 1), span(p, p + k - 1)) += (sum_exp_Z * sum_exp_beta.t() / sum_exp - sum_exp_Z_beta) / sum_exp;
        hs(span(0, p - 1), span(p + k, p + k + k - 1)) += (sum_exp_Z * sum_exp_alpha_beta.t() / sum_exp - sum_exp_Z_alpha_beta) / sum_exp;
        hs(span(p, p + k - 1), span(p, p + k - 1)) += (sum_exp_beta * sum_exp_beta.t() / sum_exp - sum_exp_beta2) / sum_exp;
        hs(span(p, p + k - 1), span(p + k, p + k + k - 1)) += (sum_exp_beta * sum_exp_alpha_beta.t() / sum_exp - sum_exp_alpha_beta_beta) / sum_exp;
        hs(span(p + k, p + k + k - 1), span(p + k, p + k + k - 1)) += (sum_exp_alpha_beta * sum_exp_alpha_beta.t() / sum_exp - sum_exp_alpha_beta2) / sum_exp;
        for(int k0 = 0; k0 < k; k0++) {
          hs(p + k0, p + k + k0) -= as_scalar((E.row(id(i)) == k0 + 1) * ((time(i) - H.row(id(i)).t()) % exp(-betam.row(id(i)).t() % (time(i) - H.row(id(i)).t()))));
          hs(p + k + k0, p + k + k0) += as_scalar((E.row(id(i)) == k0 + 1) * (alpham.row(id(i)).t() % (time(i) - H.row(id(i)).t()) % (time(i) - H.row(id(i)).t()) % exp(-betam.row(id(i)).t() % (time(i) - H.row(id(i)).t()))));
        }
        hs = symmatu(hs);
      }
      
      if(m > 1) {
        H0(id(i), span(1, m - 1)) = H0(id(i), span(0, m - 2));
        E0(id(i), span(1, m - 1)) = E0(id(i), span(0, m - 2));
        alpham0(id(i), span(1, m - 1)) = alpham0(id(i), span(0, m - 2));
        betam0(id(i), span(1, m - 1)) = betam0(id(i), span(0, m - 2));
      }
      H0(id(i), 0) = time(i);
      E0(id(i), 0) = delta(i);
      alpham0(id(i), 0) = alpha(delta(i) - 1);
      betam0(id(i), 0) = beta(delta(i) - 1);
    }
    Y0(id(i)) = inout(i);
    Z0.row(id(i)) = X.row(i);
  }
  return hs;
}

// Cumulative baseline hazard
/*** R
cH <- function(e, par, cov.list, m, k, data) {
  p <- length(cov.list)
  gamma <- par[1:p]
  alpha <- par[(p+1):(p+k)]
  beta <- par[(p+k+1):(p+k*2)]
  data <- data %>%
    mutate(id = id %>% factor() %>% as.integer() - 1)
  X <- data %>%
    select(all_of(cov.list)) %>%
    as.matrix()
  with(data, cH_c(e, m, gamma, alpha, beta, id, time, delta, inout, X))
}
*/
// [[Rcpp::export]]
mat cH_c(int e, int m, vec gamma, vec alpha, vec beta, vec id, vec time, vec delta, vec inout, mat X)  {
  int n = max(id) + 1;
  int p = gamma.n_elem;
  // Sort all observations by time.
  uvec idx = stable_sort_index(time);
  
  double t = 0;
  mat Z = zeros(n, p);
  mat Z0 = Z;
  mat H = zeros(n, m);
  mat H0 = H;
  mat alpham = H;
  mat alpham0 = H;
  mat betam = H;
  mat betam0 = H;
  vec Y = zeros(n);
  vec Y0 = Y;
  
  mat cH = zeros(2, sum(delta == e));
  double cH0 = 0;
  int pos = 0;
  for(uword i : idx) {
    if(time(i) != t) {
      Y = Y0;
      Z = Z0;
      H = H0;
      alpham = alpham0;
      betam = betam0;
      t = time(i);
    }
    if(delta(i) > 0) {
      if(delta(i) == e) {
        double sum_exp = 0;
        for(int j = 0; j < n; j++) {
          if(Y(j) == 1) {
            sum_exp += exp(as_scalar(Z.row(j) * gamma + alpham.row(j) * exp(-betam.row(j).t() % (time(i) - H.row(j).t()))));
          }
        }
        cH0 += 1 / sum_exp;
        cH(0, pos) = time(i);
        cH(1, pos) = cH0;
        pos += 1;
      }
      
      if(m > 1) {
        H0(id(i), span(1, m - 1)) = H0(id(i), span(0, m - 2));
        alpham0(id(i), span(1, m - 1)) = alpham0(id(i), span(0, m - 2));
        betam0(id(i), span(1, m - 1)) = betam0(id(i), span(0, m - 2));
      }
      H0(id(i), 0) = time(i);
      alpham0(id(i), 0) = alpha(delta(i) - 1);
      betam0(id(i), 0) = beta(delta(i) - 1);
    }
    Y0(id(i)) = inout(i);
    Z0.row(id(i)) = X.row(i);
  }
  return cH;
}
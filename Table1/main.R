#### Table 1
# Simulate data with p=1 (g1=0.5, B(1,0.5)), m=3 (a1=a2=a3=0.5, b1=b2=b3=0.5)
# estimate with gamma1, a1, a2, a3, b1, b2, b3
# estimate with gamma1, a1, a2, a3, a4, b1=b2=b3=b4

library(parallel)
library(dplyr)

main <- function(seed, node, N, n) {
  cl <- makeCluster(node)
  clusterExport(cl, c("N", "n"), envir = environment())
  clusterSetRNGStream(cl, seed)
  
  result <- bind_rows(clusterEvalQ(cl, {
    library(Rcpp)
    library(RcppArmadillo)
    library(dplyr)
    
    sourceCpp("~/paper1/functions_st.cpp")
    
    est1 <- function(p, m, data) {
      res <- try(nlminb(
        rep(0.1, p + m * 2),
        function(par) -pl(par, paste0("X", 1:p), m, data),
        function(par) -gr(par, paste0("X", 1:p), m, data)
      ), silent = T)
      if(inherits(res, "try-error")) {
        res$par <- rep(NA, p + m * 2)
        sd <- rep(NA, p + m * 2)
        res$convergence <- -1
        res$objective <- NA
      }
      else {
        V <- -hs(res$par, paste0("X", 1:p), m, data)
        sd <- try(sqrt(diag(solve(V))), silent = T)
        if(inherits(sd, "try-error")) {
          sd <- rep(NA, p + m * 2)
          res$convergence <- -2
        }
      }
      cum_H <- cH(res$par, paste0("X", 1:p), m, data)
      c(res$par, sd, max(cum_H[2, cum_H[1, ] <= 1]), max(cum_H[2, cum_H[1, ] <= 2]), max(cum_H[2, cum_H[1, ] <= 3]), res$convergence, -res$objective)
    }
    
    est2 <- function(p, m, data) {
      res <- try(nlminb(
        rep(0.1, p + m + 1),
        function(par) -pl(c(par[-(p+m+1)], rep(par[p+m+1], m)), paste0("X", 1:p), m, data),
        function(par) {gr0 <- gr(c(par[-(p+m+1)], rep(par[p+m+1], m)), paste0("X", 1:p), m, data); -c(gr0[1:(p+m)], sum(gr0[(p+m+1):(p+m+m)]))}
      ), silent = T)
      if(inherits(res, "try-error")) {
        res$par <- rep(NA, p + m + 1)
        sd <- rep(NA, p + m + 1)
        res$convergence <- -1
        res$objective <- NA
      }
      else {
        hs0 <- -hs(c(res$par[-(p+m+1)], rep(res$par[p+m+1], m)), paste0("X", 1:p), m, data)
        if(m == 1) V <- hs0
        else {
          hs00 <- cbind(hs0[1:(p+m+m), 1:(p+m)], rowSums(hs0[1:(p+m+m), (p+m+1):(p+m+m)]))
          V <- rbind(hs00[1:(p+m), 1:(p+m+1)], colSums(hs00[(p+m+1):(p+m+m), 1:(p+m+1)]))
        }
        sd <- try(sqrt(diag(solve(V))), silent = T)
        if(inherits(sd, "try-error")) {
          sd <- rep(NA, p + m + 1)
          res$convergence <- -2
        }
      }
      cum_H <- cH(c(res$par[-(p+m+1)], rep(res$par[p+m+1], m)), paste0("X", 1:p), m, data)
      c(res$par, sd, max(cum_H[2, cum_H[1, ] <= 1]), max(cum_H[2, cum_H[1, ] <= 2]), max(cum_H[2, cum_H[1, ] <= 3]), res$convergence, -res$objective)
    }
    
    sim <- function() {
      simdata <- gen(n, 1, 3, 0.5, 0.5, 0.5)
      c(n, sum(simdata$delta), sum(simdata$n > 1), sum(simdata$n > 2), sum(simdata$n > 3), est1(1, 3, simdata), est2(1, 4, simdata))
    }
    
    result <- t(replicate(N, sim()))
    rownames(result) <- NULL
    colnames(result) <- c("n_sub", paste0("n_event", 0:3),
                          "gamma1", paste0("alpha1", 1:3), paste0("beta1", 1:3),
                          paste0(c("gamma1", paste0("alpha1", 1:3), paste0("beta1", 1:3)), "_sd"),
                          paste0("cH1", 1:3), "conv1", "pl1",
                          "gamma2", paste0("alpha2", 1:4), "beta2",
                          paste0(c("gamma2", paste0("alpha2", 1:4), "beta2"), "_sd"),
                          paste0("cH2", 1:3), "conv2", "pl2")
    
    data.frame(result)
  }))
  
  stopCluster(cl)
  return(result)
}

#test <- main(1, 2, 2, 50)

args <- commandArgs(trailingOnly = T)
seed <- as.numeric(args[1])
node <- as.numeric(args[2])
N <- as.numeric(args[3])
n <- as.numeric(args[4])
outdir <- args[5]
char <- paste0("n",n,"_",seed)

cat(args, "\n")

assign(char, main(seed, node, N, n))
save.image(paste0(outdir, "/", char, ".RData"))
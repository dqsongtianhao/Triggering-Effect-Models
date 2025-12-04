#### Table 2
# Simulate data with p=1 (g1=0.5, B(1,0.5)), m=3 (a1=a2=a3=0.5, b1=b2=b3=0.5), k=3 (symmetric)
# estimate with gamma1, a1, a2, a3, b1, b2, b3

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
    
    #sourceCpp("C:/Users/Tianhao Song/OneDrive/Research/Recurrent Event Mar22/paper tables/functions_ct.cpp")
    sourceCpp("/nas/longleaf/home/dqsth/paper1/functions_ct.cpp")
    
    est <- function(e, p, m, k, data) {
      res <- try(nlminb(
        rep(0.1, p + k * 2),
        function(par) -pl(e, par, paste0("X", 1:p), m, k, data),
        function(par) -gr(e, par, paste0("X", 1:p), m, k, data),
        function(par) -hs(e, par, paste0("X", 1:p), m, k, data)
      ), silent = T)
      if(inherits(res, "try-error")) {
        res$par <- rep(NA, p + k * 2)
        sd <- rep(NA, p + k * 2)
        res$convergence <- -1
        res$objective <- NA
      }
      else {
        V <- -hs(e, res$par, paste0("X", 1:p), m, k, data)
        sd <- try(sqrt(diag(solve(V))), silent = T)
        if(inherits(sd, "try-error")) {
          sd <- rep(NA, p + k * 2)
          res$convergence <- -2
        }
      }
      cum_H <- cH(e, res$par, paste0("X", 1:p), m, k, data)
      c(res$par, sd, max(cum_H[2, cum_H[1, ] <= 1]), max(cum_H[2, cum_H[1, ] <= 2]), max(cum_H[2, cum_H[1, ] <= 3]), res$convergence, -res$objective)
    }
    
    sim <- function() {
      simdata <- gen(n, 1, 3, 3, 0.5, 0.5, 0.5)
      c(n, simdata %>% filter(delta == 1) %>% nrow(), find_ss(1, 1, 3, simdata), find_ss(1, 2, 3, simdata), find_ss(1, 3, 3, simdata), est(1, 1, 3, 3, simdata))
    }
    
    result <- t(replicate(N, sim()))
    rownames(result) <- NULL
    colnames(result) <- c("n_sub", paste0("n_event", 0:3),
                          "gamma", paste0("alpha", 1:3), paste0("beta", 1:3), paste0(c("gamma", paste0("alpha", 1:3), paste0("beta", 1:3)), "_sd"), paste0("cH", 1:3), "conv", "pl")
    
    data.frame(result)
  }))
  
  stopCluster(cl)
  return(result)
}

#test <- main(2, 2, 2, 20)

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
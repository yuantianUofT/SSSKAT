# SS-STAAR
R package for Semi-Supervised STAAR and Semi-Supervised gene-based association tests

# Installation and example
```r
  devtools::install_github(repo = "https://github.com/yuantianUofT/SSSKAT")

  # generate covariates
   X <- cbind(rnorm(100000), rbinom(100000, 1, 0.5))
  # generate phenotype
   y.p <- -1.7  + 0.5*X[, 1] + 0.5*X[, 2]
   Y <- rbinom(100000, 1, plogis(y.p))
  # generate genotype
   maf <- readRDS("~/SS-simulation/Code/parameter files/maf_1.rds")
   G <- lapply(1:length(maf), function(i) {
      geno <- rbinom(100000, 2, maf[i])
      return(geno)
    })
   G <- do.call(cbind, G)
  # generate surrogate
   S <- S_normal(100000, y=Y, m1 = 1, m0 = 0, s1 = 0.5, s0 = 0.5)
   id.t = sample(1:100000, 200)

  # Estimate SS-parameter under the null
   para_est <- SSSKAT::ssl_theta(Y = Y, X = X, S = S, Z = cbind(1, X), 
                                  id.t = id.t, weights = NULL, full_eval = TRUE, NULL_nlog_like, nit, "normal")

  # Obtain parametric boostrap results for test statistics variance estimation
   para_results <- data.frame()
   for (i in 1:5) {
     simseed <- i
     set.seed(simseed)
     out_i <- SSSKAT::para_func(nn = ntotal, theta = para_est$final_est, Y = Y,
                                 X = X, S = S, Z = cbind(1, X), id.t = id.t, distri = "normal")
     para_results <- rbind(para_results, out_i)
     }

  # Obtain SS-RV test pvalues
   SS_result <- SSSKAT::SS_test(Y = Y, para_results = para_results, X = X, G = G, S = S, id.t = id.t,
                                  wBurden = rep(1, 50), wSKAT = rep(1, 50), wACAT = rep(1, 50),
                                  weights.beta = NULL, mac.thresh = 10,
                                  full_NR_evaluation = TRUE, nit = NULL, NULL_nlog_like,
                                  testtype = "all", boot = T, theta = para_est$final_est, distri = normal)
  # Combine SS-RV test results to SS-STAAR result
   SS_p_STAARO <- ACAT(Pvals = c(SS_result$pvalue))

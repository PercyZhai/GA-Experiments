# ---------------------------------------------------------
# functions/bootstrap_sample.R
# ---------------------------------------------------------
# Block bootstrap draws for the covariance and precision estimation errors
# (paper Sec. VII; block length l, draws scaled by l^(-1/2)).
# Used for the bootstrap draws in all simulation figures/tables
#   toeplitz: Figs. 1-9, 14-16, Tables I-VI;  banded: Figs. 17-28, Tables VIII-XIII
# Usage:   source("functions/bootstrap_sample.R")
# Sourced by simulation/lowdim_sim.R and simulation/highdim_sim.R; not run on its own.
#
#   cov.error.boot(X, l)  : one draw approximating sqrt(n) * (Sigma.hat - Sigma)
#   prec.error.boot(X, l) : one draw approximating sqrt(n) * (Omega.hat - Omega)
#
# Inputs:  X, an n x p data matrix (rows = time points); l, the block length.
# Outputs: a p x p matrix (the simulation scripts keep its max absolute entry).
# Each call samples ONE random block of l consecutive rows, ending at a time
# index drawn uniformly from {l, ..., n} (one sample() call per function call).
# The data are demeaned globally (full-sample mean) before the block sum is formed,
# and Sigma.hat is the centered estimator cov(X) * (n - 1) / n.
# ---------------------------------------------------------

# 1. Block bootstrap for covariance
cov.error.boot <- function(X, l){
  n <- nrow(X)

  # Center with the full-sample mean, so that the block sum of X_k X_k^T
  # is comparable with Sigma.hat
  X.cent <- scale(X, center = TRUE, scale = FALSE)
  Sigma.hat <- cov(X) * (n - 1) / n

  # Block of l consecutive rows ending at a random index i in {l, ..., n}
  i <- sample(seq(l, n), 1)
  X_block <- X.cent[(i - l + 1):i, , drop = FALSE]

  # B.check = sum_{k in block} X_k X_k^T - l * Sigma.hat
  block_sum <- crossprod(X_block)
  B.check.i <- block_sum - (l * Sigma.hat)

  output <- B.check.i / sqrt(l)

  return(output)
}

# 2. Block bootstrap for precision
prec.error.boot <- function(X, l){
  n <- nrow(X)

  # Center with the full-sample mean
  X.cent <- scale(X, center = TRUE, scale = FALSE)
  Sigma.hat <- cov(X) * (n - 1) / n
  Omega.hat <- solve(Sigma.hat)

  # Block of l consecutive rows ending at a random index i in {l, ..., n}
  i <- sample(seq(l, n), 1)
  X_block <- X.cent[(i - l + 1):i, , drop = FALSE]

  # B.check^Omega = Omega.hat (sum_{k in block} X_k X_k^T - l * Sigma.hat) Omega.hat
  inner_term <- crossprod(X_block) - (l * Sigma.hat)
  B.check.i <- Omega.hat %*% inner_term %*% Omega.hat

  output <- B.check.i / sqrt(l)

  return(output)
}

# Entrywise block-bootstrap confidence intervals for a precision matrix.
#
# Used for the ABIDE conditional independence graphs of Fig. 12 (paper Sec. VIII)
# by real-data/03_cig_bootstrap.R.
#
#   prec.CI.boot(X, l, epsilon = 0.5, alpha = 0.05)
#     X      n x p data matrix (rows are time points); the second-moment estimate
#            t(X) %*% X / n is used as Sigma.hat, i.e. X is taken as mean zero
#            (the ABIDE ROI series are already demeaned)
#     l      block length
#     alpha  level: each interval Omega.hat[s,t] -/+ width[s,t] has nominal
#            coverage 1 - alpha
#     epsilon  not used
#   Returns list(Omega.top, Omega.btm, Omega.hat, width, lB.check.st): upper and
#   lower confidence bounds, Omega.hat = solve(Sigma.hat), the half-widths and the
#   bootstrap replicates of the last entry (s = t = p).
#
# Usage:   source("functions/bootstrap_ci.R")   (from the repository root)
# Inputs:  none
# Outputs: none (defines prec.CI.boot())

# Block bootstrap for precision
prec.CI.boot <- function(X, l, epsilon=0.5, alpha = 0.05){
  n <- nrow(X)
  p <- ncol(X)
  Sigma.hat <- t(X) %*% X / n
  Omega.hat <- solve(Sigma.hat)

  # 1. linearized summands Omega.hat (X_j X_j' - Sigma.hat) Omega.hat, j = 1, ..., n
  X.check <- list()
  for(j in 1:n){
    X.check[[j]] <- Omega.hat %*% (X[j,] %*% t(X[j,]) - Sigma.hat) %*% Omega.hat
  }

  # 2. block sums B.check_{i,l} = sum_{k=i-l+1}^{i} X.check_k for i in l:n
  B.check <- list() #for B.check.i, it is the i-l+1 th entry in the list
  for(i in l:n){
    B.check[[i-l+1]] <- matrix(0, p, p)
    for(k in (i-l+1):i){
      B.check[[i-l+1]] <- B.check[[i-l+1]] + X.check[[k]]
    }
  }

  # 3. for each entry, the half-width is the floor((1-alpha)*(n-l+1))-th smallest of
  #    |B.check_{i,l}[s,t]| / sqrt(l) over the n-l+1 blocks, divided by sqrt(n)
  width <- matrix(NA, p, p)
  for(s in 1:p){
    for(t in s:p){
      lB.check.st <- rep(NA, n-l+1)
      for(i in 1:(n-l+1)){
        lB.check.st[i] <- B.check[[i]][s,t] / sqrt(l)
      }
      lB.check.st.abs <- sort(abs(lB.check.st))
      width[s,t] <- lB.check.st.abs[(1-alpha) * (n-l+1)] / sqrt(n)
      width[t,s] <- width[s,t]
    }
  }

  Omega.top <- Omega.hat + width
  Omega.btm <- Omega.hat - width

  output <- list(Omega.top, Omega.btm, Omega.hat, width, lB.check.st)

  return(output)
}

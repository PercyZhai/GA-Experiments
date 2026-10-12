# ---------------------------------------------------------
# simulation/lowdim_sim.R
# ---------------------------------------------------------
# Low-dimensional simulation (paper Sec. VII-A; banded setup in Appendix I).
# For one (n, p, beta, setup) cell, generates 200 copies of the Gaussian linear
# process X_i = sum_{t=0}^{N-1} A_t eps_{i-t}, A_t = (t+1)^(-beta) A, N = n^2,
# with the FFT method of Appendix F, and stores 200 values each of
#   sqrt(n) |Sigma.hat - Sigma|_inf, sqrt(n) |Omega.hat - Omega|_inf  (true errors),
#   |Z|_inf, |Z^S|_inf                                       (Gaussian approximations),
#   one block bootstrap draw for covariance and for precision, l = floor(n^(2/3)).
#
# Paper figures/tables (drawn by analysis/lowdim/ from these results):
#   toeplitz: Figs. 1-3 (QQ cov), Figs. 4-6 (QQ prec), Figs. 14, 15 (ECDF cov, prec),
#             Tables I, II (KS cov, prec), Tables IV, V (W1 cov, prec)
#   banded:   Figs. 17-19 (QQ cov), Fig. 20 (ECDF cov), Figs. 21-23 (QQ prec), Fig. 24 (ECDF prec),
#             Tables VIII, X (KS cov, prec), Tables IX, XI (W1 cov, prec)
#
# Usage (from the repository root):
#   Rscript simulation/lowdim_sim.R <n> <p> <beta> [toeplitz|banded]
# Paper grid: n in {200, 500, 1000, 2000, 4000, 10000}, p in {2, 10, 30, 100},
# beta in {2, 0.9, 0.55}; see simulation/run_grid.sh and simulation/submit.sh.
# n >= 200 is required: the 200 copies are consecutive length-n segments of one
# simulated series of length N = n^2.
#
# Inputs:  functions/bootstrap_sample.R
# Outputs: results/lowdim/E2quick_<n>_<p>_<beta>_toeplitz.RData or
#          results/lowdim/E2Bquick_<n>_<p>_<beta>_banded.RData, containing
#          Sigma.error.maxs, Omega.error.maxs, Z.max.norms, Z.Omega.max.norms,
#          Sigma.boot.maxs, Omega.boot.maxs (numeric vectors of length 200)
#
# Setup:   toeplitz (default): A_jk = (|j-k|+1)^(-2)
#          banded:             A_jk = (|j-k|+1)^(-2) * 1{|j-k| <= 2}
# With set.seed(42), the runs reproduce the published result files up to
# floating-point differences (a different R/FFT/BLAS may change the last digits).
# The largest cells (n = 10000) work with FFTs of length n^2 = 1e8.
# ---------------------------------------------------------

if (!file.exists("functions/bootstrap_sample.R")) stop("Run this script from the root of the GA-Experiments repository.")

# 1. Parse Command Line Arguments
args <- commandArgs(trailingOnly = TRUE)

# Strict Argument Checking
if(length(args) < 3){
  stop("ERROR: Missing arguments. Usage: Rscript simulation/lowdim_sim.R <n> <p> <beta> [toeplitz|banded]")
}

n <- as.numeric(args[1])
p <- as.numeric(args[2])
beta <- as.numeric(args[3])
setup <- if (length(args) >= 4) args[4] else "toeplitz"
stopifnot(setup %in% c("toeplitz", "banded"))

# Derived parameters
m <- n^2               # truncation length N of the linear process
epochs <- 200          # number of copies of {X_i}_{i=1}^n
l <- floor(n^(2/3))    # bootstrap block length

print(paste("Running with: n=", n, ", p=", p, ", beta=", beta, ", setup=", setup))

save.path <- "results/lowdim"
dir.create(save.path, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
source("functions/bootstrap_sample.R")


######################################
##        Start of the task         ##
######################################
time.start <- proc.time()

######################################
# Step 1. Generate A Matrix          #
######################################
# p x p Toeplitz matrix, A_jk = (|j-k|+1)^(-2)
toeplitz <- function(p){
  result <- matrix(NA, p, p)
  for(i in 1:p){
    for(j in 1:i){
      result[i,j] <- 1 / (abs(i-j) + 1)^2
      result[j,i] <- 1 / (abs(i-j) + 1)^2
    }
  }
  return(result)
}

A <- toeplitz(p)

# Banded setup (Appendix I): keep only the entries with |j-k| <= 2
if (setup == "banded") A <- A * (abs(row(A) - col(A)) <= 2)

# Temporal weights: factor[t] = (t+1)^(-beta) is the weight of lag t, t = 1, ..., m-1;
# factor[m] = 1 is the lag-0 weight (the FFT convolution below is circular)
factor <- rep(0, m)
for(t in 1:(m-1)){
  factor[t] <- (t+1)^(-beta)
}
factor[m] <- 1

######################################
# Step 2. True Covariance            #
######################################
# True covariance and precision matrices; A.Omega = Omega A gives the
# intermediate Gaussian approximation Z^S of the precision error (Step 5)

Sigma <- sum(factor^2) * A %*% t(A)
Omega <- solve(Sigma)
A.Omega <- Omega %*% A

######################################
# Step 3 & 4. FFT Generation         #
# (Appendix F; Sequential and        #
#  Memory Safe)                      #
######################################

cat("Generating data using Optimized (Low-Memory) Method...\n")

# 1. Pre-calculate FFT of the time decay factor
if(length(factor) < m) {
  factor_padded <- c(factor, rep(0, m - length(factor)))
} else {
  factor_padded <- factor
}

# Transform filter to frequency domain
b_vec <- fft(factor_padded) / sqrt(m)

# 2. Prepare Storage
# We only store the rows we actually need
rows_needed <- n * epochs
Z_trunc <- matrix(NA, nrow = rows_needed, ncol = p)

# 3. Generate p INDEPENDENT temporally correlated series
cat(paste("Generating", p, "series sequentially...\n"))

for(j in 1:p){
  # A. Generate white noise
  e_j <- rnorm(m)

  # B. FFT of noise
  f_j <- fft(e_j) / sqrt(m)

  # C. Apply filter in frequency domain
  w_j_complex <- fft(b_vec * f_j, inverse = TRUE)

  # D. Store ONLY the real part and ONLY the truncated rows
  Z_trunc[, j] <- Re(w_j_complex)[1:rows_needed]

  # E. Aggressive Garbage Collection
  rm(e_j, f_j, w_j_complex)
  gc()

  if(j %% 10 == 0) cat(paste("Finished series", j, "of", p, "\n"))
}

# 4. Apply Spatial Correlation: row i of W.multi is X_i = A Z_i
cat("Applying spatial correlation matrix A...\n")
W.multi <- Z_trunc %*% t(A)

# Clean up
rm(Z_trunc, b_vec, factor_padded, factor)
gc()

time.end <- proc.time()
runtime <- time.end[3] - time.start[3]
print(paste("Data Generation Runtime:", round(runtime,3)))


#####################################
# Step 5. Oracle GA (FFT Method)    #
#####################################
# Sigma.Z.factor = sum_{|k| < n} (1 - |k|/n) gamma_k^2, where gamma_k is the lag-k
# autocovariance of the scalar weights a_t = (t+1)^(-beta) (computed by FFT).
# Together with A it gives the covariance of Z (Appendix B-F, eq. (A34)).

Sigma.Z.factor.comp <- function(n, beta, m) {
  t_max <- m + n
  a <- (1:t_max)^(-beta)

  # Use FFT to compute autocorrelation
  N_fft <- 2^ceiling(log2(t_max + n))
  a_padded <- c(a, rep(0, N_fft - length(a)))

  F_a <- fft(a_padded)
  Power <- Mod(F_a)^2
  Autocorr <- Re(fft(Power, inverse = TRUE)) / N_fft

  gammas <- Autocorr[1:n]

  gamma_0_sq <- gammas[1]^2

  if (n > 1) {
    k_vec <- 1:(n - 1)
    weights <- (n - k_vec) / n
    gamma_rest_sq <- gammas[2:n]^2
    total_sum <- gamma_0_sq + 2 * sum(weights * gamma_rest_sq)
  } else {
    total_sum <- gamma_0_sq
  }

  return(total_sum)
}

cat("Computing Sigma Z factor (FFT method)\n")
Sigma.Z.factor <- Sigma.Z.factor.comp(n, beta, m)

# Generate Z = A G A^T and Z^S = (Omega A) G (Omega A)^T, where G is a symmetric
# Gaussian matrix with Var(G_jk) = Sigma.Z.factor * (1 + 1{j = k})
Z.max.norms <- rep(NA, epochs)
Z.Omega.max.norms <- rep(NA, epochs)

cat("Generating Gaussian Approximation Z\n")

A_small <- A
A_small_Omega <- A.Omega

for(epoch in 1:epochs){
  X <- matrix(rnorm(p * p), nrow = p, ncol = p)

  Y <- A_small %*% X %*% t(A_small)
  Z_mat <- sqrt(Sigma.Z.factor / 2) * (Y + t(Y))

  Y.Omega <- A_small_Omega %*% X %*% t(A_small_Omega)
  Z_mat.Omega <- sqrt(Sigma.Z.factor / 2) * (Y.Omega + t(Y.Omega))

  Z.max.norms[epoch] <- max(abs(Z_mat))
  Z.Omega.max.norms[epoch] <- max(abs(Z_mat.Omega))
}


#####################################
# Step 6. Bootstrap and true error  #
#####################################

cat("Computing Sigma and Omega errors\n")

Sigma.error.maxs <- rep(NA, epochs)
Omega.error.maxs <- rep(NA, epochs)
Sigma.boot.maxs <- rep(NA, epochs)
Omega.boot.maxs <- rep(NA, epochs)

for(epoch in 1:epochs){

  # Rows of copy number 'epoch'
  # Since W.multi is already truncated to (n*epochs), we index it directly
  idx_start <- (epoch-1)*n + 1
  idx_end <- epoch*n
  X <- W.multi[idx_start:idx_end, ]

  Sigma.hat <- cov(X) * (n - 1) / n
  Sigma.error.maxs[epoch] <- max(abs(Sigma.hat - Sigma)) * sqrt(n)
  Sigma.boot.maxs[epoch] <- max(abs(cov.error.boot(X, l)))

  Omega.hat <- solve(Sigma.hat)
  Omega.error.maxs[epoch] <- max(abs(Omega.hat - Omega)) * sqrt(n)
  Omega.boot.maxs[epoch] <- max(abs(prec.error.boot(X, l)))

}

time.end <- proc.time()
runtime <- time.end[3] - time.start[3]
print(paste("Total Runtime:", round(runtime,3)))


#####################################
# Step 7. Save output               #
#####################################

out.prefix <- if (setup == "toeplitz") "E2quick_" else "E2Bquick_"
out.file.name <- paste0(out.prefix, n, "_", p, "_", beta, "_", setup, ".RData")
full.save.path <- file.path(save.path, out.file.name)

save(Sigma.error.maxs, Omega.error.maxs,
     Z.max.norms, Z.Omega.max.norms,
     Sigma.boot.maxs, Omega.boot.maxs,
     file = full.save.path)

print(paste("Saved results to:", full.save.path))

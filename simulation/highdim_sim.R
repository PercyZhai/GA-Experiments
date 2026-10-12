# ---------------------------------------------------------
# simulation/highdim_sim.R
# ---------------------------------------------------------
# High-dimensional simulation, covariance only (paper Sec. VII-B; banded setup in Appendix I).
# For one (n, p, beta, setup) cell, generates 200 copies of the Gaussian linear
# process X_i = sum_{t=0}^{N-1} A_t eps_{i-t}, A_t = (t+1)^(-beta) A, N = n^2,
# with the FFT method of Appendix F, and stores 200 values each of
#   sqrt(n) |Sigma.hat - Sigma|_inf   (true error),
#   |Z|_inf                           (Gaussian approximation),
#   one block bootstrap draw for covariance, l = floor(n^(2/3)).
#
# Paper figures/tables (drawn by analysis/highdim/ from these results):
#   toeplitz: Figs. 7-9 (QQ cov), Fig. 16 (ECDF cov), Table III (KS cov), Table VI (W1 cov)
#   banded:   Figs. 25-27 (QQ cov), Fig. 28 (ECDF cov), Table XII (KS cov), Table XIII (W1 cov)
#
# Usage (from the repository root):
#   Rscript simulation/highdim_sim.R <n> <p> <beta> [toeplitz|banded]
# Paper grid: n in {200, 400, 800}, p in {2000, 5000, 8000}, beta in {2, 0.9, 0.55};
# see simulation/run_grid.sh and simulation/submit.sh.
# n >= 200 is required: the 200 copies are consecutive length-n segments of one
# simulated series of length N = n^2.
#
# Inputs:  functions/bootstrap_sample.R
# Outputs: results/highdim/E3quick_<n>_<p>_<beta>_toeplitz.RData or
#          results/highdim/E3Bquick_<n>_<p>_<beta>_banded.RData, containing
#          Sigma.error.maxs, Z.max.norms, Sigma.boot.maxs (numeric vectors of length 200)
#
# Setup:   toeplitz (default): A_jk = (|j-k|+1)^(-2)
#          banded:             A_jk = (|j-k|+1)^(-2) * 1{|j-k| <= 2}
#
# The Gaussian approximation draws and the error / block bootstrap computations run in
# parallel (foreach + doParallel) on n_cores workers: SLURM_CPUS_PER_TASK if set, otherwise
# all cores but one. The data are generated with set.seed(42), so the true errors reproduce
# the published ones (up to floating-point differences); the workers' random numbers are
# not seeded, so the Gaussian approximation and bootstrap draws differ between runs.
#
# Memory and runtime: the data matrix holds n * 200 x p doubles (about 10 GB for n = 800,
# p = 8000; it is held twice while A is applied), each worker holds p x p matrices
# (0.5 GB each for p = 8000), and each Gaussian draw needs two dense p x p products. The
# published runs used 28 CPUs per cell; an optimized BLAS is recommended.
# ---------------------------------------------------------

if (!file.exists("functions/bootstrap_sample.R")) stop("Run this script from the root of the GA-Experiments repository.")

library(foreach)
library(doParallel)
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = max(1, parallel::detectCores() - 1)))

# 1. Parse Command Line Arguments
args <- commandArgs(trailingOnly = TRUE)

# Strict Argument Checking
if(length(args) < 3){
  stop("ERROR: Missing arguments. Usage: Rscript simulation/highdim_sim.R <n> <p> <beta> [toeplitz|banded]")
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

print(paste("Running E3 (Cov Only) with: n=", n, ", p=", p, ", beta=", beta, ", setup=", setup))

save.path <- "results/highdim"
dir.create(save.path, recursive = TRUE, showWarnings = FALSE)

set.seed(42)
source("functions/bootstrap_sample.R")


######################################
##         Start of the task        ##
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

Sigma <- sum(factor^2) * A %*% t(A)

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

# Generate Z = A G A^T, where G is a symmetric Gaussian matrix with
# Var(G_jk) = Sigma.Z.factor * (1 + 1{j = k})

A_small <- A

cat(paste("Generating Gaussian Approximation Z using", n_cores, "cores...\n"))

# Start Cluster
cl <- makeCluster(n_cores, outfile="")
registerDoParallel(cl)

Z.max.norms <- foreach(epoch = 1:epochs,
                       .combine = c,
                       .export = c("A_small", "p", "Sigma.Z.factor")) %dopar% {

                         # Print Progress
                         cat("Epoch #", epoch, "...\n", sep="")

                         # Generate Random Matrix
                         X <- matrix(rnorm(p * p), nrow = p, ncol = p)

                         # Dense Matrix Multiplication
                         Y <- A_small %*% X %*% t(A_small)

                         # Symmetrize and Scale
                         Z_mat <- sqrt(Sigma.Z.factor / 2) * (Y + t(Y))

                         # Return Max Norm
                         max(abs(Z_mat))
                       }

stopCluster(cl)


#####################################
# Step 6. Bootstrap and true error  #
# (Parallelized + Memory Optimized) #
#####################################

cat(paste("Computing Sigma errors using", n_cores, "cores...\n"))

# 1. Split W.multi into the 200 copies, so that each worker only receives
#    the copy it needs instead of the whole data matrix
cat("Splitting W.multi to reduce parallel memory overhead...\n")
X_list <- vector("list", epochs)
for(i in 1:epochs){
  idx_start <- (i-1)*n + 1
  idx_end <- i*n
  X_list[[i]] <- W.multi[idx_start:idx_end, ]
}

# 2. Free the full data matrix
rm(W.multi)
gc()

# 3. Start Cluster
# outfile="" allows worker 'cat' output to show in the main log
cl <- makeCluster(n_cores, outfile = "")
registerDoParallel(cl)

# 4. Run Loop over the copies
results_matrix <- foreach(X = X_list, epoch = 1:epochs,
                          .combine = rbind,
                          # Export Sigma (needed for diff) and bootstrap function
                          .export = c("n", "Sigma", "l", "cov.error.boot")) %dopar% {

                            # Print Progress
                            cat("Error Calc Epoch #", epoch, "\n", sep="")

                            # A. Empirical Error
                            Sigma.hat <- cov(X) * (n - 1) / n
                            err_val <- max(abs(Sigma.hat - Sigma)) * sqrt(n)

                            # B. Bootstrap Error
                            boot_val <- max(abs(cov.error.boot(X, l)))

                            # Return vector
                            c(err_val, boot_val)
                          }

stopCluster(cl)

# 5. Unpack Results
Sigma.error.maxs <- results_matrix[, 1]
Sigma.boot.maxs  <- results_matrix[, 2]

time.end <- proc.time()
runtime <- time.end[3] - time.start[3]
print(paste("Total Runtime:", round(runtime,3)))


#####################################
# Step 7. Save output               #
#####################################

out.prefix <- if (setup == "toeplitz") "E3quick_" else "E3Bquick_"
out.file.name <- paste0(out.prefix, n, "_", p, "_", beta, "_", setup, ".RData")
full.save.path <- file.path(save.path, out.file.name)

save(Sigma.error.maxs,
     Z.max.norms,
     Sigma.boot.maxs,
     file = full.save.path)

print(paste("Saved results to:", full.save.path))

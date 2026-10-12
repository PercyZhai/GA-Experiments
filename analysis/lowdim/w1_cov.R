# ---------------------------------------------------------
# Wasserstein-1 distances, low-dimensional simulation, covariance matrix.
# W1 distances between the empirical error max-norms (Sigma.error.maxs),
# the Gaussian approximation (Z.max.norms) and the bootstrap (Sigma.boot.maxs)
# for every (n, p, beta).
#
# Paper:   toeplitz: Table IV
#          banded:   Table IX, Appendix I
# Usage:   Rscript analysis/lowdim/w1_cov.R [toeplitz|banded]
# Inputs:  results/lowdim/<prefix>quick_<n>_<p>_<beta>_<setup>.RData
# Outputs: output/tables/<prefix>.W1.cov.matrix.csv,
#          prefix = "E2" (toeplitz) or "E2B" (banded)
# ---------------------------------------------------------
if (!file.exists("functions/bootstrap_sample.R")) stop("Run this script from the root of the GA-Experiments repository.")

library(transport)

args <- commandArgs(trailingOnly = TRUE)
setup <- if (length(args) >= 1) args[1] else "toeplitz"
stopifnot(setup %in% c("toeplitz", "banded"))
prefix <- if (setup == "toeplitz") "E2" else "E2B"

# ---------------------------------------------------------
# 1. Setup Parameters
# ---------------------------------------------------------
n.vec <- c(10000, 4000, 2000, 1000, 500, 200)
p.vec <- c(2, 10, 30, 100)
beta.vec <- c(2, 0.9, 0.55)

# Metrics definitions
metric_labels <- c("Err vs Z", "Err vs Boot", "Boot vs Z")

base.path <- paste0("results/lowdim/", prefix, "quick_")  # E2quick_ (toeplitz) or E2Bquick_ (banded)
out.file <- paste0("output/tables/", prefix, ".W1.cov.matrix.csv")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# 2. Initialize the Output Matrix
# ---------------------------------------------------------
# Rows: 3 metrics for each n; columns: each (p, beta) pair
results_mat <- matrix(NA, nrow = length(n.vec) * 3, ncol = length(p.vec) * length(beta.vec))

# Generate Row Names (e.g., "n=10000 (Err vs Z)")
row_names <- c()
for(n in n.vec){
  for(met in metric_labels){
    row_names <- c(row_names, paste0("n=", n, " (", met, ")"))
  }
}
rownames(results_mat) <- row_names

# Generate Column Names (e.g., "p=2 b=2")
col_names <- c()
for(p in p.vec){
  for(b in beta.vec){
    col_names <- c(col_names, paste0("p=", p, " b=", b))
  }
}
colnames(results_mat) <- col_names

# ---------------------------------------------------------
# 3. Helper Function to Load Data
# ---------------------------------------------------------
get_data <- function(n, p, beta) {
  fname <- paste0(base.path, n, "_", p, "_", beta, "_", setup, ".RData")
  if (file.exists(fname)) {
    e <- new.env()
    load(fname, envir = e)
    return(e)
  } else {
    warning("Missing input file: ", fname, immediate. = TRUE)
    return(NULL)
  }
}

# ---------------------------------------------------------
# 4. Main Computation Loop
# ---------------------------------------------------------
cat("Calculating Wasserstein distances...\n")

# Row Index Counter
row_idx_start <- 1

for(i in 1:length(n.vec)) {
  n_val <- n.vec[i]

  # Each n fills a block of 3 rows (one per metric); the column index
  # runs over the (p, beta) pairs and resets for every n
  col_idx <- 1

  for(j in 1:length(p.vec)) {
    p_val <- p.vec[j]

    for(k in 1:length(beta.vec)) {
      beta_val <- beta.vec[k]

      # Load data
      d <- get_data(n_val, p_val, beta_val)

      # Cells stay NA if the input file is missing
      if(!is.null(d)) {
        # 1. Err vs Z (Row 1 of this n-block)
        w1_z <- wasserstein1d(d$Sigma.error.maxs, d$Z.max.norms, p=1)
        results_mat[row_idx_start, col_idx] <- w1_z

        # 2. Err vs Boot (Row 2 of this n-block)
        w1_boot <- wasserstein1d(d$Sigma.error.maxs, d$Sigma.boot.maxs, p=1)
        results_mat[row_idx_start + 1, col_idx] <- w1_boot

        # 3. Boot vs Z (Row 3 of this n-block)
        w1_bz <- wasserstein1d(d$Sigma.boot.maxs, d$Z.max.norms, p=1)
        results_mat[row_idx_start + 2, col_idx] <- w1_bz
      }

      col_idx <- col_idx + 1
    }
  }

  # Move row pointer down by 3 for the next n
  row_idx_start <- row_idx_start + 3
}

# ---------------------------------------------------------
# 5. Output Results
# ---------------------------------------------------------

# Print to console formatted (rounded)
print(round(results_mat, 2))

# Write to CSV (unrounded)
write.csv(results_mat, file = out.file)
cat("\nSaved table to", out.file, "\n")

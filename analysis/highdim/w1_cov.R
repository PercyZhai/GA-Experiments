# ---------------------------------------------------------
# analysis/highdim/w1_cov.R
# ---------------------------------------------------------
# Description:
#   Computes Wasserstein-1 distances for Covariance Matrix
#   distributions (E3 high-dimensional experiment):
#   Err vs Z (Gaussian approximation), Err vs Boot, Boot vs Z.
#   Uses transport::wasserstein1d. Rows: n * 3 metrics; Cols: p * beta.
#
# Paper:
#   toeplitz: Table VI
#   banded:   Table XIII (Appendix I)
#
# Usage (from the repository root):
#   Rscript analysis/highdim/w1_cov.R [toeplitz|banded]   (default: toeplitz)
#
# Inputs:
#   results/highdim/<prefix>quick_<n>_<p>_<beta>_<setup>.RData
#   (written by simulation/highdim_sim.R; uses Sigma.error.maxs,
#   Z.max.norms, Sigma.boot.maxs)
#
# Outputs:
#   output/tables/<prefix>.W1.cov.matrix.csv  (unrounded),
#   prefix = "E3" (toeplitz) or "E3B" (banded)
# ---------------------------------------------------------
if (!file.exists("functions/bootstrap_sample.R")) stop("Run this script from the root of the GA-Experiments repository.")

args <- commandArgs(trailingOnly = TRUE)
setup <- if (length(args) >= 1) args[1] else "toeplitz"
stopifnot(setup %in% c("toeplitz", "banded"))
prefix <- if (setup == "toeplitz") "E3" else "E3B"

# Load required library
library(transport)

# ---------------------------------------------------------
# 1. Setup Parameters (E3)
# ---------------------------------------------------------
n.vec <- c(800, 400, 200)
p.vec <- c(2000, 5000, 8000)
beta.vec <- c(2, 0.9, 0.55)

# Metrics definitions
metric_labels <- c("Err vs Z", "Err vs Boot", "Boot vs Z")

# Input files (relative to the repository root)
base.path <- paste0("results/highdim/", prefix, "quick_")  # E3quick_ (toeplitz) or E3Bquick_ (banded)

# Output file
table.dir <- "output/tables"
dir.create(table.dir, recursive = TRUE, showWarnings = FALSE)
out.file <- paste0(table.dir, "/", prefix, ".W1.cov.matrix.csv")

# ---------------------------------------------------------
# 2. Initialize the Output Matrix
# ---------------------------------------------------------
# Rows: (Length of n) * 3 metrics
# Cols: (Length of p) * (Length of beta)
results_mat <- matrix(NA, nrow = length(n.vec) * 3, ncol = length(p.vec) * length(beta.vec))

# Generate Row Names
row_names <- c()
for(n in n.vec){
  for(met in metric_labels){
    row_names <- c(row_names, paste0("n=", n, " (", met, ")"))
  }
}
rownames(results_mat) <- row_names

# Generate Column Names
col_names <- c()
for(p in p.vec){
  for(b in beta.vec){
    col_names <- c(col_names, paste0("p=", p, " b=", b))
  }
}
colnames(results_mat) <- col_names

# ---------------------------------------------------------
# 3. Helper Functions
# ---------------------------------------------------------
result_file <- function(n, p, beta) {
  paste0(base.path, n, "_", p, "_", beta, "_", setup, ".RData")
}

get_data <- function(n, p, beta) {
  fname <- result_file(n, p, beta)
  if (file.exists(fname)) {
    e <- new.env()
    load(fname, envir = e)
    return(e)
  } else {
    warning("Result file not found: ", fname, call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
}

# Warn when a distance cannot be computed; the table cell is then left as NA
warn_failed <- function(met, n, p, beta, err) {
  warning("wasserstein1d (", met, ") failed for ", result_file(n, p, beta), ": ",
          conditionMessage(attr(err, "condition")), call. = FALSE, immediate. = TRUE)
}

# ---------------------------------------------------------
# 4. Main Computation Loop
# ---------------------------------------------------------
cat("Calculating Wasserstein-1 distances (Covariance) for", prefix, "...\n")

row_idx_start <- 1

for(i in 1:length(n.vec)) {
  n_val <- n.vec[i]
  col_idx <- 1

  for(j in 1:length(p.vec)) {
    p_val <- p.vec[j]

    for(k in 1:length(beta.vec)) {
      beta_val <- beta.vec[k]

      d <- get_data(n_val, p_val, beta_val)

      if(!is.null(d)) {
        # --- Covariance W1 Distances ---

        # 1. Err vs Z (Row 1)
        w1_z <- try(wasserstein1d(d$Sigma.error.maxs, d$Z.max.norms, p=1), silent=TRUE)
        if(!inherits(w1_z, "try-error")) {
          results_mat[row_idx_start, col_idx] <- w1_z
        } else {
          warn_failed(metric_labels[1], n_val, p_val, beta_val, w1_z)
        }

        # 2. Err vs Boot (Row 2)
        w1_boot <- try(wasserstein1d(d$Sigma.error.maxs, d$Sigma.boot.maxs, p=1), silent=TRUE)
        if(!inherits(w1_boot, "try-error")) {
          results_mat[row_idx_start + 1, col_idx] <- w1_boot
        } else {
          warn_failed(metric_labels[2], n_val, p_val, beta_val, w1_boot)
        }

        # 3. Boot vs Z (Row 3)
        w1_bz <- try(wasserstein1d(d$Sigma.boot.maxs, d$Z.max.norms, p=1), silent=TRUE)
        if(!inherits(w1_bz, "try-error")) {
          results_mat[row_idx_start + 2, col_idx] <- w1_bz
        } else {
          warn_failed(metric_labels[3], n_val, p_val, beta_val, w1_bz)
        }
      }

      col_idx <- col_idx + 1
    }
  }
  row_idx_start <- row_idx_start + 3
}

# ---------------------------------------------------------
# 5. Output Results
# ---------------------------------------------------------

# Print to console formatted (rounded)
print(round(results_mat, 4))

# Write to CSV (unrounded)
write.csv(results_mat, file=out.file)
cat("\nSaved table to", out.file, "\n")

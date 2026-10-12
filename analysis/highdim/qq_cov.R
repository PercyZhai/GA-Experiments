# ---------------------------------------------------------
# analysis/highdim/qq_cov.R
# ---------------------------------------------------------
# Description:
#   QQ plots for the high-dimensional covariance experiment (E3).
#   For each (n, p) the left panel compares the max-norm error of the
#   sample covariance with its Gaussian approximation (GA) and the right
#   panel compares it with the block bootstrap. One figure per beta.
#   Layout: header row (p) + 3 rows (n) x 6 cols (p * {GA, Boot}).
#
# Paper:
#   toeplitz: Figs. 7-9   (beta = 2, 0.9, 0.55)
#   banded:   Figs. 25-27 (beta = 2, 0.9, 0.55), Appendix I
#
# Usage (from the repository root):
#   Rscript analysis/highdim/qq_cov.R [toeplitz|banded]   (default: toeplitz)
#
# Inputs:
#   results/highdim/<prefix>quick_<n>_<p>_<beta>_<setup>.RData
#   (written by simulation/highdim_sim.R; uses Sigma.error.maxs,
#   Z.max.norms, Sigma.boot.maxs)
#
# Outputs:
#   output/figures/<prefix>.QQ.cov.<beta>.pdf  (10 x 5 in),
#   prefix = "E3" (toeplitz) or "E3B" (banded)
# ---------------------------------------------------------
if (!file.exists("functions/bootstrap_sample.R")) stop("Run this script from the root of the GA-Experiments repository.")

args <- commandArgs(trailingOnly = TRUE)
setup <- if (length(args) >= 1) args[1] else "toeplitz"
stopifnot(setup %in% c("toeplitz", "banded"))
prefix <- if (setup == "toeplitz") "E3" else "E3B"

# ---------------------------------------------------------
# Parameters
# ---------------------------------------------------------
n.vec <- c(800, 400, 200)
p.vec <- c(2000, 5000, 8000)
beta.vec <- c(2, 0.9, 0.55)

# Input files (relative to the repository root)
base.path <- paste0("results/highdim/", prefix, "quick_")  # E3quick_ (toeplitz) or E3Bquick_ (banded)

# Output directory
fig.dir <- "output/figures"
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# Layout Setup
# ---------------------------------------------------------
# Layout Matrix Definition:
# 0  1  1  2  2  3  3   (Top Header Row: labels for p)
# 4  7  8  9 10 11 12   (Row 1: n=800)
# 5 13 14 15 16 17 18   (Row 2: n=400)
# 6 19 20 21 22 23 24   (Row 3: n=200)

layout_mat <- matrix(c(
  0, 1, 1, 2, 2, 3, 3,
  4, 7, 8, 9, 10,11,12,
  5, 13,14,15,16,17,18,
  6, 19,20,21,22,23,24
), nrow=4, byrow=TRUE)

# ---------------------------------------------------------
# Helper Function to Load Data
# ---------------------------------------------------------
get_data <- function(n, p, beta) {
  fname <- paste0(base.path, n, "_", p, "_", beta, "_", setup, ".RData")
  if (file.exists(fname)) {
    e <- new.env()
    load(fname, envir = e)
    return(e)
  } else {
    warning("Result file not found: ", fname, call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
}

# ---------------------------------------------------------
# One figure per beta
# ---------------------------------------------------------
for(beta in beta.vec) {

  out.file <- paste0(fig.dir, "/", prefix, ".QQ.cov.", beta, ".pdf")
  pdf(out.file, width = 10, height = 5)

  layout(layout_mat, widths=c(0.3, 1,1,1,1,1,1), heights=c(0.2, 1,1,1))

  # Minimize internal margins for the plots
  par(mar=c(2, 2, 1, 1), oma=c(0,0,0,0))

  # ---------------------------------------------------------
  # 1. Draw Top Labels (p values)
  # ---------------------------------------------------------
  for(p_val in p.vec) {
    par(mar=c(0,0,0,0))
    plot(0,0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1))
    text(0, 0, paste0("p = ", p_val, "\n(GA vs Boot)"), cex=1.4, font=2)
  }

  # ---------------------------------------------------------
  # 2. Draw Left Labels (n values)
  # ---------------------------------------------------------
  for(n_val in n.vec) {
    par(mar=c(0,0,0,0))
    plot(0,0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1))
    text(0, 0, paste0("n =", n_val), cex=1.4, font=2, srt=90)
  }

  # ---------------------------------------------------------
  # 3. Draw the Plots (3 Rows x 6 Cols)
  # ---------------------------------------------------------
  # Reset margins for plotting
  par(mar=c(2, 2, 1, 1))

  for(n_val in n.vec) {
    for(p_val in p.vec) {

      # Load specific dataset
      d <- get_data(n_val, p_val, beta)

      if(!is.null(d)) {

        # --- Plot Left: Gaussian Approximation (Covariance) ---
        x_vals <- d$Z.max.norms
        y_vals <- d$Sigma.error.maxs

        # Calculate common range to ensure diagonal y=x line
        lims <- range(c(x_vals, y_vals), na.rm=TRUE)

        qqplot(x_vals, y_vals, pch=20, cex=0.8,
               xlim=lims, ylim=lims, # Force identical axes
               main="", xlab="", ylab="", axes=TRUE)
        abline(a=0, b=1, col="red", lwd=1.5)
        grid()
        box()

        # --- Plot Right: Bootstrap Validity (Covariance) ---
        x_vals <- d$Sigma.boot.maxs
        y_vals <- d$Sigma.error.maxs

        # Calculate common range
        lims <- range(c(x_vals, y_vals), na.rm=TRUE)

        qqplot(x_vals, y_vals, pch=20, cex=0.8,
               xlim=lims, ylim=lims, # Force identical axes
               main="", xlab="", ylab="", axes=TRUE)
        abline(a=0, b=1, col="red", lwd=1.5)
        grid()
        box()

      } else {
        # Placeholder if data is missing
        plot(1, type="n", axes=FALSE, main="")
        text(1, 1, "File Not Found")
        box()

        plot(1, type="n", axes=FALSE, main="")
        text(1, 1, "File Not Found")
        box()
      }
    }
  }

  dev.off()
  cat("Saved figure to", out.file, "\n")
}

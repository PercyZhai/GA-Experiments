# ---------------------------------------------------------
# Q-Q plots, low-dimensional simulation, precision matrix.
# For every (n, p) two panels are drawn: the empirical error max-norms
# (Omega.error.maxs) against the Gaussian approximation (Z.Omega.max.norms,
# "GA") and against the bootstrap (Omega.boot.maxs, "Boot"). One figure per beta.
#
# Paper:   toeplitz: Figs. 4-6 (beta = 2, 0.9, 0.55)
#          banded:   Figs. 21-23 (beta = 2, 0.9, 0.55), Appendix I
# Usage:   Rscript analysis/lowdim/qq_prec.R [toeplitz|banded]
# Inputs:  results/lowdim/<prefix>quick_<n>_<p>_<beta>_<setup>.RData
# Outputs: output/figures/<prefix>.QQ.prec.<beta>.pdf,
#          prefix = "E2" (toeplitz) or "E2B" (banded)
# ---------------------------------------------------------
if (!file.exists("functions/bootstrap_sample.R")) stop("Run this script from the root of the GA-Experiments repository.")

args <- commandArgs(trailingOnly = TRUE)
setup <- if (length(args) >= 1) args[1] else "toeplitz"
stopifnot(setup %in% c("toeplitz", "banded"))
prefix <- if (setup == "toeplitz") "E2" else "E2B"

# ---------------------------------------------------------
# Parameters
# ---------------------------------------------------------
n.vec <- c(10000, 2000, 200)
p.vec <- c(2, 10, 100)
beta.vec <- c(2, 0.9, 0.55)

base.path <- paste0("results/lowdim/", prefix, "quick_")  # E2quick_ (toeplitz) or E2Bquick_ (banded)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------
# Layout Setup
# ---------------------------------------------------------
# Layout Matrix Definition:
# 0  1  1  2  2  3  3   (Top Header Row: labels for p)
# 4  7  8  9 10 11 12   (Row 1: n=10000)
# 5 13 14 15 16 17 18   (Row 2: n=2000)
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
    warning("Missing input file: ", fname, immediate. = TRUE)
    return(NULL)
  }
}

for(beta in beta.vec) {

  pdf(paste0("output/figures/", prefix, ".QQ.prec.", beta, ".pdf"), width = 10, height = 5)

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

        # --- Plot Left: Gaussian Approximation (Precision) ---
        x_vals <- d$Z.Omega.max.norms
        y_vals <- d$Omega.error.maxs

        # Calculate common range to ensure diagonal y=x line
        lims <- range(c(x_vals, y_vals), na.rm=TRUE)

        qqplot(x_vals, y_vals, pch=20, cex=0.8,
               xlim=lims, ylim=lims, # Force identical axes
               main="", xlab="", ylab="", axes=TRUE)
        abline(a=0, b=1, col="red", lwd=1.5)
        grid()
        box()

        # --- Plot Right: Bootstrap Validity (Precision) ---
        x_vals <- d$Omega.boot.maxs
        y_vals <- d$Omega.error.maxs

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

  invisible(dev.off())
}

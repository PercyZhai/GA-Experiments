# ---------------------------------------------------------
# Empirical CDFs, low-dimensional simulation, covariance matrix.
# For every (n, p, beta) the ECDFs of the empirical error max-norms
# (Sigma.error.maxs), the Gaussian approximation (Z.max.norms) and the
# bootstrap (Sigma.boot.maxs) are overlaid, for beta = 2 and 0.9.
#
# Paper:   toeplitz: Fig. 14
#          banded:   Fig. 20, Appendix I
# Usage:   Rscript analysis/lowdim/ecdf_cov.R [toeplitz|banded]
# Inputs:  results/lowdim/<prefix>quick_<n>_<p>_<beta>_<setup>.RData
# Outputs: output/figures/<prefix>.ECDF.cov.pdf,
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
beta.vec <- c(2, 0.9)

base.path <- paste0("results/lowdim/", prefix, "quick_")  # E2quick_ (toeplitz) or E2Bquick_ (banded)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

pdf(paste0("output/figures/", prefix, ".ECDF.cov.pdf"), width = 10, height = 6)

# ---------------------------------------------------------
# Layout Setup
# ---------------------------------------------------------
# Header row (p labels), 3 rows of plots (n = 10000, 2000, 200), legend row
layout_mat <- matrix(c(
  0, 1, 1, 2, 2, 3, 3,
  4, 7, 8, 9, 10,11,12,
  5, 13,14,15,16,17,18,
  6, 19,20,21,22,23,24,
  25,25,25,25,25,25,25
), nrow=5, byrow=TRUE)

# Heights: Header=0.2, Plots=1, Legend=0.4
layout(layout_mat, widths=c(0.3, 1,1,1,1,1,1), heights=c(0.2, 1,1,1, 0.4))

# Minimize margins
par(mar=c(2, 2, 1, 1), oma=c(0,0,0,0))

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

# ---------------------------------------------------------
# 1. Draw Top Labels (p values)
# ---------------------------------------------------------
for(p_val in p.vec) {
  par(mar=c(0,0,0,0))
  plot(0,0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1))
  text(0, 0, paste0("p = ", p_val), cex=1.4, font=2)
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
# 3. Draw the ECDF Plots (3 Rows x 6 Cols)
# ---------------------------------------------------------
par(mar=c(2, 2, 2, 1))

for(n_val in n.vec) {
  for(p_val in p.vec) {
    for(beta_val in beta.vec) {

      d <- get_data(n_val, p_val, beta_val)

      if(!is.null(d)) {

        err_dat <- d$Sigma.error.maxs
        gaus_dat <- d$Z.max.norms
        boot_dat <- d$Sigma.boot.maxs

        full_range <- range(c(err_dat, gaus_dat, boot_dat), na.rm=TRUE)

        if(n_val == n.vec[1]) {
          plot_title <- paste0("Beta = ", beta_val)
        } else {
          plot_title <- ""
        }

        # Plot 1: Empirical (Black)
        plot(ecdf(err_dat),
             main = plot_title,
             col = "black",
             lwd = 2,
             verticals = TRUE,
             do.points = FALSE,
             xlim = full_range,
             xlab = "", ylab = "")

        # Plot 2: Gaussian (Red)
        lines(ecdf(gaus_dat),
              col = "red",
              lwd = 2,
              verticals = TRUE,
              do.points = FALSE)

        # Plot 3: Bootstrap (Blue)
        lines(ecdf(boot_dat),
              col = "blue",
              lwd = 2,
              lty = 1,
              verticals = TRUE,
              do.points = FALSE)

        grid()

      } else {
        plot(1, type="n", axes=FALSE, main="")
        text(1, 1, "Missing", cex=0.8)
        box()
      }
    }
  }
}

# ---------------------------------------------------------
# 4. Draw the Legend (Bottom Row)
# ---------------------------------------------------------
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend("center",
       legend = c("Empirical Error", "Gaussian Approx", "Bootstrap Approx"),
       col = c("black", "red", "blue"),
       lty = c(1, 1, 1),
       lwd = 3,
       cex = 1.4,           # Text size
       seg.len = 1.2,       # Shorter lines
       x.intersp = 0.5,     # Less space between line and text
       horiz = TRUE,
       bty = "n")

invisible(dev.off())

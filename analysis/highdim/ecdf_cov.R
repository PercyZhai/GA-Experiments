# ---------------------------------------------------------
# analysis/highdim/ecdf_cov.R
# ---------------------------------------------------------
# Description:
#   ECDF plots for the high-dimensional covariance experiment (E3).
#   Comparison: Empirical Error (Black) vs Gaussian Limit (Red) vs Bootstrap (Blue).
#   Layout: 3 Rows (n) x 6 Cols (p * beta), beta = 2, 0.9.
#
# Paper:
#   toeplitz: Fig. 16
#   banded:   Fig. 28 (Appendix I)
#
# Usage (from the repository root):
#   Rscript analysis/highdim/ecdf_cov.R [toeplitz|banded]   (default: toeplitz)
#
# Inputs:
#   results/highdim/<prefix>quick_<n>_<p>_<beta>_<setup>.RData
#   (written by simulation/highdim_sim.R; uses Sigma.error.maxs,
#   Z.max.norms, Sigma.boot.maxs)
#
# Outputs:
#   output/figures/<prefix>.ECDF.cov.pdf  (10 x 6 in),
#   prefix = "E3" (toeplitz) or "E3B" (banded)
# ---------------------------------------------------------
if (!file.exists("functions/bootstrap_sample.R")) stop("Run this script from the root of the GA-Experiments repository.")

args <- commandArgs(trailingOnly = TRUE)
setup <- if (length(args) >= 1) args[1] else "toeplitz"
stopifnot(setup %in% c("toeplitz", "banded"))
prefix <- if (setup == "toeplitz") "E3" else "E3B"

# 1. Parameters
n.vec <- c(800, 400, 200)
p.vec <- c(2000, 5000, 8000)
beta.vec <- c(2, 0.9)

# Input files (relative to the repository root)
base.path <- paste0("results/highdim/", prefix, "quick_")  # E3quick_ (toeplitz) or E3Bquick_ (banded)

# Output file
fig.dir <- "output/figures"
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
out.file <- paste0(fig.dir, "/", prefix, ".ECDF.cov.pdf")
pdf(out.file, width = 10, height = 6)

# ---------------------------------------------------------
# 2. Layout Setup
# ---------------------------------------------------------
# Grid:
# Row 1: Header (p-values) -> Spans 2 cols each
# Row 2: n=800 -> Label + 6 plots
# Row 3: n=400 -> Label + 6 plots
# Row 4: n=200 -> Label + 6 plots
# Row 5: Legend

layout_mat <- matrix(c(
  0, 1, 1, 2, 2, 3, 3,    # Header
  4, 7, 8, 9, 10,11,12,   # n = 800
  5, 13,14,15,16,17,18,   # n = 400
  6, 19,20,21,22,23,24,   # n = 200
  25,25,25,25,25,25,25    # Legend
), nrow=5, byrow=TRUE)

# Heights: Header=0.2, Plots=1, Legend=0.4
# Widths: Label=0.3, Plots=1 each
layout(layout_mat, widths=c(0.3, 1,1,1,1,1,1), heights=c(0.2, 1,1,1, 0.4))

# Minimize margins
par(mar=c(2, 2, 1, 1), oma=c(0,0,0,0))

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
    warning("Result file not found: ", fname, call. = FALSE, immediate. = TRUE)
    return(NULL)
  }
}

# ---------------------------------------------------------
# 4. Draw Top Labels (p values)
# ---------------------------------------------------------
# These span 2 columns each (for beta=2 and beta=0.9)
for(p_val in p.vec) {
  par(mar=c(0,0,0,0))
  plot(0,0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1))
  text(0, 0, paste0("p = ", p_val), cex=1.4, font=2)
}

# ---------------------------------------------------------
# 5. Draw Left Labels (n values)
# ---------------------------------------------------------
for(n_val in n.vec) {
  par(mar=c(0,0,0,0))
  plot(0,0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1))
  text(0, 0, paste0("n =", n_val), cex=1.4, font=2, srt=90)
}

# ---------------------------------------------------------
# 6. Draw the ECDF Plots (Row by Row)
# ---------------------------------------------------------
par(mar=c(2, 2, 2, 1))

for(n_val in n.vec) {
  for(p_val in p.vec) {
    for(beta_val in beta.vec) {

      d <- get_data(n_val, p_val, beta_val)

      if(!is.null(d)) {

        err_dat  <- d$Sigma.error.maxs
        gaus_dat <- d$Z.max.norms
        boot_dat <- d$Sigma.boot.maxs

        # Determine global range for x-axis to keep scales consistent within plot
        full_range <- range(c(err_dat, gaus_dat, boot_dat), na.rm=TRUE)

        # Title only on first row (n=800) to show Beta
        if(n_val == n.vec[1]) {
          plot_title <- paste0("Beta = ", beta_val)
        } else {
          plot_title <- ""
        }

        # 1. Empirical CDF (Black)
        plot(ecdf(err_dat),
             main = plot_title,
             col = "black",
             lwd = 2,
             verticals = TRUE,
             do.points = FALSE,
             xlim = full_range,
             xlab = "", ylab = "")

        # 2. Gaussian CDF (Red)
        lines(ecdf(gaus_dat),
              col = "red",
              lwd = 2,
              verticals = TRUE,
              do.points = FALSE)

        # 3. Bootstrap CDF (Blue)
        lines(ecdf(boot_dat),
              col = "blue",
              lwd = 2,
              lty = 1,
              verticals = TRUE,
              do.points = FALSE)

        grid()

      } else {
        # Placeholder for missing file
        plot(1, type="n", axes=FALSE, main="")
        text(1, 1, "Missing", cex=0.8)
        box()
      }
    }
  }
}

# ---------------------------------------------------------
# 7. Draw Legend (Bottom Row)
# ---------------------------------------------------------
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

legend("center",
       legend = c("Empirical Error", "Gaussian Approx", "Bootstrap Approx"),
       col = c("black", "red", "blue"),
       lty = c(1, 1, 1),
       lwd = 3,
       cex = 1.4,
       seg.len = 1.2,
       x.intersp = 0.5,
       horiz = TRUE,
       bty = "n")

invisible(dev.off())
cat("Saved figure to", out.file, "\n")

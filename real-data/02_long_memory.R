# Long-range dependence of the ABIDE NYU fMRI time series (paper Sec. VIII).
#
#   Fig. 10: histograms of the Hurst exponents of all time series (all subjects and
#            ROIs) in the ASD and in the control group
#   Fig. 11: ACF plots (lags 0-100) of four time series of the first ASD subject
#   Text:    share of time series with at least one significant ACF (|ACF| > 1.96/sqrt(n),
#            n = 176) at a lag between 21 and 100 (86.2% ASD, 85.9% control)
#
# Usage:   Rscript real-data/02_long_memory.R   (from the repository root)
# Inputs:  data/time.series.ABIDE.RData (written by real-data/01_read_abide.R)
# Outputs: output/figures/Hurst.ASD.pdf, output/figures/Hurst.ctrl.pdf,
#          output/figures/acf.1.pdf, ..., output/figures/acf.4.pdf;
#          the share of time series with H > 0.5 and the ACF statistic are printed

if (!file.exists("functions/bootstrap_sample.R")) stop("Run this script from the root of the GA-Experiments repository.")

library(pracma)
source("functions/load_abide.R")
load_abide() # h.1 (ASD), h.2 (control), n.1, n.2, p

fig.path <- "output/figures/"
dir.create(fig.path, recursive = TRUE, showWarnings = FALSE)

# Hurst exponent of one time series (simple R/S estimate Hs of pracma::hurstexp)
long.memory <- function(time_series) {
  # Check if the input is a numeric vector
  if (!is.numeric(time_series)) {
    stop("The input must be a numeric vector representing a time series.")
  }

  # Estimate the Hurst exponent using the 'hurstexp' function from pracma
  hurst_result <- hurstexp(time_series, display = FALSE)

  # Extract the estimated Hurst exponent H
  H <- hurst_result$Hs

  # Return the Hurst exponent as a measure of temporal dependence
  return(H)
}


# Fig. 10: Histogram of Hurst Exponent
# ASD group

H.Autism <- rep(NA, n.1 * p)
for(i in 1:n.1){
  for(j in 1:p){
    H.Autism[(i-1)*p+j] <- long.memory(h.1[i,j,])
  }
}
pdf(paste0(fig.path, "Hurst.ASD.pdf"), width = 6, height = 6)
hist(H.Autism, main="ASD Group", xlab="Hurst Exponent")
abline(v=0.5, col="red", lty=2, lwd=2)
invisible(dev.off())

# Control group

H.ctrl <- rep(NA, n.2 * p)
for(i in 1:n.2){
  for(j in 1:p){
    H.ctrl[(i-1)*p+j] <- long.memory(h.2[i,j,])
  }
}
pdf(paste0(fig.path, "Hurst.ctrl.pdf"), width = 6, height = 6)
hist(H.ctrl, main="Control Group", xlab="Hurst Exponent")
abline(v=0.5, col="red", lty=2, lwd=2)
invisible(dev.off())

cat("Share of time series with H > 0.5: ASD", round(mean(H.Autism > 0.5), 3),
    ", control", round(mean(H.ctrl > 0.5), 3), "\n")


# ACF (lags 0 to max_lag) of every time series [i, j, ] of a 3D array;
# column (i-1)*p + j holds subject i, ROI j
acf_matrix <- function(h, max_lag = 50) {
  # Dimensions of the 3D array
  n <- dim(h)[1]
  p <- dim(h)[2]

  # Initialize a matrix to store ACF values
  acf_values <- matrix(0, nrow = max_lag + 1, ncol = n * p)

  # Loop through each series (each [i, j, ] slice) and compute ACF
  count <- 1
  for (i in 1:n) {
    for (j in 1:p) {
      acf_values[, count] <- acf(h[i, j, ], plot = FALSE, lag.max = max_lag)$acf
      count <- count + 1
    }
  }

  return(acf_values)
}

acf.mat.1 <- acf_matrix(h.1, 100)
acf.mat.2 <- acf_matrix(h.2, 100)

# how many acfs have exceeded 1.96/sqrt(T) at some lag in 21:100 (rows 22:101; row 1 is lag 0),
# T = 176 time points
T.len <- dim(h.1)[3]
# ASD group: n.1 * p = 71 * 116 = 8236 time series
acf.thres.1 <- abs(acf.mat.1) > 1.96/sqrt(T.len)
lrd.1 <- 0
for(i in 1:ncol(acf.thres.1)){
  if(any(acf.thres.1[22:101, i] == 1)){
    lrd.1 <- lrd.1 + 1/ncol(acf.thres.1)
  }
}

# Control group: n.2 * p = 95 * 116 = 11020 time series
acf.thres.2 <- abs(acf.mat.2) > 1.96/sqrt(T.len)
lrd.2 <- 0
for(i in 1:ncol(acf.thres.2)){
  if(any(acf.thres.2[22:101, i] == 1)){
    lrd.2 <- lrd.2 + 1/ncol(acf.thres.2)
  }
}

cat("Share of time series with a significant ACF at a lag in 21-100: ASD",
    sprintf("%.1f%%", 100 * lrd.1), ", control", sprintf("%.1f%%", 100 * lrd.2), "\n")


# Fig. 11: ACF plots of ROIs 3, 52, 4 and 23 (AAL labels 2101, 5202, 2102, 2601)
# of the first ASD subject
acf.roi <- c(3, 52, 4, 23)
for(k in 1:4){
  pdf(paste0(fig.path, "acf.", k, ".pdf"), width = 6, height = 6)
  acf(h.1[1,acf.roi[k],], 100, main="")
  dev.off()
}

# Loading of the ABIDE NYU time series for the real-data analysis (paper Sec. VIII).
#
#   load_abide(file = "data/time.series.ABIDE.RData", envir = parent.frame())
#     Loads h.1 (ASD group) and h.2 (control group), arrays indexed
#     [subject, ROI, time] of dimension 73 x 116 x 176 and 98 x 116 x 176, and
#     removes ASD subjects 61, 71 and control subjects 1, 62, 78: each of them has
#     a constant ROI time series, so its sample covariance matrix is singular.
#     This leaves the 71 ASD subjects and 95 controls analysed in the paper.
#     Assigns h.1, h.2, n.1 (= 71), n.2 (= 95) and p (= 116) into envir.
#
# Used by real-data/02_long_memory.R and real-data/03_cig_bootstrap.R.
#
# Usage:   source("functions/load_abide.R"); load_abide()   (from the repository root)
# Inputs:  data/time.series.ABIDE.RData (written by real-data/01_read_abide.R)
# Outputs: none

load_abide <- function(file = "data/time.series.ABIDE.RData", envir = parent.frame()){
  if(!file.exists(file)) stop("Cannot find ", file, "; run real-data/01_read_abide.R first.")
  load(file) # h.1, h.2

  # remove subjects with a constant ROI series (singular sample covariance)
  h.1 <- h.1[-c(61,71),,]
  h.2 <- h.2[-c(1,62,78),,]

  assign("h.1", h.1, envir = envir)
  assign("h.2", h.2, envir = envir)
  assign("n.1", dim(h.1)[1], envir = envir)
  assign("n.2", dim(h.2)[1], envir = envir)
  assign("p", dim(h.1)[2], envir = envir)
  invisible(NULL)
}

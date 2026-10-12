# Brain connectome of the ABIDE NYU data by block bootstrap (paper Sec. VIII).
#
# For every subject, entrywise block-bootstrap confidence intervals (l = floor(n^(2/3)) = 31,
# alpha = 0.05) for the precision matrix give an individual conditional independence
# graph (CIG); the individual graphs are averaged within each group and the top 2%
# of the off-diagonal entries are kept.
#
#   Fig. 12: aggregated CIGs of the ASD group (Autism_02_95.pdf) and of the control
#            group (Control_02_95.pdf)
#   Text:    change of the number of connections per ROI from the control to the
#            ASD group (e.g. 2311: 6 -> 1, 5021: 1 -> 4)
#
# Usage:   Rscript real-data/03_cig_bootstrap.R   (from the repository root)
# Inputs:  data/time.series.ABIDE.RData (written by real-data/01_read_abide.R)
# Outputs: output/figures/Autism_02_95.pdf, output/figures/Control_02_95.pdf;
#          the degree-change table is printed

if (!file.exists("functions/bootstrap_sample.R")) stop("Run this script from the root of the GA-Experiments repository.")

source("functions/load_abide.R")
source("functions/bootstrap_ci.R")
source("functions/graph_plotter.R") # also loads igraph and coordinate.list, name.tag.vec
load_abide() # h.1 (ASD), h.2 (control), n.1, n.2, p

fig.path <- "output/figures/"
dir.create(fig.path, recursive = TRUE, showWarnings = FALSE)

# Analyzing h.1 (Autism group)
# Autism.graph[j,k] = share of subjects whose CI for Omega_jk does not contain zero
cat("ASD group:", n.1, "subjects\n")
Autism.graph <- matrix(0,p,p)
for(i in 1:n.1){
  X <- t(h.1[i,,]) # Now X is a 176 * 116 matrix
  n <- nrow(X)
  p <- ncol(X)

  l <- as.integer(n^(2/3))

  CI.boot <- prec.CI.boot(X, l, alpha=0.05)
  graph <- CI.boot[[1]] < 0 | CI.boot[[2]] > 0

  Autism.graph <- Autism.graph + graph / n.1
}

# Analyzing h.2 (Control group)
cat("Control group:", n.2, "subjects\n")
CTRL.graph <- matrix(0,p,p)
for(i in 1:n.2){
  X <- t(h.2[i,,]) # Now X is a 176 * 116 matrix
  n <- nrow(X)
  p <- ncol(X)

  l <- as.integer(n^(2/3))

  CI.boot <- prec.CI.boot(X, l, alpha=0.05)
  graph <- CI.boot[[1]] < 0 | CI.boot[[2]] > 0

  CTRL.graph <- CTRL.graph + graph / n.2
}

Autism.graph.undiag <- Autism.graph - diag(diag(Autism.graph))
CTRL.graph.undiag <- CTRL.graph - diag(diag(CTRL.graph))

# Keep the edges whose share is strictly larger than the 0.02*p*(p-1)-th largest
# off-diagonal value (sparsity at most 2%; entries tied with the cutoff are dropped)
Autism.thres <- Autism.graph.undiag > sort(Autism.graph.undiag, decreasing=T)[0.02 * p*(p-1)]
CTRL.thres <- CTRL.graph.undiag > sort(CTRL.graph.undiag, decreasing=T)[0.02 * p*(p-1)]

cat("Graph sparsity: ASD", round(sum(Autism.thres) / (p*(p-1)), 4),
    ", control", round(sum(CTRL.thres) / (p*(p-1)), 4), "\n")

# Fig. 12
pdf(paste0(fig.path, "Autism_02_95.pdf"), width = 5, height = 5)
invisible(graph.plotter(Autism.thres, coordinate.list, name.tag.vec, "yz", "Autism", "blue", 0.8))
invisible(dev.off())
pdf(paste0(fig.path, "Control_02_95.pdf"), width = 5, height = 5)
invisible(graph.plotter(CTRL.thres, coordinate.list, name.tag.vec, "yz", "Control", "black", 0.8))
invisible(dev.off())


# Number of connections of each ROI and its change from the control to the ASD group
Autism.connect <- apply(Autism.thres, 2, sum)
CTRL.connect <- apply(CTRL.thres, 2, sum)
node.change <- Autism.connect - CTRL.connect

degree.change <- data.frame(ROI = name.tag.vec, ASD = Autism.connect,
                            Control = CTRL.connect, change = node.change)
degree.change <- degree.change[node.change != 0, ]
degree.change <- degree.change[order(degree.change$change), ]
cat("ROIs whose number of connections changed (ASD - control), most decreased first:\n")
print(degree.change, row.names = FALSE)

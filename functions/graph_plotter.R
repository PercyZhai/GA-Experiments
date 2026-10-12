# Plotting of conditional independence graphs (CIGs) over the 116 AAL ROIs.
#
# Used for Fig. 12 (Autism_02_95.pdf, Control_02_95.pdf) by real-data/03_cig_bootstrap.R.
# Loads igraph and sources functions/aal_info.R (coordinate.list, name.tag.vec).
#
#   roi_layout(coords, option)  2-D node coordinates of the ROIs: projection of
#                               coords (coordinate.list) onto the "xy", "xz" or "yz" plane
#   graph.plotter(adj.mat, coordinate.list, name.tag.vec, option, title, color, curve)
#                             draws the graph with adjacency matrix adj.mat on the
#                             current device; ROIs are placed by their (y, z)
#                             coordinates (sagittal view) and the node size grows
#                             with the node degree. 'option' and 'title' are not used.
#
# Usage:   source("functions/graph_plotter.R")   (from the repository root)
# Inputs:  none
# Outputs: none (defines roi_layout() and graph.plotter())

library(igraph)
source("functions/aal_info.R")

# Named roi_layout (not layout) so that graphics::layout is not masked.
roi_layout <- function(coords, option=c("xy","xz","yz")){
  option <- match.arg(option)
  p <- length(coords)
  result <- matrix(NA, nrow=p, ncol=2)
  coord.vec <- unlist(coords) # length is 3p
  x.vec <- rep(NA, p)
  y.vec <- rep(NA, p)
  z.vec <- rep(NA, p)

  for(j in 1:p){
    x.vec[j] <- coord.vec[3*(j-1) + 1]
    y.vec[j] <- coord.vec[3*(j-1) + 2]
    z.vec[j] <- coord.vec[3*(j-1) + 3]
  }

  if(option=="xy"){
    result[,1] <- x.vec
    result[,2] <- y.vec
  }
  if(option=="xz"){
    result[,1] <- x.vec
    result[,2] <- z.vec
  }
  if(option=="yz"){
    result[,1] <- y.vec
    result[,2] <- z.vec
  }

  return(result)
}

graph.plotter <- function(adj.mat, coordinate.list, name.tag.vec,
                          option=c("xy","xz","yz"), title, color, curve=0.8){
  # compute number of neighbors of each node
  nbr.count <- colSums(adj.mat)

  net <- graph_from_adjacency_matrix(adj.mat, mode="undirected") %>%
    set_vertex_attr("label", value = name.tag.vec)

  plot(net, vertex.label=name.tag.vec, vertex.shape = "circle",
       vertex.color = "red", vertex.frame.color = "red",
       vertex.label.font = 1, vertex.label.family = "Helvetica",
       vertex.label.cex = 0.3, vertex.label.dist = 0.4, vertex.label.degree = pi/4,
       edge.color=color,  edge.width = 1.2, vertex.size = 0.8*(nbr.count+1),
       layout=roi_layout(coordinate.list, "yz"), edge.curved = curve)
  return(0)
}

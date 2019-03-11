#' Get Overlaping Modules
#'
#' @keywords internal
#'
#' @import data.table
#'
#' @return
getToMerge <- function(nodes) {
  overlapping <- unique(nodes[nodes, on = .(nodeID), allow.cartesian = TRUE][, .(moduleID, i.moduleID)])
  overlapGraph <- graph_from_edgelist(as.matrix(overlapping), directed = FALSE)
  cluster <- clusters(overlapGraph)
  return(data.table(cluster = cluster$membership, moduleID = V(overlapGraph)))
}

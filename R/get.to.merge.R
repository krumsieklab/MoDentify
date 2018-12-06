#' Get Overlaping Modules
#' 
#' @keywords internal
#'
#' @import data.table
#'
#' @return
getToMerge<-function(nodes){
  overlapping<-unique(nodes[nodes, on=.(nodeID), allow.cartesian=TRUE][,.(moduleID, i.moduleID)])
  overlap_graph<-graph_from_edgelist(as.matrix(overlapping), directed=F)
  cluster=clusters(overlap_graph)
  return(data.table(cluster=cluster$membership, moduleID=V(overlap_graph)))
}

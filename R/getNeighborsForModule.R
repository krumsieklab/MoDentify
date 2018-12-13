#' Get Neighbors for Module
#' 
#' @keywords internal
#'
#' @param graph an \code{\link[igraph]{igraph}} object
#' @param module a vector containing nodes corresponding to the graph
#'
#' @import igraph
#' @export
#' @return a vector containing all neighbors of the module
#' @examples
#' data(qmdiab.data)
#' data(qmdiab.annos)
#' 
#' net.graph<-generateNetwork(data=qmdiab.data, annotations=qmdiab.annos)
#' module<-c(3,4)
#' neighbors<-getNeighborsForModule(graph = net.graph, module = module)
getNeighborsForModule<-function(graph, module){
  neighbors<-unique(unname(unlist(lapply(module, neighbors, graph=graph))))
  
  return(neighbors[!neighbors %in% module])
}

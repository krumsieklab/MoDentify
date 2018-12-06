#' Delete Duplicates from Modules
#' 
#' @keywords internal
#'
#' @param nodes
#' @param modules
#' @param graph
#'
#' @import data.table
#' @import igraph
#'
#' @return

deleteDuplicates<-function(nodes, modules, graph){
  DT<-modules[nodes, on=.(moduleID)]
  DTUnique<-DT
  moduleIDs<-modules$moduleID
  for(i in 1:(length(moduleIDs)-1)){
    for(j in (i+1):length(moduleIDs)){
      v1<-V(graph)[DT[ moduleID == moduleIDs[i],nodeID]]
      v2 <- V(graph)[DT[ moduleID == moduleIDs[j],nodeID]]
      if(length(difference(v1, v2)) == 0 & length(difference(v2, v1)) ==0){
        DTUnique<-DTUnique[moduleID !=moduleIDs[j]]
      }
    }
  }
  return(DTUnique)
}

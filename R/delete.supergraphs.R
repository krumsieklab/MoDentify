#' Delete Supergraphs from Modules
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

delete.supergraphs<-function(nodes, modules, graph){
  DT<-modules[nodes, on=.(moduleID)]
  DT_max<-DT
  moduleIDs<-modules$moduleID
  for(i in 1:(length(moduleIDs))){
    for(j in (i):length(moduleIDs)){

      if(i != j){
        v1<-sort(V(graph)[DT[ moduleID == moduleIDs[i],nodeID]])
        v2<-sort(V(graph)[DT[ moduleID == moduleIDs[j],nodeID]])

        intersect<-sort(igraph::intersection(v1, v2))

        if(length(intersect) == length(v1) && all(intersect == v1)){
          if(modules[moduleID == moduleIDs[i], module.score] >=
             modules[moduleID == moduleIDs[j], module.score]){
            DT_max<-DT_max[moduleID !=moduleIDs[j]]
          }else{
            DT_max<-DT_max[moduleID !=moduleIDs[i]]
          }
        }else if(length(intersect) == length(v2) && all(intersect == v2)){
          if(modules[moduleID ==moduleIDs[j], module.score] >=
             modules[moduleID ==moduleIDs[i], module.score]){
            DT_max<-DT_max[moduleID !=moduleIDs[i]]
          }else{
            DT_max<-DT_max[moduleID !=moduleIDs[j]]
          }
        }
      }
    }
  }
  DT_max$i.module.beta <- NULL
  DT_max$i.adjusted.score <- NULL
  return(DT_max)
}

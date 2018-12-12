#' Plotting of the edges between attribute types
#'
#'
#' @param graph an \code{\link[igraph]{igraph}} object
#' @param name the name of the attribute
#' @param rm.unknown remove unknowns
#' @import ggplot2
#' @export
#' @import data.table
#' @import igraph
#' @importFrom stats setNames
#' @examples
#' data(qmdiab.data)
#' data(qmdiab.annos)
#' 
#' net.graph<-generateNetwork(data = qmdiab.data, annotations = qmdiab.annos)
#' ggEdgesBetweenAttributes(net.graph, "Super.pathway")
ggEdgesBetweenAttributes<-function(graph, name, rm.unknown=FALSE){
  l<-lapply(E(graph), function(x)
    vertex_attr(graph= graph, name = name, index=as.vector(ends(graph = graph, es = x))))

  DT<-as.data.table(setNames(do.call(rbind.data.frame, l), c("att1", "att2")))
  
  if (rm.unknown){
    DT<-DT[!grepl("Unknown", DT$att1),]
    DT<-DT[!grepl("Unknown", DT$att2),]}
  
  DT[, att1:=as.character(att1)]
  DT[, att2:=as.character(att2)]
  DT<-rbind(DT,DT[att1 != att2, .(att1=att2, att2=att1)])
  p<-ggplot(DT, aes(att1, att2))+ geom_bin2d()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dat<-ggplot_build(p)$data[[1]]
  p+geom_text(data=dat, aes((xmin + xmax)/2, (ymin + ymax)/2, label=count), col="white", size=3)+ylab(name)+xlab(name)
}

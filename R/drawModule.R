#' Draw Single Module
#'
#'
#' @param moduleNR number of the module, that should be drawn.
#' @param graph an \code{\link[igraph]{igraph}} object, which can be generated with
#' \code{\link{generateNetwork}}.
#' @param title a name for the Cytoscape session.
#' @param nodes \code{\link[data.table]{data.table}} with data containing the nodes.
#' It is contained in the output list of \code{\link{identifyModules}}.
#' @param colors a colour palette as returned from \code{\link[grDevices]{rainbow}} 
#' for colouring the different modules.
#' @param save.image TRUE, if the modules should be saved as png files.
#' @param close.cycnets.afterwards TRUE, if the windows in the cytoscape environment should be closed
#' after drawing. This might be useful repeated function call.
#'
#' @import RCy3
#' @import data.table
#' @import igraph
#' @keywords internal
#' @references
#' \insertRef{Shannon2003}{MoDentify}
#' @references
#' \insertRef{Smoot2011}{MoDentify}
#' @references
#' \insertRef{RCy3}{MoDentify}

drawModule<-function(moduleNR, graph, title="", nodes, colors, save.image=TRUE,
                     close.cycnets.afterwards = FALSE){
    
    
    
   
  nodes<-nodes[moduleID==moduleNR]
  module_graph<-induced_subgraph(graph, nodes$nodeID)
  for(i in seq_len(dim(nodes)[1])){
    if("Fluid" %in% vertex_attr_names(graph)){
      module_graph<-set_vertex_attr(module_graph, "label", nodes[i]$name,
                                    paste0(
                                      toupper(substr(vertex_attr(graph, "Fluid" ,nodes[i]$name),0,1)),
                                      "::", nodes[i]$label))
    }else{
      module_graph<-set_vertex_attr(module_graph, "label", nodes[i]$name, nodes[i]$label)
    }

    module_graph<-set_vertex_attr(module_graph, "module.name", nodes[i]$name, moduleNR)
    module_graph<-set_vertex_attr(module_graph, "p.value", nodes[i]$name, nodes[i]$node.pval)
    module_graph<-set_vertex_attr(module_graph, "is.significant", nodes[i]$name, nodes[i]$is.significant)
  }


  netSUID <- createNetworkFromIgraph(igraph = module_graph,title = 
                                         paste0(title, "_", "module", moduleNR))
 
  setNodeColorMapping("module.name", names(colors), c(substr(colors, 0, 7)), mapping.type = "d")  
  setNodeSizeMapping("p.value", c(0.05/vcount(graph), 0.1), c(100, 75, 40, 25), default.size = c(20))
  setNodeShapeMapping("is.significant", c("TRUE", "FALSE"), c("diamond", "ellipse"))
  layoutNetwork(layout.name = "kamada-kawai")
  
  if(save.image){
      exportImage(filename =  paste0(getwd(), "/", title, "_", "module", moduleNR, ".png"),
                  type = "png", resolution = 600, height = 800)
  }

  
  if(close.cycnets.afterwards){
      closeSession()
  }

}

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
#'
#' @import data.table
#' @import igraph
#' @keywords internal
#' @references
#' \insertRef{Shannon2003}{MoDentify}
#' @references
#' \insertRef{Smoot2011}{MoDentify}
#' @references
#' \insertRef{RCy3}{MoDentify}

drawModule<-function(moduleNR, graph, title="", nodes, colors, save.image=TRUE){
    
    if (!requireNamespace("RCy3", quietly=TRUE)){
        stop("drawModule() requires 'RCy3' package")
    }
    
   
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
 
  
  
  
  # setNodeColorRule(cw,"module.name" , names(colors), c(substr(colors, 0, 7)), mode="lookup", default.color='#D8D0CE')
  # 
  # setWindowSize (cw, 1200, 1200)
  # layoutNetwork (cw, "kamada-kawai")
  # 
  # 
  # if(save.image){
  #   saveImage(cw, paste0(getwd(), "/", title, "_", "module", moduleNR, ".png"), "png", 1.0)
  # }


}

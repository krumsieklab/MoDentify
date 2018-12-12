#' Module Identification
#'
#'
#' @param graph an \code{\link[igraph]{igraph}} object, which can be generated with \code{\link{generateNetwork}}.
#' The ID of the nodes must correspond to the name of the variables.
#' @param data either a matrix, where the columns correspond to the variables and the rows to the observations.
#' Or a \code{\link[data.table]{data.table}} with three columns: name, sampleID and value.
#' @param phenotype a vector with the values for a phenotype of interest.
#' It must have the same number of observations as data.
#' @param covars a \code{\link[data.table]{data.table}} containing the covariates
#' to correct for.
#' The rows for the samples must be in the same order as in the phenotype vector.
#' @param annotations a \code{\link[data.table]{data.table}} containing annotations for the variables.
#' The columns correspond to the different annotations, the rows to the variables.
#' @param merge.overlapping if \code{TRUE}, overlapping modules will be merged.
#' @param better.than.components if \code{TRUE}, modules will only be enlarged and accepted,
#' if they are better than all of their components.
#' @param alpha significance level for accepting the modules.
#' @param level Must be set to the name of the column to be used, if modules should be calculated for pathways.
#' @param pathway.column
#' @param representative.method the method, that is used for the calculation of the module representation.
#' Currently implemented: "eigenmetabolite" and "average"
#' @param correction.method the method that used for multiple testing correction ("bonferroni", "BH", "BY", "fdr", "holm", "hochberg", "hommel", "none").
#' Default is set to bonferroni. See \code{\link[stats]{p.adjust}}.
#'
#' @references
#' \insertRef{Do2017}{MoDentify}
#' @references
#' \insertRef{Chuang2007}{MoDentify}
#' @import data.table
#' @import igraph
#' @export identifyModules
#' @usage identifyModules(graph, data, phenotype, covars = NULL, annotations,
#' merge.overlapping=FALSE, better.than.components= TRUE, alpha=0.05,
#' level=NULL, representative.method="average", correction.method="bonferroni")
#' @return a list consisting of four elements.
#' @examples
#' data(qmdiab.data)
#' data(qmdiab.annos)
#' data(qmdiab.phenos)
#' 
#' data<-qmdiab.data[, 1:75]
#' annotations<-qmdiab.annos[1:75]
#'
#' net.graph<-generateNetwork(data=data, annotations=annotations)
#' mods<-identifyModules(graph=net.graph, data=data, annotations = annotations
#' , phenotype = qmdiab.phenos$T2D, alpha = 0.05)
#'
#' pathway.graph<-generatePathwaysNetwork(data=data, annotations=annotations)
#' 
#' pathway.modules<-identifyModules(graph=pathway.graph$network, data=data,
#' phenotype = qmdiab.phenos$T2D, level = pathway.graph$level, annotations = annotations,
#' alpha = 0.05)
identifyModules<-function(graph, data, phenotype, covars = NULL,
                          annotations,
                          merge.overlapping = FALSE,
                          better.than.components = TRUE,
                          alpha = 0.05,
                          level = NULL,
                          representative.method = "average",
                          correction.method = "bonferroni"){
  
  
  # for incomplete data only average approach possible
  if (sum(is.na(data))>0){
    if(representative.method=="eigenmetabolite"){
      stop("Data matrix contains missing values.\n Module identification with eigenmetabolite approach not possible.\n")
    }else{
      warning("Data matrix contains missing values.\n Only complete cases were used for module representatives (average).\n")
    }
  }


  if(!is.data.table(data)){
    if(is.data.frame(data)){
      data<-as.data.table(data)
      data <- cbind(sampleID = paste0("sample", seq_len(dim(data)[1])), data)
    }else if(is.matrix(data)){
      data<-data.table(sampleID= paste0("sample", seq_len(dim(data)[1])), data)
    }
    data<-melt(data=data, id.vars = "sampleID", variable.name = "name")
  }

  annotations <- as.data.table(annotations)
  
  #Calculate the z-scores for variables
  data[, z.score:=scale(value), by=.(name)]
  if(!is.null(level)){

    if(!"Fluid" %in% colnames(annotations)){
      annotations[, s.id:=annotations[, get(level)]]
    }else{
      annotations[, s.id:=paste(annotations[, get(level)], Fluid)]
    }

    data<-data[annotations[, .(name, s.id)], on="name"]
    data[, met.name:=name]
    data[, name:= s.id]
    data[, s.id := NULL]
  }


  if(!is.null(covars)){
    if(dim(covars)[1]!=length(phenotype)){stop("Covars and Phenotype have a different amount of samples")}
    if(dim(covars)[1]!=length(unique(data$sampleID))){
      stop("Variables and covariates have a different number of observations.")}
  }else{
    if(length(phenotype)!=length(unique(data$sampleID))){
      stop("Variables and covariates have a different number of observations.")}
  }


  if(!all(as.character(unique(data$name)) %in% vertex_attr(graph, "name"))){
    warning("Not all of your variables are represented in the graph.")
  }
  if(!all(vertex_attr(graph, "name") %in% as.character(unique(data$name)))){
    stop("Not all of your nodes are represented in the data data.table")
  }
  
  message("The function identifyModules could take a few minutes.")

  modules<-data.table(moduleID=integer(), module.score=numeric(),
                      module.beta=numeric(), adjusted.score=numeric())
  nodes<-data.table(moduleID=integer(), nodeID=integer(), name=character(), label=character(),
                      order.added=integer(), score.after.adding=numeric())
  already_calculated<-data.table(key.value=character(), score=numeric(), beta=numeric(),
                                 times.accessed=numeric())
  seed.scores<-c()
  seed.betas<-c()    
  
  message("Proceeding node... ")

  for(v in V(graph)){
    #message(cat(paste0(v, " ")))
    cat(paste(v," "))
    module<-greedyModuleSelection(graph, v, data, phenotype, covars, alpha,
                                    already_calculated, better.than.components, representative.method=representative.method)
    already_calculated<-module$already_calculated
    seed.scores<-c(seed.scores, module$seed.score)
    seed.betas<-c(seed.betas, unname(module$seed.beta))

    adjusted.score<-p.adjust(p = module$module.score, n = vcount(graph), method = correction.method)
    if(length(module$module) > 1 &  adjusted.score < alpha){
      newID<-dim(modules)[1]+1

      modules<-rbind(modules, data.table(moduleID=newID, module.score=module$module.score,
                                         module.beta=module$beta, adjusted.score=adjusted.score))

      nodes<-rbind(nodes, data.table(moduleID=newID, nodeID=module$module,
                                         name=vertex_attr(graph, "name", module$module),
                                         label=vertex_attr(graph, "label", module$module),
                                         order.added=seq_len(length(module$module)),
                                         score.after.adding=module$score.sequence))

    }
  }
  seed_scores_DT<-data.table(nodeID=V(graph), seed.score=seed.scores, seed.beta=seed.betas)
  
  if (dim(modules)[1]>0){

    message("\nDeleting duplicate modules.")
    tmp_DT<-deleteDuplicates(nodes, modules, graph)
    modules<-unique(tmp_DT[, .(moduleID, module.score, module.beta, adjusted.score)])
    tmp_DT[,module.score := NULL]
    nodes<-tmp_DT
  
    if(better.than.components){
      message("Deleting less significant contained modules.")
      tmp_DT<-deleteSupergraphs(nodes, modules, graph)
      modules<-unique(tmp_DT[, .(moduleID, module.score, module.beta, adjusted.score)])
      tmp_DT[,module.score := NULL]
      nodes<-tmp_DT
    }
  
    if(merge.overlapping){
      message("Merging overlapping modules.")
      toMerge<-getToMerge(nodes)
      nodes<-toMerge[nodes, on=.(moduleID)]
      nodes[, moduleID := NULL]
      nodes[, module.beta:= NULL]
      nodes[, adjusted.score:= NULL]
      nodes[, adjusted.pval:= NULL]
      setnames(nodes, c("moduleID", "nodeID", "name", "label", "order.added", "score.after.adding"))
      modules<-getMergedModules(graph, data, phenotype, covars, nodes)
      modules$adjusted.score<-p.adjust(p = modules$module.score, n = vcount(graph), method = correction.method)
    }
  
    nodes<-nodes[seed_scores_DT, on=.(nodeID)]
    nodes[, name:=vertex_attr(graph, "name", nodeID)]
    nodes[, label:=vertex_attr(graph, "label", nodeID)]
    message(paste0("Number of modules found: ", dim(modules)[1]))
  }else(message("No modules found."))
  return(list(modules=modules, nodes=nodes, seeds=seed_scores_DT, cache=already_calculated))
}

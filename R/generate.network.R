#' Network inference with Gaussian graphical models (GGMs) or Pearson correlation
#'
#'
#' @param data a \code{\link[data.table]{data.table}} or matrix containing the data.
#' The columns correspond to the variables (e.g. metabolites), the rows to the observations.
#' @param covars a \code{\link[data.table]{data.table}} containing covariates to correct for.
#' The columns correspond to the different covariates, the rows to the observations.
#' @param annotations a \code{\link[data.table]{data.table}} containing annotations for the variables (e.g. pathway annotations).
#' The columns correspond to the different annotations, the rows to the variables
#' @param correlation.type type of correlation to be estimated. Can either be "pearson", or "partial".
#' @param alpha significance level (type 1 error) for multiple testing correction.
#' @param correction.method the method that should be used for multiple testing correction ("bonferroni", "BH", "BY", "fdr", "holm", "hochberg", "hommel", "none").
#' Default is bonferroni. See \code{\link[stats]{p.adjust}}.
#' @import GeneNet
#' @import data.table
#' @import igraph
#' @importFrom Hmisc rcorr
#' @importFrom stats complete.cases
#' @importFrom stats p.adjust
#' @export
#' @examples
#' net.graph<-generate.network(data = qmdiab.data,
#' annotations = qmdiab.annos, alpha=0.05, correction.method = "bonferroni")
#' @return a network containing the variables as nodes as an \code{\link[igraph]{igraph}} object.
#' @references
#' \insertRef{Krumsiek2011}{MoDentify}

generate.network <- function(data, covars = NULL, annotations, correlation.type = "partial", alpha=0.05,
                           correction.method = "bonferroni"){

  
  # warning that there are missing values and that only complete cases are taken
  if (sum(is.na(data))>0){
    warning("Data matrix contains missing values.\n Only complete cases were used for network inference.\n")
  }
  
  # calculating Pearson correlation
  C_pear<-rcorr(as.matrix(data), type = "pearson")
  colnames(C_pear$P)<-1:dim(C_pear$P)[1]
  rownames(C_pear$P)<-1:dim(C_pear$P)[1]
  
  pearp <- C_pear$P
  pearp[lower.tri(pearp,diag = T)] <- NA
  
  net_DT<-data.table(node1= as.integer(rownames(pearp)), pearp)
  net_DT<-melt(net_DT, id.vars = "node1", variable.name = "node2", value.name = "p")
  net_DT[, node2 := as.integer(node2)]
  net_DT[, pearp.adjust:=p.adjust(p, method = correction.method, n=choose(dim(data)[2], 2))]
  
  # partial correlation (GGM) or only Pearson correlation
  if(correlation.type=="partial"){
    
    
    # combining data and covariates
    if(!is.null(covars)){
      mat<-cbind(data, covars)
    }else{
      mat<-as.matrix(data)
    }
    mat<-as.matrix(mat[complete.cases(mat),])
    
    ggm<-ggm.estimate.pcor(mat, method = "static")
    ggm<-ggm[1:dim(data)[2], 1:dim(data)[2]]
    ggm_pvals <- as.data.table(network.test.edges(ggm, direct=FALSE, fdr=TRUE, plot=FALSE))
    ggm_pvals[, c("qval", "prob") := NULL]
    
    # only get complete cases
    net_DT<-net_DT[complete.cases(net_DT)]
    
    net_DT<-net_DT[ggm_pvals, on=c("node1", "node2")]
    
    net_DT$p<-net_DT$pval
    net_DT[, pval.adjust:=p.adjust(pval, method = correction.method, n=choose(dim(data)[2], 2))]
    names(net_DT)[5] <- "cor"
    net_DT[,"pval" := NULL]
    
    
  }else{
    C_pear_DT<-data.table(node1= as.integer(rownames(C_pear$P)), C_pear$r)
    colnames(C_pear_DT) <- c("node1", colnames(C_pear$P))
    C_pear_DT<-melt(C_pear_DT, id.vars = "node1", variable.name = "node2", value.name = "cor")
    C_pear_DT[, node2 := as.integer(node2)]
    
    
    net_DT<-net_DT[C_pear_DT, on=c("node1", "node2")]
    net_DT$pval.adjust <- rep(0,dim(net_DT)[1])
    # only get complete cases
    net_DT<-net_DT[complete.cases(net_DT)]
    
    }

  # assigning variable labels
  net_DT[, from:=colnames(C_pear$r)[node1]]
  net_DT[, to:=colnames(C_pear$r)[node2]]
  
  # in annotations?
  net_DT<-net_DT[to %in% annotations$name & from %in% annotations$name]
  
  # significance filtering
  net_DT<-net_DT[pval.adjust <= alpha & pearp.adjust <= alpha]
  
  
  net_DT[, c("node1", "node2", "pearp.adjust", "pval.adjust") := NULL]
  setcolorder(net_DT, c("from", "to", "cor", "p"))



  if("Labels" %in% colnames(annotations)){
    annotations[, label:=Labels]
  }
  
  if(!"label" %in% colnames(annotations)){
    annotations$label <- annotations$name
  }
  
  # make sure that column name is first column in annotations
  tmpnames <- annotations$name
  annotations$name <- NULL
  annotations <- cbind(name=tmpnames, annotations)

  net_graph<-graph_from_data_frame(d=net_DT, directed = F, vertices = annotations)

  return(net_graph)
}

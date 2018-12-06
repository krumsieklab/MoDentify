#' Subnetwork Creation function
#'
#'
#' @param data a \code{\link[data.table]{data.table}} or matrix containing data.
#' The columns correspond to the variables, the rows to the observations.
#' @param covars a \code{\link[data.table]{data.table}} containing covariates to correct for.
#' The columns correspond to the different covariates, the rows to the observations.
#' @param annotations  a \code{\link[data.table]{data.table}} containing annotations for the variables.
#' The columns correspond to the different annotations, the rows to the variables.
#' @param level name of the column that should be used for grouping the variables. The default is "SUB_PATHWAY".
#' @param alpha significance threshold (type 1 error) for multiple testing correction.
#' @param correction.method the method that should be used for multiple testing correction ("bonferroni", "BH", "BY", "fdr", "holm", "hochberg", "hommel", "none").
#' Default is bonferroni. See \code{\link[stats]{p.adjust}}.
#' @param rm.unknown remove variables, for which the group is unknown. The corresponding grouping variable name should be "Unknown".
#' Default is \code{TRUE}.
#' @param representative.method the method, that is used for calculation of the eigenmetabolites.
#' Currently implemented: "eigenmetabolite" and "average"
#' @import igraph
#' @export
#' @examples
#' pathway.graph<-generate.pathways.network(data = qmdiab.data,
#'  annotations = qmdiab.annos, level = "Sub.pathway")
#' @return a list with a network containing the variables as nodes as
#' an \code{\link[igraph]{igraph}} object and a list with the module representatives
#' including their explained variances.
#' @references
#' \insertRef{Krumsiek2011}{MoDentify}
#' @references
#' \insertRef{Do2017}{MoDentify}
generate.pathways.network <- function(data, covars=NULL, annotations,
                         level="Sub.pathway", 
                         correlation.type="partial",
                         alpha=0.05,
                         correction.method = "bonferroni", rm.unknown=TRUE,
                         representative.method = "eigenmetabolite"){
  
  
  annotations <- as.data.table(annotations)
  annotations<-order.annotation(annotations, colnames(data))

  data<-scale(data)
  subannotations<-copy(annotations)
  subannotations[,label:=get(level)]
  if(!"Fluid" %in% colnames(annotations)){
    subannotations[, name:=subannotations[, get(level)]]
  }else{
    subannotations[, name:=paste(subannotations[, get(level)], Fluid)]
  }
  
  if(representative.method=="eigenmetabolite"){
    if (sum(is.na(data))>0){
      stop("Data matrix contains missing values.\n Pathway network via eigenemtabolite approach not possible.\n")
    }else{
      eigen.metabolites<-eigen.metabolites(data=data, group=subannotations[, name], method = "PCA")
      mat<-as.matrix(eigen.metabolites$M)
      representatives <- eigen.metabolites
    }
  
  }else if(representative.method=="average"){
    representatives <- NULL
    mat<-as.matrix(aggregated.mean(data=data, group=subannotations[, name]))
    representatives$M <- mat
    representatives$expvar <- NULL
  }else{stop("\n Method not supported.\n")}


  subannotations<-subset(subannotations, !duplicated(name))

  if(rm.unknown){
    subannotations<-subannotations[!grepl("Unknown|NA", get(level)),]
    mat<-mat[,!grepl("Unknown|NA", colnames(mat))]
  }

  return(list(network=generate.network(mat, covars = covars, annotations = subannotations, correlation.type = correlation.type, alpha = alpha),
              representatives=representatives, level=level))
}

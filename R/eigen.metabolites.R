#' Calculate eigenmetabolites as representatives
#'
#'
#' @param data a \code{\link[data.table]{data.table}} containing the data, where columns correspond to variables and
#' rows to observations.
#' @param group A list containing the grouping (e.g. pathways) of the variables. If this is \code{NULL} all
#' variables will be treated as they were in the same group and only one representation will
#' be calculated
#' @param method the method for calculating eigenmetabolites. Currently PCA and SVD are supported.
#' @references
#' \insertRef{Langfelder2007}{MoDentify}
#' @export
#' @import data.table
#' @importFrom stats prcomp
#' @return list of two, M: \code{\link[data.table]{data.table}} of eigenscores,
#' expvar: full lists of explained variances
#' @examples
#' eigen.data<-eigen.data(data=qmdiab.data,
#'  group = qmdiab.annos$Sub.pathway, method = "PCA")
eigen.metabolites = function(data ,group=NULL, method="PCA") {
  if(is.null(group)){
    group<-rep("Group", dim(data)[2])
  }
  unique.group=unique(group)
  if(method == "PCA"){
    res=sapply(unique.group,function(g){

      pca = prcomp(as.data.table(data[,group==g]))
      list(pc1=pca$x[,1],
           expvar=(pca$sdev)^2 / sum(pca$sdev^2) )
    }, simplify=F)
    # assemble matrix
    M = as.data.table(sapply(res, function(x){x$pc1}))
    # assemble expvars
    expvar = sapply(res, function(x){x$expvar}, simplify=F)
    # return
    return(list(M=M,expvar=expvar))
  }else if(method == "SVD"){
    res=sapply(unique.group,function(g){
      svd = svd(t(as.data.table(data[,group==g])))
      list(v1=svd$v[,1])
    }, simplify=F)
    # assemble matrix
    M = as.data.table(sapply(res, function(x){x$v1}))
    # return
    return(list(M=M))
  }

}

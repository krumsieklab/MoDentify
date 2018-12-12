#' Calculate Module Representative
#'
#'
#' @param data a \code{\link[data.table]{data.table}} with three columns: name, 
#' sampleID and z-score.
#' @param representative.method the method, that is used for calculation of the 
#' representatives. Currently implemented:"eigenmetabolite" and "average"
#'
#' @importFrom stats prcomp
#' @import data.table
#' @references
#' \insertRef{Langfelder2007}{MoDentify}
#' @references
#' \insertRef{Chuang2007}{MoDentify}
#' @return a \code{\link[data.table]{data.table}} containing the module 
#' representative for each sample
calculate.module.representatives<-function(data, representative.method="average"){
  if(representative.method=="average"){
    repdata<-data[,.(representative=mean(z.score)),by=.(sampleID)]
  }else if(representative.method=="eigenmetabolite"){
    if("met.name" %in% colnames(data)){
      data<-dcast(data, sampleID ~ met.name, value.var = "z.score")
    }else{
      data<-dcast(data, sampleID ~ name, value.var = "z.score")
    }
    rownames(data) <- data$sampleID
    data$sampleID <- NULL
    data <- as.data.frame(data)
    eigendata <- eigenMetabolites(data=data, method = "PCA")
    repdata<-data.table(representative=eigendata$M$Group)
    
  }else{stop("\n Method not supported.\n")}
  
  return(repdata)
}

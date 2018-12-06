#' Aggregated Mean Scores
#'
#'
#' @param data a \code{\link[data.table]{data.table}} containing the data, where the columns correspond to the variables and
#' the rows to the samples
#' @param group A list containing the grouping corresponding to the variables. If this is \code{NULL} all
#' variables will be treated as they were in the same group and only one representative will
#' be calculated
#'
#' @export
#' @return matrix of aggregated mean values for each group x sample
#' @examples
#' data(qmdiab.data)
#' data(qmdiab.annos)
#' scaled.data<-scale(qmdiab.data)
#' aggregated.z.scores<-aggregated.mean(data=scaled.data,
#' group = qmdiab.annos$Sub.pathway)
aggregated.mean = function(data, group) {
  unique.group=unique(group)
  res=sapply(unique.group,function(g){apply(as.data.frame(data[,group==g]),1,mean, na.rm=T  )})
  colnames(res)=unique.group
  return(res)
}

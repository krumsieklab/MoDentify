#' Order Annotations
#' 
#' @keywords internal
#'
#'
#' @import data.table
#'
#' @return
orderAnnotation<-function(annotSubnetDT, names){

  newDT<-annotSubnetDT[name==names[1]]
  for(i in 2:length(names)){
    newDT<-rbind(newDT, annotSubnetDT[name==names[i]])
  }
  return(newDT)
}

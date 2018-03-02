#' Order Annotations
#' 
#' @keywords internal
#'
#' @param annot_subnet_DT
#' @param names
#'
#' @import data.table
#'
#' @return
order.annotation<-function(annot_subnet_DT, names){

  newDT<-annot_subnet_DT[name==names[1]]
  for(i in 2:length(names)){
    newDT<-rbind(newDT, annot_subnet_DT[name==names[i]])
  }
  return(newDT)
}

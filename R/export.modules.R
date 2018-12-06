#' Export Modules
#'
#'
#' @param modules a list as returned by \code{\link{identify.modules}}
#' @param file a name for the txt file that should be stored
#' @importFrom utils write.table
#' @export
export.modules <- function(modules, file){
  exportdat<-modules$modules[modules$nodes, on=.(moduleID), nomatch=0]
  exportdat[, module.pval:= exp(-as.numeric(module.score))]
  exportdat[, node.pval:= seed.score]
  exportdat[, module.pval.after.adding.node:= exp(-score.after.adding)]
  exportdat<-exportdat[order(moduleID, order.added)]
  write.table(exportdat, file, sep = "\t", row.names = FALSE)
}
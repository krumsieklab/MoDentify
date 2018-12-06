#' Phenotype-driven module identification
#'
#' @docType package
#' @references
#' \insertRef{Do2017}{MoDentify}
#' 
#' \insertRef{Do2018}{MoDentify}
#' @import Rdpack
"_PACKAGE"
utils::globalVariables(c("node2", "from", "node1", "to", "pval.adjust", "p", "cor",
                         "ppear.adjust", "label", "Labels", ".",  "z.score", "sampleID", "DF",
                         "representative", "Fluid", "moduleID", "nodeID", "module.score", "net.graph",
                         "order.added", "seed.score", "score.after.adding", "is.significant",
                         "node.pval", "i.moduleID", "att1", "att2", "xmin", "xmax", "ymin", "ymax",
                         "count", "principal.component.number", "explained.variance", "score",
                         "key.value", "times.accessed", "name", "value", "s.id", "met.name",
                         "module.beta", "module.pval", "module.pval.after.adding.node",
                         "pearp.adjust", "pval", "adjusted.pval", "module.pval", 
                         "module.pval.after.adding.node"),
                        package = "MoDentify", add=FALSE)


#' qmdiab annotations
#'@docType data
#'@format a \code{\link[data.table]{data.table}} with 1524 rows and 14 variables
#'@usage data(qmdiab.annos)
"qmdiab.annos"

#' qmdiab data
#'@docType data
#'@format a \code{\link[data.table]{data.table}} with 310 rows and 1524 colums
#'@usage data(qmdiab.data)
"qmdiab.data"

#' qmdiab phenotypes
#'@docType data
#'@format a \code{\link[data.table]{data.table}} with 310 rows and 4 columns
#'@usage data(qmdiab.phenos)
"qmdiab.phenos"


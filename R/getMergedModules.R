#' Get Module scores for Merged Modules
#'
#'
#' @param graph an \code{\link[igraph]{igraph}} object, which can be generated
#' with \code{\link{generateNetwork}}. The ID of the nodes must correspond to
#' the name of the variables.
#' @param data either a matrix, where the columns correspond to the variables
#' and the rows to the observations. Or a \code{\link[data.table]{data.table}}
#' with three columns: name, sampleID and value.
#' @param phenotype  a vector with the values for a phenotype of interest.
#' It must have the same number of observations as in data.
#' @param covars a \code{\link[data.table]{data.table}} containing the
#' covariates to correct for.
#' @param nodes an \code{\link[igraph]{igraph}} object, containing the
#' information for the nodes and their new modules.
#' @param scoringFunction a scoring function accepting parameters 
#' moduleRepresentatives, phenotype and covars. See \code{\link[MoDentify]{linearScoring}}
#'
#' @import data.table
#' @import igraph
#'
#' @return
getMergedModules <- function(graph, data, phenotype, covars, nodes,
                             scoringFunction=linearScoring) {
  modules <- nodes[, unique(moduleID)]
  module.scores <- c()
  for (module in modules) {
    tmp <- nodes[moduleID == module, nodeID]
    module.scores <- rbind(module.scores, calculateModuleScore(graph, tmp, data, 
                                                               phenotype, covars, 
                                                               scoringFunction=scoringFunction))
  }

  modules_DT <- data.table(moduleID = modules, module.score = module.scores)
  setnames(modules_DT, c("moduleID", "module.score", "module.beta"))
  return(modules_DT)
}

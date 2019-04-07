#' Calculate Module-score
#'
#'
#' @param graph an igraph object, which can be generated with \code{\link{generateNetwork}}.
#' The ID of the nodes must correspond to the name of the variables.
#' @param nodes a vector containing the ID of the nodes contained in the module
#' @param data a \code{\link[data.table]{data.table}} with three columns: name, sampleID and z-score.
#' @param phenotype a vector with the values for a phenotype of interest.
#' It must have the same number of observations as in data.
#' @param covars a \code{\link[data.table]{data.table}} containing the covariates.
#' The rows for the observations must be in the same order as in the phenotype vector.
#' @param representative.method the method, that is used for the calculation of the eigenmetabolites.
#' Currently implemented: "eigenmetabolite" and "average"
#' @param scoringFunction a scoring function accepting parameters 
#' moduleRepresentatives, phenotype and covars. See \code{\link[MoDentify]{linearScoring}}
#'
#' @import data.table
#'
#' @references
#' \insertRef{Do2017}{MoDentify}
#' @references
#' \insertRef{Chuang2007}{MoDentify}
#' @export
#' @importFrom stats lm
#' @import igraph
#' @return a list containing the module score and the regression coefficient
#' @examples
#' data(qmdiab.data)
#' data(qmdiab.annos)
#' data(qmdiab.phenos)
#' 
#' net.graph <- generateNetwork(data = qmdiab.data, annotations = qmdiab.annos)
#' 
#' data <- data.table(
#'   sampleID = paste0("sample", 1:dim(qmdiab.data)[1]),
#'   qmdiab.data
#' )
#' data <- melt(data = data, id.vars = "sampleID", variable.name = "name")
#' data[, z.score := scale(value), by = .(name)]
#' 
#' module.nodes <- c(3, neighbors(net.graph, 3))
#' 
#' module.score <- calculateModuleScore(
#'   graph = net.graph, nodes = module.nodes,
#'   data = data, phenotype = qmdiab.phenos$T2D
#' )
calculateModuleScore <- function(graph, nodes, data, phenotype, covars = NULL,
                                 representative.method = "average", scoringFunction=linearScoring) {
  # Subsetting the variables for the ones contained in the module
  data <- data[name %in% vertex_attr(graph, "name", nodes)]

  # Calculate the pathway representative
  moduleRepresentatives <- calculateModuleRepresentatives(data, representative.method)
  
  return(scoringFunction(moduleRepresentatives, phenotype, covars))

}

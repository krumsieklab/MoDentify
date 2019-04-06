#' Greedy Module Selection for a Seed
#'
#'
#' @param graph an igraph object, which can be generated with \code{\link{generateNetwork}}.
#' The ID of the nodes must correspond to the name of the variables.
#' @param nodeNr the number of the node, which should be used as a seed.
#' @param data a \code{\link[data.table]{data.table}} with three columns: name,
#' sampleID and z-score.
#' @param phenotype a vector with the values for a phenotype of interest.
#' It must have the same number of samples as in data.
#' @param covars a \code{\link[data.table]{data.table}} containing the covariates to correct for.
#' The rows for the observations must be in the same order as in the phenotype vector.
#' @param alpha significance level (type 1 error) for accepting the modules.
#' @param better.than.components if \code{TRUE}, modules will only be enlarged and accepted,
#' if they are better than all of their components.
#' @param representative.method the method used for the calculation of the module representatives.
#' Currently implemented: "eigenmetabolite" and "average"
#' @param cacheFolder location of the temporary files, where the cache is stored
#'
#' @references
#' \insertRef{Do2017}{MoDentify}
#' @import data.table
#' @import igraph
#' @importFrom digest digest
#' @export
#' @return a list containing the members of the module, the module-score, its regression coefficient
#' for the given phenotype, the score and regression coefficient for the seed, the cache
#' and the consecutive module-scores after adding each new node.
#' @examples
#' data(qmdiab.data)
#' data(qmdiab.annos)
#' data(qmdiab.phenos)
#' 
#' net.graph <- generateNetwork(data = qmdiab.data, annotations = qmdiab.annos)
#' data <- data.table(
#'   sampleID = paste0("sample", 1:dim(qmdiab.data)[1]),
#'   qmdiab.data
#' )
#' data <- melt(data = data, id.vars = "sampleID", variable.name = "name")
#' data[, z.score := scale(value), by = .(name)]
#' 
#' module <- greedyModuleSelection(
#'   graph = net.graph, nodeNr = 51, data = data,
#'   phenotype = qmdiab.phenos$T2D
#' )
greedyModuleSelection <- function(nodeNr, graph, data, phenotype, covars = NULL,
                                  alpha = 0.05, moduleCache = NULL, cacheFolder = NULL, 
                                  better.than.components = TRUE, 
                                  representative.method = "average") {
    
    
    module <- c(nodeNr)
    ## Calculate score for seed
    
    seed <- calculateModuleScore(graph, nodeNr, data, phenotype, covars, 
                                 representative.method = representative.method)
    seed.score <- seed$score
    high.score <- seed.score
    beta <- 0
    old_module <- c()
    score.sequence <- c(seed.score)
    neighbor.index <- 2
    highest.not.signifikant <- seed.score
    
    while (length(old_module) != length(module) & high.score == highest.not.signifikant) {
        old_module <- module
        ## get all neighbors of current module
        neighbors <- getNeighborsForModule(graph, module)
        for (neighbor in neighbors) {
            new_module <- unique(c(old_module, neighbor))
            
            
            neighbor.result <- calculateModuleScore(graph, neighbor, data, phenotype, covars, representative.method = representative.method)
            neighbor.score <- neighbor.result$score
            
            currentDigest <- digest(new_module)
            if(!is.null(cacheFolder) && currentDigest %in% dir(cacheFolder)){
                
                tmp <- as.numeric(read.csv(file = paste0(cacheFolder, currentDigest), 
                                           header = FALSE))
                current.score <- tmp[1]
                current.beta <- tmp[2]
                rm(tmp)
            }else{
                current <- calculateModuleScore(graph, new_module, data, phenotype, covars, representative.method = representative.method)
                current.score <- current$score
                current.beta <- current$beta
                # moduleCache<-rbind(moduleCache, data.table(key.value = asString, 
                #                                            score = current.score, 
                #                                            beta=current$beta, times.accessed=0))
                
                if(!is.null(cacheFolder)){
                    fileName <- paste0(cacheFolder, "/", currentDigest)
                    write(x= c(current.score, current.beta), file = fileName, sep = ",")
                }
            }
            
            
            
            if (better.than.components) {
                if (current.score < high.score & current.score < highest.not.signifikant) {
                    highest.not.signifikant <- current.score
                    if (current.score < neighbor.score) {
                        score.sequence[neighbor.index] <- current.score
                        high.score <- current.score
                        beta <- current.beta
                        module <- new_module
                    }
                }
            } else {
                if (current.score < high.score) {
                    score.sequence[neighbor.index] <- current.score
                    high.score <- current.score
                    highest.not.signifikant <- high.score
                    beta <- current.beta
                    module <- new_module
                }
            }
        }
        neighbor.index <- neighbor.index + 1
    }
    
    
    
    return(list(module = module, module.score = high.score, beta = beta,  
                moduleCache = moduleCache,seed.score = seed.score, 
                seed.beta = seed$beta,  score.sequence = score.sequence))
}
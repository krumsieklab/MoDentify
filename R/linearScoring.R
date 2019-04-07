#' Linear Scoring Function
#' 
#' @description scores the module in regard to a given phenotype using a linear regression
#'
#' @param moduleRepresentatives a \code{\link[data.table]{data.table}} containing the module
#' representative for each sample as returned by \code{\link[MoDentify]{calculateModuleRepresentatives}}
#' @param phenotype a vector with the values for a phenotype of interest.
#' It must have the same number of observations as in data.
#' @param covars a \code{\link[data.table]{data.table}} containing the covariates.
#' The rows for the observations must be in the same order as in the phenotype vector.
#'
#' @return a list containing the module score and the regression coefficient
#' 
#' @seealso \code{\link[stats]{lm}}
#' @export linearScoring
#'
#' @examples
#' #TODO
linearScoring <- function(moduleRepresentatives, phenotype, covars) {
    if (is.null(covars)) {
        DT <- data.table(y = moduleRepresentatives[, representative], phenotype = phenotype)
    } else {
        DT <- data.table(y = moduleRepresentatives[, representative], phenotype = phenotype, covars)
    }
    lm.Module <- lm(formula = y ~ ., data = DT)
    
    sum <- summary(lm.Module)
    p_value <- sum$coefficients[2, "Pr(>|t|)"]
    return(list(score = p_value, beta = lm.Module$coefficients[2]))
}
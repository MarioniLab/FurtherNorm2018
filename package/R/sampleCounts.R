#' Sample counts
#' 
#' Sample count data from a negative binomial (NB) distribution.
#'
#' @param means A matrix of observation-specific means.
#' @param dispersions A vector of gene-specific dispersions, or a matrix of observation-specific dispersions.
#'
#' @details
#' This simulates count data from a NB distribution with parameters estimated from \code{\link{generateRawMeans}}.
#'
#' @seealso
#' \code{\link{generateRawMeans}} to obtain parameters.
#'
#' @return
#' A numeric matrix of sampled counts with the same dimensions as \code{means}.
#'
#' @author
#' Aaron Lun
#'
#' @export
#' @examples
#' out <- generateRawMeans()
#' counts <- sampleCounts(out$fitted, out$dispersions)
sampleCounts <- function(means, dispersions) {
    matrix(rnbinom(length(means), mu=means, size=1/dispersions), ncol=ncol(means))
}


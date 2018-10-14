#' Spawn simulation scenarios
#'
#' Creates multiple listings of simulation parameter combinations.
#'
#' @param ... Named arguments. Each argument is usually a vector of possible values for a single simulation parameter.
#'
#' @details
#' The output of this function is intended to be passed to a \code{\link{bpmapply}} call.
#' This allows easy parallel execution across multiple simulation scenarios.
#'
#' @return
#' A named list of vectors, where each vector corresponds to a parameter.
#' When processed in parallel, the vectors contain all possible combinations of parameter settings provided in \code{...}.
#'
#' @author
#' Aaron Lun
#' 
#' @export
#' @examples
#' spawnScenarios(ncells=c(100, 200, 500), effect=c(2, 5))
spawnScenarios <- function(...) 
{
    components <- list(...)
    N <- prod(lengths(components))

    output <- components
    niter <- 1L
    per.iter <- N
    for (i in seq_along(components)) {
        current <- components[[i]]
        curN <- length(current)
        output[[i]] <- rep(rep(current, each=per.iter/curN), niter)
        niter <- niter * curN
        per.iter <- per.iter / curN
    }

    return(output)
}



#' Cluster-specific composition bias
#'
#' Compute the composition bias that remains after normalization, based on comparisons between clusters.
#'
#' @param counts Numeric matrix, count data for genes (in rows) and cells/samples (in columns).
#' @param sf Numeric vector, size factors for all cells.
#' @param clust Factor, cluster identity for each cell.
#' @param threshold Numeric scalar, threshold on the average count for bias estimation.
#'
#' @details
#' We apply \code{sf} to \code{counts} to obtain normalized counts for all cells,
#' and then compute the average normalized count for each level of \code{clust}.
#' If cell-specific composition biases were successfully removed by \code{sf}, there should not be any composition biases between clusters.
#'
#' We use \code{\link{calcNormFactors}} to estimate the remaining bias between clusters.
#' This is possible as the cluster averages should have few zeroes \emph{and} be very precise (thus reducing the bias of the trimmed mean when DE is imbalanced).
#' Ideally, the returned values should be equal to unity, indicating that there are no biases in the cluster-specific averages.
#' 
#' Note that \code{lib.size} is set to a constant (of 1 million) to avoid re-introducing composition biases in \code{\link{calcNormFactors}}.
#'
#' @return 
#' Normalization factors representing the scaling bias between clusters.
#'
#' @author
#' Aaron Lun
#'
#' @export
#' @importFrom scater normalizeCounts calcAverage
#' @importFrom Matrix rowMeans
#' @importFrom edgeR calcNormFactors
#'
#' @examples
#' counts <- matrix(rpois(100000, lambda=2), ncol=100)
#' sf <- scater::librarySizeFactors(counts)
#' clust <- sample(3, ncol(counts), replace=TRUE)
#' clusterBias(counts, sf, clust)
clusterBias <- function(counts, sf, clust, threshold=0.1) {
    out <- normalizeCounts(counts, size_factors=sf, return_log=FALSE)

    U <- unique(clust)
    collected <- vector("list", length(U))
    names(collected) <- U
    for (x in U) {
        chosen <- clust==x
        collected[[x]] <- rowMeans(out[,chosen,drop=FALSE])
    }
    by.clust <- do.call(cbind, collected)

    keep <- calcAverage(by.clust) > threshold
    calcNormFactors(by.clust[keep,], lib.size=rep(mean(colSums(by.clust)), ncol(by.clust)))
}


#' Cluster-specific composition bias
#'
#' Compute the composition bias that remains after normalization, based on comparisons between clusters.
#'
#' @param counts Numeric matrix, count data for genes (in rows) and cells/samples (in columns).
#' @param sf Numeric vector, size factors for all cells.
#' @param clust Factor, cluster identity for each cell.
#' @param groups Character vector of cluster identities for which to compute bias.
#' Clusters with low numbers of cells should be ignored.
#' @param threshold Numeric scalar, threshold on the average count for bias estimation.
#'
#' @details
#' We apply \code{sf} to \code{counts} to obtain normalized counts for all cells,
#' and then compute the average normalized count for each level of \code{clust}.
#' If cell-specific composition biases were successfully removed by \code{sf}, there should not be any composition biases between clusters.
#'
#' We use a median-based approach to estimate the remaining bias between clusters relative to a single reference.
#' This combines the pairwise approach of TMM normalization (with weaker assumptions) and the relative robustness of \pkg{DESeq} normalization.
#' It is possible as the cluster averages should have few zeroes \emph{and} be very precise (thus reducing the bias of the trimmed mean when DE is imbalanced).
#' 
#' Ideally, the returned values should be equal to unity, indicating that there are no biases in the cluster-specific averages.
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
#'
#' @examples
#' counts <- matrix(rpois(100000, lambda=2), ncol=100)
#' sf <- scater::librarySizeFactors(counts)
#' clust <- sample(3, ncol(counts), replace=TRUE)
#' clusterBias(counts, sf, clust)
clusterBias <- function(counts, sf, clust, groups=NULL, threshold=0, ref=1) {
    clust <- as.factor(clust)
    U <- levels(clust)
    if (!is.null(groups)) {
        keep <- clust %in% groups
        counts <- counts[,keep,drop=FALSE]
        sf <- sf[keep]
        clust <- clust[keep]
        U <- intersect(U, groups)
    }

    sf <- sf/mean(sf)
    out <- normalizeCounts(counts, size_factors=sf, return_log=FALSE)

    by.clust <- vector("list", length(U))
    names(by.clust) <- U
    for (x in U) {
        chosen <- clust==x
        by.clust[[x]] <- rowMeans(out[,chosen,drop=FALSE])
    }

    sf <- numeric(length(by.clust))
    names(sf) <- names(by.clust)
    ref.clust <- by.clust[[ref]]
    for (x in U) {
        current <- by.clust[[x]]
        keep <- calcAverage(cbind(current, ref.clust)) > threshold
        sf[x] <- median(current[keep]/ref.clust[keep], na.rm=TRUE)
    }
    sf
}


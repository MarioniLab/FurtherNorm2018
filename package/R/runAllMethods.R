#' Run all methods
#'
#' Compute size factors with a variety of normalization methods.
#'
#' @param counts Numeric matrix, count data with genes in rows and cells/samples in columns.
#' @param threshold Numeric scalar, threshold on the mean below which low-abundance genes are to be removed.
#'
#' @details
#' We apply a variety of methods:
#' \itemize{
#' \item TMM normalization with \code{\link{calcNormFactors}}.
#' \item TMM normalization using an averaged expression profile as the reference sample.
#' \item DESeq normalization with \code{\link{estimateSizeFactorsForMatrix}}, where all zeroes are ignored during calculation of the geometric mean.
#' \item DESeq normalization using the arithmetic mean.
#' \item DESeq normalization with an added pseudo-count of 1.
#' \item DESeq normalization with a library size-adjusted pseudo-count 
#' \item Library size normalization with \code{\link{librarySizeFactors}}.
#' \item Deconvolution with \code{\link{computeSumFactors}} and without any pre-clustering.
#' \item Deconvolution with \code{\link{computeSumFactors}} and pre-clustering with \code{\link{quickCluster}}.
#' \item The nearest-neighbor approach in \code{\link{simpleSumFactors}}.
#' }
#' For the TMM and DESeq approaches, filtering on \code{threshold} is performed based on the average abundances computed from \code{\link{calcAverage}}.
#' For the \pkg{scran}-based methods, filtering is done by passing \code{threshold} to the \code{min.mean} arguments.
#' Library size normalization is unaffected by \code{threshold}.
#'
#' @return
#' A named list containing numeric vectors of size factors computed from each method.
#'
#' @author Aaron Lun
#'
#' @export
#' @importFrom edgeR calcNormFactors
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @importFrom scran quickCluster computeSumFactors simpleSumFactors
#' @importFrom scater calcAverage librarySizeFactors
#'
#' @examples
#' counts <- matrix(rpois(100000, lambda=2), ncol=100)
#' out <- runAllMethods(counts)
#' str(out)
runAllMethods <- function(counts, threshold=0) {
{
    subcounts <- counts[calcAverage(counts) > threshold,]
    subcounts <- as.matrix(subcounts)

    # TMM with raw counts:
    tmm.sf <- calcNormFactors(subcounts) * colSums(subcounts)

    # TMM with averaged counts:
    combined <- cbind(rowSums(subcounts), subcounts)
    tmm2.sf <- calcNormFactors(combined, refColumn=1) * colSums(combined)
    tmm2.sf <- tmm2.sf[-1]

    # Size factors with raw counts (must counter zeros for both geometric mean and for each library):
    logvals <- log(subcounts)
    logvals[is.infinite(logvals)] <- NA_real_
    gm <- exp(rowMeans(logvals, na.rm=TRUE))
    size.sf <- estimateSizeFactorsForMatrix(subcounts, geoMeans=gm)

    # Size factors with averaged counts (still removes zeros in each library):
    size2.sf <- estimateSizeFactorsForMatrix(subcounts, geoMeans=rowMeans(subcounts))

    # Size factors with an added pseudo-count.
    sizeP.sf <- estimateSizeFactorsForMatrix(subcounts+1)
                                
    # Size factors with a library size-adjusted pseudo-count.
    lib.size <- colSums(subcounts)
    pcounts <- t(t(subcounts) + lib.size/mean(lib.size))
    sizeP2.sf <- estimateSizeFactorsForMatrix(pcounts)
}

# Library size normalization.
lib.sf <- librarySizeFactors(counts)

# Size factors with summation, no abundance filtering for simplicity.
sizes <- seq(21, 101, by=5)
sizes <- sizes[sizes <= ncol(counts)]
final.sf <- computeSumFactors(counts, sizes=sizes, clusters=NULL, min.mean=threshold)

# Size factors with clustering prior to summation:
if (ncol(counts) >= 200) {
    emp.clusters <- quickCluster(counts, method="igraph")
    final2.sf <- computeSumFactors(counts, sizes=sizes, clusters=emp.clusters, min.mean=threshold)
} else {
    final2.sf <- final.sf
}

# Using a knn-based method for summation.
final3.sf <- simpleSumFactors(counts, min.mean=threshold, approximate=TRUE)

    # Reporting all methods.
    return(list(TMM=tmm.sf, TMM.ave=tmm2.sf,
                DESeq.geo=size.sf, DESeq.ave=size2.sf, DESeq.pseudo=sizeP.sf, DESeq.pseudo.lib=sizeP2.sf,
                Lib=lib.sf,
                Deconv=final.sf, Deconv.clust=final2.sf, QuickSum=final3.sf))
}


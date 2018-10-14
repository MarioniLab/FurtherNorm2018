#' Test for DE
#'
#' Use \pkg{edgeR} to test for differential expression between two groups.
#'
#' @param counts Numeric matrix, count data for genes (in rows) and cells/samples (in columns).
#' @param sf Numeric vector, size factors for all cells.
#' @param group1 Integer vector, identities of cells in the first group.
#' @param group2 Integer vector, identities of cells in the second group.
#' @param lfc Numeric scalar, threshold on the log-fold change for \code{\link{glmTreat}}.
#'
#' @details
#' This function uses the TREAT methods in \pkg{edgeR} to test for differential expression between groups.
#' The idea is to identify changes in the detected DEGs when the size factors change, due to more-or-less effective removal of composition biases between groups.
#'
#' @return
#' A \code{\link{DGELRT}} object with the results of the comparison between groups.
#' Note that log-fold changes refer to the second group over the first group.
#'
#' @author Aaron Lun
#'
#' @export
#' @importFrom stats model.matrix
#' @importFrom edgeR DGEList aveLogCPM estimateDisp glmFit glmTreat 
#'
#' @examples
#' counts <- matrix(rpois(100000, lambda=2), ncol=100)
#' sf <- scater::librarySizeFactors(counts)
#' testDE(counts, sf, 1:50, 51:100)
testDE <- function(counts, sf, group1, group2, lfc=1) {
    to.use <- as.matrix(counts[,c(group1, group2)])
    sf.use <- sf[c(group1, group2)]
    nf <- sf.use/colSums(to.use)
    nf <- nf/mean(nf)
    y <- DGEList(to.use, norm.factors=nf)

    g <- factor(rep(1:2, c(length(group1), length(group2))))
    design <- model.matrix(~g)
    y <- estimateDisp(y, design)
    fit <- glmFit(y, design)
    glmTreat(fit, lfc=lfc)
}


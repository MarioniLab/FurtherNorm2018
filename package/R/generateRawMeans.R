#' Generate raw means
#'
#' Compute the mean expression of each gene across a population of cells, assuming no DE genes.
#'
#' @param ngenes Integer, number of genes.
#' @param ncells Integer, number of cells.
#' @param mode String, type of scRNA-seq data.
#'
#' @details
#' All counts are presumed to be sampled from a negative binomial (NB) distribution.
#' The mean of the distribution for any observation is the product of the gene-specific mean and the cell-specific size factor.
#' The size factor for each cell is sampled from a log-normal distribution.
#' 
#' If \code{mode="UMI"}, we use low gene-specific means and low dispersions.
#' This mimics what is observed in most UMI data sets, with parameters taken from eyeballing plots of the original inDrop data set.
#'
#' If \code{mode="read"}, we use high gene-specific means and large dispersions with a decreasing mean-dispersion trend.
#' This mimics what is observed in many read-based data sets, with parameters taken from eyeballing plots of the ESpresso data set.
#'
#' @return
#' A list containing \code{fitted}, a numeric matrix of means for each observation;
#' \code{mean}, a numeric vector of the gene-specific means;
#' \code{dispersion}, a numeric vector of the gene-specific dispersions;
#' and \code{sf}, a numeric vector of the true cell-specific size factors.
#'
#' @author Aaron Lun
#' @export
#' @importFrom stats rgamma runif rnorm
#'
#' @examples
#' out <- generateRawMeans()
#' str(out)
generateRawMeans <- function(ngenes=10000, ncells=200, mode=c("UMI", "read")) {
    mode <- match.arg(mode)
    if (mode=="UMI") {
        # Guestimated from inDrop data.
        true.means <- rgamma(ngenes, 2, 2)
        dispersions <- rep(0.1, ngenes)
    } else {
        # Guestimated from ESpresso data (see PlateEffect supplementaries).
        true.means <- 2^runif(ngenes, 0, 12)
        dispersions <- 100/true.means^0.6
    }

    all.facs <- 2^rnorm(ncells, sd=1)
    return(list(fitted=outer(true.means, all.facs),
                mean=true.means, dispersion=dispersions, sf=all.facs))
}
